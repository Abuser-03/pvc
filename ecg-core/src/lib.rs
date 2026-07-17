
//! Интерфейс потоковый: `Analyzer::push(sample)` принимает ОДИН отсчёт и,
//! когда удар подтверждён, возвращает отчёт по нему. «Банчи» — это забота
//! вызывающего кода (просто цикл по чанку); состояние живёт между вызовами,
//! поэтому граница банча на результат не влияет.

#![no_std]


// Параметры под частоту дискретизации 250 Гц (1 отсчёт = 4 мс).
// Всё, что зависит от fs, посчитано здесь в ОТСЧЁТАХ.
pub const FS: u32 = 250;

const MWI_WIN: usize = 38; // окно интегрирования ~150 мс (0.15 * 250)
const REFRACTORY: u32 = 45; // рефрактерный период ~180 мс (беты при 300 уд/мин идут через 50 отсчётов)
const TWAVE_WIN: u32 = 90; // окно дискриминации T-волны ~360 мс
const RLOC_WIN: usize = 48; // окно поиска истинного R по полосовому сигналу (~192 мс)
const LEARN_SAMPLES: u32 = 500; // фаза обучения порогов ~2 с

// Пороги ритма (в отсчётах RR):
const BRADY_RR: i32 = 250; // RR > 1000 мс  => ЧСС < 60  => брадикардия
const TACHY_RR: i32 = 150; // RR <  600 мс  => ЧСС > 100 => тахикардия

// Преждевременность: удар считаем PVC, если RR < PVC_PCT% от опорного RR.
const PVC_PCT: i32 = 85;

// Бигеминия: столько строгих чередований типов (норма/PVC) подряд => бигеминия.
const BIGEM_MIN_ALT: u32 = 6;

// Кольцевой буфер фиксированного размера. Без деления (МК зачастую без HW-divide):
// перенос индекса делаем сравнением, а не `%`.

struct Ring<T: Copy + Default, const N: usize> {
    buf: [T; N],
    pos: usize, // индекс последнего записанного элемента
}

impl<T: Copy + Default, const N: usize> Ring<T, N> {
    fn new() -> Self {
        Ring { buf: [T::default(); N], pos: 0 }
    }
    #[inline]
    fn push(&mut self, v: T) {
        self.pos += 1;
        if self.pos >= N {
            self.pos = 0;
        }
        self.buf[self.pos] = v;
    }
    /// `d` отсчётов назад: get(0) — только что записанный, get(1) — предыдущий, ...
    #[inline]
    fn get(&self, d: usize) -> T {
        let mut i = self.pos as isize - d as isize;
        if i < 0 {
            i += N as isize;
        }
        self.buf[i as usize]
    }
}

// Ступени препроцессинга Pan–Tompkins (все целочисленные).
//   вход -> ФНЧ -> ФВЧ(=полоса 5..15 Гц) -> производная -> квадрат -> интегрирование
// Возвращает на каждый отсчёт (mwi, bandpass) — интеграл и полосовой сигнал.

struct Preprocessor {
    // ФНЧ: y[n] = 2y[n-1] - y[n-2] + x[n] - 2x[n-6] + x[n-12]
    lp_x: Ring<i32, 13>,
    lp_y1: i32,
    lp_y2: i32,
    // ФВЧ: ma[n]=ma[n-1]+lp[n]-lp[n-32];  hp[n]=lp[n-16]-ma[n]/32
    hp_x: Ring<i32, 33>,
    hp_ma: i32,
    // Производная: d = (2h[n] + h[n-1] - h[n-3] - 2h[n-4]) / 8
    dv_x: Ring<i32, 5>,
    // Интегрирование: скользящая сумма квадратов
    mwi_sq: Ring<i32, MWI_WIN>,
    mwi_sum: i64, // i64 только у аккумулятора — гарантированно без переполнения
}

impl Preprocessor {
    fn new() -> Self {
        Preprocessor {
            lp_x: Ring::new(),
            lp_y1: 0,
            lp_y2: 0,
            hp_x: Ring::new(),
            hp_ma: 0,
            dv_x: Ring::new(),
            mwi_sq: Ring::new(),
            mwi_sum: 0,
        }
    }

    #[inline]
    fn step(&mut self, x: i16) -> (i32, i32) {
        let x = x as i32;

        // --- ФНЧ ---
        self.lp_x.push(x);
        let lp = 2 * self.lp_y1 - self.lp_y2 + self.lp_x.get(0)
            - 2 * self.lp_x.get(6)
            + self.lp_x.get(12);
        self.lp_y2 = self.lp_y1;
        self.lp_y1 = lp;

        // --- ФВЧ (полосовой на выходе) ---
        self.hp_x.push(lp);
        self.hp_ma += self.hp_x.get(0) - self.hp_x.get(32);
        let bandpass = self.hp_x.get(16) - (self.hp_ma >> 5); // /32

        // --- производная ---
        self.dv_x.push(bandpass);
        let mut d = (2 * self.dv_x.get(0) + self.dv_x.get(1)
            - self.dv_x.get(3)
            - 2 * self.dv_x.get(4))
            >> 3; // /8

        // клип производной в диапазон i16, чтобы квадрат гарантированно влез в i32
        if d > 32767 {
            d = 32767;
        } else if d < -32768 {
            d = -32768;
        }

        // --- квадрат ---
        let sq = d * d; // <= 32768^2 < 2^31, влезает в i32

        // --- интегрирование (скользящее среднее квадратов) ---
        self.mwi_sq.push(sq);
        self.mwi_sum += sq as i64 - self.mwi_sq.get(MWI_WIN - 1) as i64;
        let mwi = (self.mwi_sum / MWI_WIN as i64) as i32;

        (mwi, bandpass)
    }
}

// Детектор QRS: адаптивные пороги + рефрактерность + дискриминация T-волны +
// локализация истинного R по полосовому сигналу.
// push() возвращает Some(индекс_R), когда удар подтверждён.

struct QrsDetector {
    pre: Preprocessor,
    bp_ring: Ring<i32, RLOC_WIN>, // полосовой сигнал для поиска пика R
    n: u32,                       // счётчик отсчётов (u32 — не usize!)

    // детекция локального максимума MWI
    mwi_prev: i32,
    rising: bool,

    // адаптивные пороги (целые)
    spki: i32, // оценка уровня "сигнала"
    npki: i32, // оценка уровня "шума"
    learning: bool,
    learn_max: i32,
    learn_sum: i64,

    last_qrs_n: u32,   // индекс последнего принятого R
    last_qrs_peak: i32, // амплитуда MWI последнего R (для T-волны)
    have_qrs: bool,
    rr_avg: i32, // среднее RR (все принятые беты) — для поиска-назад

    // Кандидат для поиска-назад: самый крупный подпороговый пик (между thr2 и thr1),
    // встреченный с момента последнего R. Если R давно не было — продвигаем его.
    sb_peak_val: i32,
    sb_peak_n: u32,
    sb_r_index: u32,
}

impl QrsDetector {
    fn new() -> Self {
        QrsDetector {
            pre: Preprocessor::new(),
            bp_ring: Ring::new(),
            n: 0,
            mwi_prev: 0,
            rising: false,
            spki: 0,
            npki: 0,
            learning: true,
            learn_max: 0,
            learn_sum: 0,
            last_qrs_n: 0,
            last_qrs_peak: 0,
            have_qrs: false,
            rr_avg: 0,
            sb_peak_val: 0,
            sb_peak_n: 0,
            sb_r_index: 0,
        }
    }

    /// Принять пик как QRS: обновить пороги/RR, вернуть индекс R.
    fn accept(&mut self, peak_val: i32, peak_n: u32, r_index: u32) -> u32 {
        self.spki = self.spki - (self.spki >> 3) + (peak_val >> 3); // 0.875*SPKI+0.125*peak
        if self.have_qrs {
            let rr = peak_n.wrapping_sub(self.last_qrs_n) as i32;
            if self.rr_avg == 0 {
                self.rr_avg = rr;
            } else {
                let band = self.rr_avg >> 2;
                if rr >= self.rr_avg - band && rr <= self.rr_avg + band {
                    self.rr_avg += (rr - self.rr_avg) >> 3;
                }
            }
        }
        self.last_qrs_n = peak_n;
        self.last_qrs_peak = peak_val;
        self.have_qrs = true;
        self.sb_peak_val = 0; // сбрасываем кандидата поиска-назад
        r_index
    }

    #[inline]
    fn threshold_i1(&self) -> i32 {
        // Классический Pan–Tompkins: NPKI + 0.25*(SPKI-NPKI)  [>>2].
        // Здесь 0.125 [>>3]: сигнал с симулятора чистый, а широкие PVC дают
        // низкий MWI — более низкий порог их достаёт; T-волны при этом
        // отсекаются отдельной дискриминацией (на 300 уд/мин ложняка нет).
        // Для шумного реального ЭКГ верни >>2.
        self.npki + ((self.spki - self.npki) >> 3)
    }

    /// Найти индекс истинного R: максимум |полосового| в окне RLOC_WIN.
    fn locate_r(&self) -> u32 {
        let mut best_d = 0usize;
        let mut best_v = -1i32;
        let mut d = 0usize;
        while d < RLOC_WIN {
            let mut v = self.bp_ring.get(d);
            if v < 0 {
                v = -v;
            }
            if v > best_v {
                best_v = v;
                best_d = d;
            }
            d += 1;
        }
        self.n - best_d as u32
    }

    #[inline]
    fn push(&mut self, x: i16) -> Option<u32> {
        let (mwi, bandpass) = self.pre.step(x);
        self.bp_ring.push(bandpass);
        self.n += 1;

        // --- фаза обучения порогов (~2 с) ---
        if self.learning {
            if mwi > self.learn_max {
                self.learn_max = mwi;
            }
            self.learn_sum += mwi as i64;
            self.mwi_prev = mwi;
            if self.n >= LEARN_SAMPLES {
                self.learning = false;
                self.spki = self.learn_max / 4; // ~0.25 от пикового
                let mean = (self.learn_sum / self.n as i64) as i32;
                self.npki = mean; // средний уровень как оценка шума
            }
            return None;
        }

        // --- детекция локального максимума MWI ---
        let mut peak_val = 0i32;
        let mut peak_n = 0u32;
        let mut have_peak = false;
        if mwi > self.mwi_prev {
            self.rising = true;
        } else if mwi < self.mwi_prev && self.rising {
            // предыдущий отсчёт был локальным пиком
            peak_val = self.mwi_prev;
            peak_n = self.n - 1;
            have_peak = true;
            self.rising = false;
        }
        self.mwi_prev = mwi;

        let thr1 = self.threshold_i1();
        let thr2 = thr1 >> 1;

        // --- обработка найденного пика ---
        if have_peak {
            let refractory_ok =
                !self.have_qrs || peak_n.wrapping_sub(self.last_qrs_n) >= REFRACTORY;
            if refractory_ok {
                // дискриминация T-волны: близкий и низкий пик = T-волна
                let is_twave = self.have_qrs
                    && peak_n.wrapping_sub(self.last_qrs_n) < TWAVE_WIN
                    && peak_val < (self.last_qrs_peak >> 1);

                if is_twave {
                    self.npki = self.npki - (self.npki >> 3) + (peak_val >> 3);
                } else if peak_val > thr1 {
                    // ---- QRS принят сразу ----
                    let r = self.locate_r();
                    return Some(self.accept(peak_val, peak_n, r));
                } else {
                    // ниже основного порога: подкручиваем оценку шума,
                    // но если пик выше половинного порога — запоминаем как
                    // кандидата для поиска-назад (возможный пропущенный удар).
                    self.npki = self.npki - (self.npki >> 3) + (peak_val >> 3);
                    if peak_val > thr2 && peak_val > self.sb_peak_val {
                        self.sb_peak_val = peak_val;
                        self.sb_peak_n = peak_n;
                        self.sb_r_index = self.locate_r();
                    }
                }
            }
        }

        // --- ретроспективный поиск-назад ---
        // Если R не было дольше ~1.6*RR — продвигаем лучшего подпорогового кандидата.
        if self.have_qrs && self.rr_avg > 0 {
            let missed = self.rr_avg + (self.rr_avg >> 1) + (self.rr_avg >> 3); // ~1.625*rr_avg
            if self.n.wrapping_sub(self.last_qrs_n) > missed as u32 && self.sb_peak_val > 0 {
                return Some(self.accept(self.sb_peak_val, self.sb_peak_n, self.sb_r_index));
            }
        }

        None
    }
}

// Классификатор ритма: по последовательности R-пиков (RR-интервалам) ставит
// метку каждому удару (норма/PVC) и выдаёт эпизоды ритма.


/// Тип удара.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum BeatKind {
    Normal,
    /// Желудочковая экстрасистола (преждевременный удар).
    Pvc,
}

/// Тип эпизода ритма, «сработавшего» на данном ударе.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum RhythmKind {
    Bradycardia,
    Tachycardia,
    /// Бигеминия: чередование норма-PVC-норма-PVC.
    Bigeminy,
    /// Куплет: 2 PVC подряд.
    Couplet,
    /// Пробежка/неустойчивая ЖТ: >= 3 PVC подряд.
    Run,
}

/// Отчёт по подтверждённому удару.
#[derive(Clone, Copy, Debug)]
pub struct BeatReport {
    /// Индекс отсчёта R-пика (0-based от начала потока).
    pub sample: u32,
    /// RR-интервал, ведущий в этот удар, в отсчётах (0 — для первого удара).
    pub rr_samples: u32,
    /// ЧСС по этому RR, уд/мин (0 — для первого удара).
    pub bpm: u32,
    pub kind: BeatKind,
    /// Эпизод ритма, если он начался/подтвердился на этом ударе.
    pub rhythm: Option<RhythmKind>,
}

struct Classifier {
    have_prev: bool,
    prev_r: u32,
    rr_ref: i32, // опорный RR "нормального" ритма (EMA по нормальным ударам)

    // бигеминия: считаем строгое чередование ТИПОВ ударов (норма/PVC)
    have_prev_kind: bool,
    prev_kind: BeatKind,
    kind_alt: u32,
    bigem_active: bool,

    pvc_run: u32,
    rate_state: i32, // -1 бради, 0 норма, +1 тахи
}

impl Classifier {
    fn new() -> Self {
        Classifier {
            have_prev: false,
            prev_r: 0,
            rr_ref: 0,
            have_prev_kind: false,
            prev_kind: BeatKind::Normal,
            kind_alt: 0,
            bigem_active: false,
            pvc_run: 0,
            rate_state: 0,
        }
    }

    fn on_beat(&mut self, r: u32) -> BeatReport {
        if !self.have_prev {
            self.have_prev = true;
            self.prev_r = r;
            return BeatReport {
                sample: r,
                rr_samples: 0,
                bpm: 0,
                kind: BeatKind::Normal,
                rhythm: None,
            };
        }

        let rr = r.wrapping_sub(self.prev_r) as i32;
        self.prev_r = r;
        let bpm = if rr > 0 { (15000 / rr) as u32 } else { 0 }; // 60000/(rr*4)

        if self.rr_ref == 0 {
            self.rr_ref = rr; // затравка
        }

        // --- метка удара по преждевременности ---
        let premature = rr * 100 < self.rr_ref * PVC_PCT;
        let kind = if premature { BeatKind::Pvc } else { BeatKind::Normal };

        // Опорный RR обновляем только нормальными ударами И только если RR
        // в разумной полосе вокруг текущего (0.75..1.25). Иначе пропущенный
        // удар (RR ~= 2x) утащил бы опорный вверх и всё залипло бы.
        if kind == BeatKind::Normal {
            let band = self.rr_ref >> 2;
            if rr >= self.rr_ref - band && rr <= self.rr_ref + band {
                self.rr_ref += (rr - self.rr_ref) >> 3;
            }
        }

        // --- бигеминия: строгое чередование норма-PVC-норма-PVC по типам ---
        if !self.have_prev_kind {
            self.have_prev_kind = true;
            self.kind_alt = 0;
        } else if kind != self.prev_kind {
            self.kind_alt += 1; // тип сменился — чередование продолжается
        } else {
            self.kind_alt = 0; // два одинаковых подряд — не бигеминия
        }
        self.prev_kind = kind;
        let entering_bigem = !self.bigem_active && self.kind_alt >= BIGEM_MIN_ALT;
        if entering_bigem {
            self.bigem_active = true;
        }
        if self.kind_alt == 0 {
            self.bigem_active = false;
        }

        // счётчик подряд идущих PVC
        if kind == BeatKind::Pvc {
            self.pvc_run += 1;
        } else {
            self.pvc_run = 0;
        }

        // --- выбор эпизода ритма (по приоритету) ---
        let mut rhythm = None;
        if self.pvc_run == 3 {
            rhythm = Some(RhythmKind::Run);
        } else if self.pvc_run == 2 {
            rhythm = Some(RhythmKind::Couplet);
        } else if entering_bigem {
            rhythm = Some(RhythmKind::Bigeminy);
        } else {
            // бради/тахи по опорному (нормальному) ритму, с гистерезисом
            if self.rr_ref > BRADY_RR && self.rate_state != -1 {
                self.rate_state = -1;
                rhythm = Some(RhythmKind::Bradycardia);
            } else if self.rr_ref < TACHY_RR && self.rate_state != 1 {
                self.rate_state = 1;
                rhythm = Some(RhythmKind::Tachycardia);
            } else if self.rr_ref <= BRADY_RR && self.rr_ref >= TACHY_RR {
                self.rate_state = 0;
            }
        }

        BeatReport { sample: r, rr_samples: rr as u32, bpm, kind, rhythm }
    }
}


// Верхнеуровневый анализатор: отсчёт -> (детектор) -> (классификатор) -> отчёт.

pub struct Analyzer {
    det: QrsDetector,
    cls: Classifier,
}

impl Analyzer {
    pub fn new() -> Self {
        Analyzer { det: QrsDetector::new(), cls: Classifier::new() }
    }

    /// Скормить один отсчёт ЭКГ. Возвращает Some(отчёт), если удар подтверждён.
    #[inline]
    pub fn push(&mut self, sample: i16) -> Option<BeatReport> {
        let r = self.det.push(sample)?;
        Some(self.cls.on_beat(r))
    }
}

impl Default for Analyzer {
    fn default() -> Self {
        Self::new()
    }
}
