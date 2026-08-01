//! Потоковый детектор QRS (Pan–Tompkins) + классификатор аритмий.
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
// Не знал как показать но что-то вроде такого должно быть под ревью должно уйти !!!! #TODO
const BRADY_RR: i32 = 250; // RR > 1000 мс  => ЧСС < 60  => брадикардия
const TACHY_RR: i32 = 150; // RR <  600 мс  => ЧСС > 100 => тахикардия

// Сколько последних RR держим для опорного значения. Опорный RR берём как
// МЕДИАНУ этой истории, а не скользящее среднее: при тригеминии треть
// интервалов короткие, и среднее уползает за ними, а медиана — нет.
const RR_HIST: usize = 8;

// Ширина комплекса (отсчёты), начиная с которой считаем его ШИРОКИМ.
// Измерено на записях: наджелудочковые 7..9, желудочковые 14 и выше;
// 12 проходит посередине с запасом в обе стороны.
const WIDE_QRS: u32 = 12;

// Устойчивая ЖТ: столько подряд широких И быстрых ударов.
// Пробежки из 5 экстрасистол дают максимум 4 подряд, поэтому 6 их не задевает.
const VT_MIN_RUN: u32 = 6;

// Фибрилляция: из последних 16 ударов столько должны быть "рваными" по RR,
// И столько же широкими, И столько же быстрыми. Три условия сразу нужны,
// потому что одного хаоса мало: у реальной записи здорового человека с
// артефактами хаос доходит до 10 из 16 — но там ритм нормальной частоты и
// комплексы узкие, а при фибрилляции всё три признака совпадают.
const FIB_MIN: u32 = 8;

// Для «широких» при фибрилляции порог мягче: комплексы там рваные и половина
// их не дотягивает до WIDE_QRS. Запас всё равно огромный — у реальной записи
// здорового человека с артефактами широких максимум 1 из 16, здесь нужно 5.
const FIB_MIN_WIDE: u32 = 5;

// Преждевременность: удар считаем PVC, если RR < PVC_PCT% от опорного RR. Нейронка от нее краевая задача
const PVC_PCT: i32 = 85;

// Бигеминия: столько строгих чередований типов (норма/PVC) подряд => бигеминия. Так же
const BIGEM_MIN_ALT: u32 = 6;

// Кольцевой буфер фиксированного размера. Без деления (МК зачастую без HW-divide):
// перенос индекса делаем сравнением, а не `%`.

struct Ring<T: Copy + Default, const N: usize> {
    buf: [T; N],
    pos: usize, // инкс последнего в памяти элемента
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
    // ФНЧ: y[n] = 2y[n-1] - y[n-2] + x[n] - 2x[n-6] + x[n-12] TODO
    lp_x: Ring<i32, 13>,
    lp_y1: i32,
    lp_y2: i32,
    // ФВЧ: ma[n]=ma[n-1]+lp[n]-lp[n-32];  hp[n]=lp[n-16]-ma[n]/32 TODO
    hp_x: Ring<i32, 33>,
    hp_ma: i32,
    // Производная: d = (2h[n] + h[n-1] - h[n-3] - 2h[n-4]) / 8 TODO
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

    /// Возвращает (интеграл, полосовой, производная).
    /// Производная нужна снаружи как мера КРУТИЗНЫ фронта: у QRS она резкая,
    /// у зубца T пологая — по одной амплитуде их не различить.
    #[inline]
    fn step(&mut self, x: i16) -> (i32, i32, i32) {
        let x = x as i32;

        //  ФНЧ
        self.lp_x.push(x);
        let lp = 2 * self.lp_y1 - self.lp_y2 + self.lp_x.get(0)
            - 2 * self.lp_x.get(6)
            + self.lp_x.get(12);
        self.lp_y2 = self.lp_y1;
        self.lp_y1 = lp;

        //  ФВЧ (полосовой на выходе)
        self.hp_x.push(lp);
        self.hp_ma += self.hp_x.get(0) - self.hp_x.get(32);
        let bandpass = self.hp_x.get(16) - (self.hp_ma >> 5); // /32

        //  производная
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


        let sq = d * d;

        //  интегрирование (скользящее среднее квадратов)
        self.mwi_sq.push(sq);
        self.mwi_sum += sq as i64 - self.mwi_sq.get(MWI_WIN - 1) as i64;
        let mwi = (self.mwi_sum / MWI_WIN as i64) as i32;

        (mwi, bandpass, d)
    }
}

// Детектор QRS: адаптивные пороги + рефрактерность + дискриминация T-волны +
// локализация истинного R по полосовому сигналу.
// push() возвращает Some(индекс_R), когда удар подтверждён.
// Должно вроде работать TODO
// Дискриминация T-волны по КРУТИЗНЕ: пик внутри окна T-волны считаем зубцом T,
// если его фронт положе TWAVE_SLOPE_PCT% от фронта предыдущего комплекса.
// Одной амплитуды мало: при подъёме сегмента ST зубец T набирает высоту и
// проходит амплитудный критерий, а вот резким он не становится.
// Порог выбран по измерениям: у ложных детекций отношение 61..64%,
// у настоящих быстрых ударов (синусовая 300/мин) ~100%.
const TWAVE_SLOPE_PCT: i32 = 70;

// Артефакт: пик интегратора, во столько раз превышающий текущую оценку
// сигнала, физиологическим быть не может (отрыв электрода, движение).
// Такой пик не даёт удара И НЕ ОБНОВЛЯЕТ пороги — иначе адаптивный порог
// взлетает и детектор глохнет на десятки секунд.
const ARTIFACT_MULT: i32 = 16;

/// Групповая задержка полосового фильтра в отсчётах: ФНЧ даёт 6, ФВЧ 16.
/// Пик, найденный по полосовому сигналу, отстоит от истинного R ровно на эту
/// величину — без компенсации все позиции R уезжали на 88 мс вперёд.
/// На RR-интервалы это не влияло (сдвиг постоянный и сокращается), но ломало
/// привязку к форме сигнала: измерение ширины и будущий R-on-T.
const GROUP_DELAY: u32 = 22;

/// Модуль целого (без ветвлений в горячем пути хватает и такого).
#[inline]
fn iabs(v: i32) -> i32 {
    if v < 0 {
        -v
    } else {
        v
    }
}

struct QrsDetector {
    pre: Preprocessor,
    bp_ring: Ring<i32, RLOC_WIN>, // полосовой сигнал для поиска пика R
    dv_ring: Ring<i32, RLOC_WIN>, // |производная| — мера крутизны фронта
    n: u32,                       // счётчик отсчётов (u32 — не usize!)

    // морфология последнего выданного удара (заполняется в emit)
    last_width: u32, // ширина комплекса в отсчётах
    last_slope: i32, // максимальная крутизна фронта
    artifacts: u32,  // счётчик отбракованных артефактов (диагностика качества)
    twaves: u32,     // счётчик отбракованных зубцов T (диагностика)

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
    last_r: u32,        // индекс ПОСЛЕДНЕГО ВЫДАННОГО R (не пика MWI)
    have_qrs: bool,
    rr_avg: i32, // среднее RR (все принятые беты) — для поиска-назад

    // Кандидат для поиска-назад: самый крупный подпороговый пик (между thr2 и thr1),
    // встреченный с момента последнего R. Если R давно не было — продвигаем его.
    sb_peak_val: i32,
    sb_peak_n: u32,
    sb_r_index: u32,
    sb_slope: i32, // крутизна кандидата — чтобы проверить его на зубец T при продвижении
}

impl QrsDetector {
    fn new() -> Self {
        QrsDetector {
            pre: Preprocessor::new(),
            bp_ring: Ring::new(),
            dv_ring: Ring::new(),
            last_width: 0,
            last_slope: 0,
            artifacts: 0,
            twaves: 0,
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
            last_r: 0,
            have_qrs: false,
            rr_avg: 0,
            sb_peak_val: 0,
            sb_peak_n: 0,
            sb_r_index: 0,
            sb_slope: 0,
        }
    }

    /// Выдать удар, только если найденный R достаточно далеко от предыдущего R.
    /// Закрывает сразу две дыры:
    ///   * два разных пика интегратора могут указать на ОДИН И ТОТ ЖЕ R
    ///     (окно поиска R шире рефрактерного периода) — был двойной счёт;
    ///   * ветка поиска-назад раньше выдавала удар вообще без проверки
    ///     рефрактерности — отсюда брались интервалы по 44 мс.
    fn emit(&mut self, peak_val: i32, peak_n: u32, r: u32) -> Option<u32> {
        if self.have_qrs && (r <= self.last_r || r - self.last_r < REFRACTORY) {
            return None;
        }
        self.last_width = self.measure_width();
        self.last_slope = self.measure_slope();
        Some(self.accept(peak_val, peak_n, r))
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
        self.last_r = r_index;
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
        // Я не ебу то как это обрабатывать так что тут вот такая обвязка
        self.npki + ((self.spki - self.npki) >> 3)
    }

    /// Ширина комплекса в отсчётах — по длительности всплеска ПРОИЗВОДНОЙ.
    /// Считаем именно по ней, а не по полосовому сигналу: полосовой фильтр
    /// «размазывает» узкие комплексы собственной импульсной характеристикой,
    /// и узкий наджелудочковый комплекс становится неотличим от широкого
    /// желудочкового. Производная же не зависит ни от изолинии, ни от усиления.
    fn measure_width(&self) -> u32 {
        let mut best_d = 0usize;
        let mut best_v = -1i32;
        let mut d = 0usize;
        while d < RLOC_WIN {
            let v = self.dv_ring.get(d);
            if v > best_v {
                best_v = v;
                best_d = d;
            }
            d += 1;
        }
        if best_v <= 0 {
            return 0;
        }
        let level = best_v >> 2; // четверть от пика крутизны

        let mut lo = best_d;
        while lo > 0 && self.dv_ring.get(lo - 1) > level {
            lo -= 1;
        }
        let mut hi = best_d;
        while hi + 1 < RLOC_WIN && self.dv_ring.get(hi + 1) > level {
            hi += 1;
        }
        (hi - lo + 1) as u32
    }

    /// Максимальная крутизна фронта в окне (макс |производная|).
    fn measure_slope(&self) -> i32 {
        let mut best = 0i32;
        let mut d = 0usize;
        while d < RLOC_WIN {
            let v = self.dv_ring.get(d);
            if v > best {
                best = v;
            }
            d += 1;
        }
        best
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
        self.n - best_d as u32 - GROUP_DELAY
    }

    #[inline]
    fn push(&mut self, x: i16) -> Option<u32> {
        let (mwi, bandpass, deriv) = self.pre.step(x);
        self.bp_ring.push(bandpass);
        self.dv_ring.push(iabs(deriv));
        self.n += 1;

        //  фаза обучения порогов (+-2 с)
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

        //  детекция локального максимума MWI
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

        //  обработка найденного пика
        if have_peak {
            let refractory_ok =
                !self.have_qrs || peak_n.wrapping_sub(self.last_qrs_n) >= REFRACTORY;
            if refractory_ok {
                let slope_now = self.measure_slope();

                // Артефакт: пик кратно выше всего, что мы видели как сигнал.
                // Ни удара, ни обновления порогов — просто пропускаем.
                let artifact = self.have_qrs
                    && self.spki > 0
                    && peak_val / ARTIFACT_MULT > self.spki;

                // Позицию R считаем заранее: окно T-волны меряем между
                // САМИМИ R, а не между пиками интегратора. У пологих волн
                // (спад приподнятого ST) интегратор запаздывает сильнее, и по
                // его пикам такой всплеск «выпадал» за окно и проскакивал.
                let r_now = self.locate_r();
                let since_r = r_now.wrapping_sub(self.last_r);

                // дискриминация T-волны: близкий пик, который либо низкий,
                // либо пологий (см. TWAVE_SLOPE_PCT)
                let is_twave = self.have_qrs
                    && since_r < TWAVE_WIN
                    && (peak_val < (self.last_qrs_peak >> 1)
                    || slope_now * 100 < self.last_slope * TWAVE_SLOPE_PCT);

                if artifact {
                    self.artifacts = self.artifacts.wrapping_add(1);
                } else if is_twave {
                    self.twaves = self.twaves.wrapping_add(1);
                    self.npki = self.npki - (self.npki >> 3) + (peak_val >> 3);
                } else if peak_val > thr1 {
                    //  QRS принят сразу
                    if let Some(r) = self.emit(peak_val, peak_n, r_now) {
                        return Some(r);
                    }
                } else {
                    // ниже основного порога: подкручиваем оценку шума,
                    // но если пик выше половинного порога — запоминаем как
                    // кандидата для поиска-назад (возможный пропущенный удар).
                    self.npki = self.npki - (self.npki >> 3) + (peak_val >> 3);
                    if peak_val > thr2 && peak_val > self.sb_peak_val {
                        self.sb_peak_val = peak_val;
                        self.sb_peak_n = peak_n;
                        self.sb_r_index = r_now;
                        self.sb_slope = slope_now;
                    }
                }
            }
        }

        //  ретроспективный поиск-назад
        // Если R не было дольше ~1.6*RR — продвигаем лучшего подпорогового кандидата.
        if self.have_qrs && self.rr_avg > 0 {
            let missed = self.rr_avg + (self.rr_avg >> 1) + (self.rr_avg >> 3); // ~1.625*rr_avg
            if self.n.wrapping_sub(self.last_qrs_n) > missed as u32 && self.sb_peak_val > 0 {
                let (v, pn, ri, sl) =
                    (self.sb_peak_val, self.sb_peak_n, self.sb_r_index, self.sb_slope);
                self.sb_peak_val = 0; // кандидат израсходован в любом случае

                // Кандидат тоже обязан пройти проверку на зубец T. Без неё
                // поиск-назад возвращал обратно ровно те пики, которые прямой
                // путь только что отбраковал как T-волну.
                let cand_twave = ri.wrapping_sub(self.last_r) < TWAVE_WIN
                    && (v < (self.last_qrs_peak >> 1)
                    || sl * 100 < self.last_slope * TWAVE_SLOPE_PCT);
                if cand_twave {
                    self.twaves = self.twaves.wrapping_add(1);
                } else if let Some(r) = self.emit(v, pn, ri) {
                    return Some(r);
                }
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
    /// Устойчивая желудочковая тахикардия: быстрый ритм ШИРОКИМИ комплексами.
    /// Отличается от синусовой тахикардии именно шириной: при синусовой
    /// комплекс узкий, преждевременности нет ни там, ни там.
    VentricularTachycardia,
    /// Фибрилляция желудочков: хаотичный + быстрый + широкий.
    /// Требует немедленной дефибрилляции — высший приоритет тревоги.
    Fibrillation,
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
    /// Ширина комплекса в отсчётах (мера "широкий/узкий QRS").
    pub width_samples: u32,
    /// Максимальная крутизна фронта (мера "резкий/пологий").
    pub slope: i32,
    /// Эпизод ритма, если он начался/подтвердился на этом ударе.
    pub rhythm: Option<RhythmKind>,
}

struct Classifier {
    have_prev: bool,
    prev_r: u32,
    rr_ref: i32, // опорный RR "нормального" ритма (медиана истории RR)
    rr_hist: [i32; RR_HIST], // последние RR (порядок не важен — берём медиану)
    rr_cnt: usize,           // сколько ячеек уже заполнено (насыщается)
    rr_pos: usize,           // куда писать следующий

    // Битовые истории последних 16 ударов (1 бит на удар) — дёшево и без деления.
    chaos_bits: u16, // RR отличается от предыдущего более чем на 50%
    wide_bits: u16,  // комплекс широкий
    fast_bits: u16,  // интервал короткий (быстрый ритм)
    prev_rr: i32,

    vt_run: u32,       // подряд идущих широких+быстрых ударов
    vt_active: bool,   // эпизод ЖТ уже объявлен
    fib_active: bool,  // эпизод фибрилляции уже объявлен

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
            rr_hist: [0; RR_HIST],
            rr_cnt: 0,
            rr_pos: 0,
            chaos_bits: 0,
            wide_bits: 0,
            fast_bits: 0,
            prev_rr: 0,
            vt_run: 0,
            vt_active: false,
            fib_active: false,
            have_prev_kind: false,
            prev_kind: BeatKind::Normal,
            kind_alt: 0,
            bigem_active: false,
            pvc_run: 0,
            rate_state: 0,
        }
    }

    /// Записать очередной RR в кольцевую историю.
    fn rr_push(&mut self, rr: i32) {
        self.rr_hist[self.rr_pos] = rr;
        self.rr_pos += 1;
        if self.rr_pos >= RR_HIST {
            self.rr_pos = 0;
        }
        if self.rr_cnt < RR_HIST {
            self.rr_cnt += 1;
        }
    }

    /// Медиана истории RR. Сортировка вставками по копии на стеке:
    /// не больше 8 элементов, без кучи и без деления.
    fn rr_median(&self) -> i32 {
        let n = self.rr_cnt;
        if n == 0 {
            return 0;
        }
        let mut tmp = [0i32; RR_HIST];
        let mut i = 0;
        while i < n {
            tmp[i] = self.rr_hist[i];
            i += 1;
        }
        let mut i = 1;
        while i < n {
            let v = tmp[i];
            let mut j = i;
            while j > 0 && tmp[j - 1] > v {
                tmp[j] = tmp[j - 1];
                j -= 1;
            }
            tmp[j] = v;
            i += 1;
        }
        tmp[n / 2]
    }

    fn on_beat(&mut self, r: u32, width: u32, slope: i32) -> BeatReport {
        if !self.have_prev {
            self.have_prev = true;
            self.prev_r = r;
            return BeatReport {
                sample: r,
                rr_samples: 0,
                bpm: 0,
                kind: BeatKind::Normal,
                width_samples: width,
                slope,
                rhythm: None,
            };
        }

        let rr = r.wrapping_sub(self.prev_r) as i32;
        self.prev_r = r;
        let bpm = if rr > 0 { (15000 / rr) as u32 } else { 0 }; // 60000/(rr*4)

        if self.rr_ref == 0 {
            self.rr_ref = rr; // затравка
        }

        //  метка удара по преждевременности (сравниваем с опорным ДО учёта
        //  текущего интервала, иначе экстрасистола сама себя оправдает)
        let premature = rr * 100 < self.rr_ref * PVC_PCT;
        let kind = if premature { BeatKind::Pvc } else { BeatKind::Normal };

        // Опорный RR = медиана последних RR_HIST интервалов. В историю кладём
        // ВСЕ интервалы, включая экстрасистолические: медиана к ним устойчива,
        // пока их меньше половины. Прежний вариант (скользящее среднее только
        // по "нормальным" в полосе 0.75..1.25) при тригеминии затравливался
        // коротким интервалом и залипал — экстрасистолы переставали находиться.
        self.rr_push(rr);
        self.rr_ref = self.rr_median();

        // Битовые истории по последним 16 ударам: рваность RR, ширина, темп.
        let jump = self.prev_rr > 0 && {
            let d = if rr > self.prev_rr { rr - self.prev_rr } else { self.prev_rr - rr };
            d * 2 > self.prev_rr // разница более 50%
        };
        self.prev_rr = rr;
        self.chaos_bits = (self.chaos_bits << 1) | (jump as u16);
        self.wide_bits = (self.wide_bits << 1) | ((width >= WIDE_QRS) as u16);
        self.fast_bits = (self.fast_bits << 1) | ((rr < TACHY_RR) as u16);

        // Устойчивая ЖТ: непрерывная серия широких И быстрых ударов.
        if width >= WIDE_QRS && rr < TACHY_RR {
            self.vt_run += 1;
        } else {
            self.vt_run = 0;
            self.vt_active = false;
        }

        //  бигеминия: строгое чередование норма-PVC-норма-PVC по типам
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

        //  выбор эпизода ритма (по приоритету)
        // Фибрилляция и устойчивая ЖТ идут ВЫШЕ всего остального: это
        // состояния, требующие немедленных действий, и объявлять на них
        // "куплет" или "тахикардию" было бы опасной недооценкой.
        let fib_now = self.chaos_bits.count_ones() >= FIB_MIN
            && self.wide_bits.count_ones() >= FIB_MIN_WIDE
            && self.fast_bits.count_ones() >= FIB_MIN;
        if !fib_now {
            self.fib_active = false;
        }

        let mut rhythm = None;
        if fib_now && !self.fib_active {
            self.fib_active = true;
            rhythm = Some(RhythmKind::Fibrillation);
        } else if fib_now {
            // эпизод уже объявлен — молчим, чтобы не заваливать тревогами
        } else if self.vt_run >= VT_MIN_RUN && !self.vt_active {
            self.vt_active = true;
            rhythm = Some(RhythmKind::VentricularTachycardia);
        } else if self.pvc_run == 3 {
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

        BeatReport {
            sample: r,
            rr_samples: rr as u32,
            bpm,
            kind,
            width_samples: width,
            slope,
            rhythm,
        }
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
        let (w, sl) = (self.det.last_width, self.det.last_slope);
        Some(self.cls.on_beat(r, w, sl))
    }
}

impl Analyzer {
    /// Диагностика качества сигнала: (отбраковано артефактов, отбраковано зубцов T).
    pub fn quality(&self) -> (u32, u32) {
        (self.det.artifacts, self.det.twaves)
    }

}

impl Default for Analyzer {
    fn default() -> Self {
        Self::new()
    }
}