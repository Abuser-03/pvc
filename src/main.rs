//! Хост-обвязка (std) вокруг no_std-ядра `ecg`.
//!
//! На микроконтроллере отсчёты капали бы с АЦП; здесь мы ИМИТИРУЕМ поток,
//! читая файл .v1 и скармливая его ядру чанками (банчами). Граница банча на
//! результат не влияет — состояние живёт в Analyzer между чанками.
//!
//! Проверяется ОДИН файл!!!!!!!!

use std::io::Write;
use std::path::Path;

use ecg::{Analyzer, BeatKind, BeatReport, RhythmKind, FS};

const INPUT: &str = "nsr30.v1"; // норм. ритм 30 уд/мин (брадикардия)
// const INPUT: &str = "nsr300.v1"; // норм. ритм 300 уд/мин (тахикардия)
// const INPUT: &str = "mf1_II.v1"; // мультифокусные PVC
// const INPUT: &str = "bigem_II.v1"; // бигеминия
// const INPUT: &str = "run5.v1"; // пробежка из 5 PVC

// Размер банча (имитация потоковой подачи). На результат не влияет.
const BATCH: usize = 64;

fn main() {
    let path = std::env::args().nth(1).unwrap_or_else(|| INPUT.to_string());

    //  читаем файл: строки "индекс, амплитуда, статус"
    let text = match std::fs::read_to_string(&path) {
        Ok(t) => t,
        Err(e) => {
            eprintln!("не могу открыть {path}: {e}");
            std::process::exit(1);
        }
    };

    let mut xs: Vec<i32> = Vec::new(); // столбец 0: индекс отсчёта из файла
    let mut ys: Vec<i16> = Vec::new(); // столбец 1: амплитуда ЭКГ
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let mut it = line.split(',');
        let x = it.next().and_then(|s| s.trim().parse::<i32>().ok());
        let y = it.next().and_then(|s| s.trim().parse::<f64>().ok()); // на всякий случай допускаем дробное
        if let (Some(x), Some(y)) = (x, y) {
            xs.push(x);
            ys.push(y as i16);
        }
    }

    if ys.is_empty() {
        eprintln!("в {path} не нашлось данных");
        std::process::exit(1);
    }

    let n = ys.len();
    let dur_s = n as f64 / FS as f64;
    println!("файл: {path}");
    println!("отсчётов: {n}  (~{dur_s:.1} с при {FS} Гц)\n");

    //  гоним поток через ядро ЧАНКАМИ
    let mut ana = Analyzer::new();
    let mut beats: Vec<BeatReport> = Vec::new();

    let mut i = 0usize;
    while i < n {
        let end = (i + BATCH).min(n);
        for &s in &ys[i..end] {
            if let Some(rep) = ana.push(s) {
                beats.push(rep);
            }
        }
        i = end;
    }

    //  сводка
    let total = beats.len();
    let pvc = beats.iter().filter(|b| b.kind == BeatKind::Pvc).count();
    let normal = total - pvc;

    // средняя ЧСС по нормальным ударам (с ненулевым RR)
    let hr_vals: Vec<u32> = beats
        .iter()
        .filter(|b| b.kind == BeatKind::Normal && b.bpm > 0)
        .map(|b| b.bpm)
        .collect();
    let mean_hr = if hr_vals.is_empty() {
        0
    } else {
        hr_vals.iter().sum::<u32>() / hr_vals.len() as u32
    };

    println!("-*10 итог -*15");
    println!("ударов найдено : {total}");
    println!("  норма        : {normal}");
    println!("  PVC          : {pvc}");
    println!("средняя ЧСС    : {mean_hr} уд/мин (по нормальным)");

    // эпизоды ритма
    let episodes: Vec<&BeatReport> = beats.iter().filter(|b| b.rhythm.is_some()).collect();
    if episodes.is_empty() {
        println!("эпизоды ритма  : ");
    } else {
        println!("эпизоды ритма  :");
        for b in &episodes {
            let t = b.sample as f64 / FS as f64;
            let name = match b.rhythm.unwrap() {
                RhythmKind::Bradycardia => "брадикардия",
                RhythmKind::Tachycardia => "тахикардия",
                RhythmKind::Bigeminy => "бигеминия",
                RhythmKind::Couplet => "куплет (2 PVC)",
                RhythmKind::Run => "пробежка ЖТ (>=3 PVC)",
            };
            println!("  t={t:6.2}с (отсчёт {:6})  {name}", b.sample);
        }
    }

    //  пишем файл маркеров для питоновского вьюера

    // формат строки: "x_raw, y_raw, kind"  (kind: 0=норма, 1=PVC)
    let stem = Path::new(&path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("out");
    let markers_path = format!("{stem}.markers");
    match std::fs::File::create(&markers_path) {
        Ok(mut f) => {
            for b in &beats {
                let idx = b.sample as usize;
                if idx < n {
                    let k = if b.kind == BeatKind::Pvc { 1 } else { 0 };
                    let _ = writeln!(f, "{}, {}, {}", xs[idx], ys[idx], k);
                }
            }
            println!("\nмаркеры записаны: {markers_path}");
            println!("посмотреть:  python3 show.py {stem}");
        }
        Err(e) => eprintln!("не смог записать {markers_path}: {e}"),
    }
}
