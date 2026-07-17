//! Прошивка детектора аритмий для STM32F4 (по умолчанию — F407 Discovery).
//!
//! Ядро (`ecg`) не меняется — прошивка лишь: даёт ему отсчёты, а на события
//! зажигает "тревогу" (светодиод) и пишет лог в USART2.
//!
//! Источник отсчётов выбирается фичей на этапе компиляции:
//!   --features test-data  зашитая запись .v1 (для Renode / без железа)
//!   --features adc         реальный сигнал с АЦП (только на живой плате)
//!
//! Вывод: USART2 (TX = PA2). В Renode виден через `showAnalyzer sysbus.usart2`.
//! Такты одного push меряются счётчиком DWT (в Renode — приблизительно,
//! точные — на реальном чипе).

#![no_std]
#![no_main]

use panic_halt as _;

use core::fmt::Write;
use cortex_m::peripheral::DWT;
use cortex_m_rt::entry;
use embedded_hal::digital::OutputPin;

use stm32f4xx_hal::{self as hal, rcc::Config};
use crate::hal::{pac, prelude::*};

use ecg::{Analyzer, BeatKind, RhythmKind};

#[cfg(feature = "test-data")]
mod ecg_data; // src/ecg_data.rs — зашитая запись

#[cfg(feature = "test-data")]
mod source_test;
#[cfg(feature = "adc")]
mod source_adc;

/// Счётчики для итоговой сводки.
struct Stats {
    n: u32,     // отсчётов обработано
    cyc: u64,   // суммарно тактов на push (для среднего)
    beats: u32, // ударов найдено
    pvc: u32,   // из них PVC
}

fn rhythm_name(r: RhythmKind) -> &'static str {
    match r {
        RhythmKind::Bradycardia => "БРАДИКАРДИЯ",
        RhythmKind::Tachycardia => "ТАХИКАРДИЯ",
        RhythmKind::Bigeminy => "БИГЕМИНИЯ",
        RhythmKind::Couplet => "КУПЛЕТ (2 PVC)",
        RhythmKind::Run => "ПРОБЕЖКА ЖТ (>=3 PVC)",
    }
}

/// Обработать один отсчёт: замерить такты push, залогировать событие,
/// зажечь тревогу на опасный паттерн. Общая для обоих источников.
fn process_sample<W: Write, P: OutputPin>(
    ana: &mut Analyzer,
    s: i16,
    tx: &mut W,
    alarm: &mut P,
    st: &mut Stats,
) {
    let t0 = DWT::cycle_count();
    let ev = ana.push(s);
    let dt = DWT::cycle_count().wrapping_sub(t0);
    st.cyc = st.cyc.wrapping_add(dt as u64);
    st.n = st.n.wrapping_add(1);

    if let Some(rep) = ev {
        st.beats += 1;
        if rep.kind == BeatKind::Pvc {
            st.pvc += 1;
        }
        if let Some(r) = rep.rhythm {
            // t в мс = номер отсчёта * 4 (fs = 250 Гц -> 4 мс/отсчёт)
            let _ = writeln!(
                tx,
                "t={:>6}ms  R@{:<6}  {}\r",
                rep.sample.wrapping_mul(4),
                rep.sample,
                rhythm_name(r)
            );
            // Тревога — на желудочковые пробежки/куплеты (опасное для реанимации).
            if matches!(r, RhythmKind::Run | RhythmKind::Couplet) {
                let _ = alarm.set_high();
            }
        }
    }
}

#[entry]
fn main() -> ! {
    let dp = pac::Peripherals::take().unwrap();
    let mut cp = cortex_m::peripheral::Peripherals::take().unwrap();

    // Тактовая: HSI по умолчанию (надёжно в Renode; на железе можно поднять
    // частоту через .sysclk(...), на такты/отсчёт это не влияет).
    let mut rcc = dp.RCC.freeze(Config::hsi());

    // Светодиод "тревога": на F407 Discovery красный на PD14.
    let gpiod = dp.GPIOD.split(&mut rcc);
    let mut alarm = gpiod.pd14.into_push_pull_output();

    // GPIOA: PA2 = USART2 TX (лог), PA1 = вход АЦП (в режиме adc).
    let gpioa = dp.GPIOA.split(&mut rcc);

    // USART2 TX (PA2) — лог. В Renode: showAnalyzer sysbus.usart2.
    let mut tx = dp.USART2.tx(gpioa.pa2, 115200.bps(), &mut rcc).unwrap();

    // Счётчик тактов ядра.
    cp.DCB.enable_trace();
    cp.DWT.enable_cycle_counter();

    let mut ana = Analyzer::new();
    let mut st = Stats { n: 0, cyc: 0, beats: 0, pvc: 0 };

    let _ = writeln!(tx, "=== ECG detector: старт (fs=250) ===\r");

    // -------- источник отсчётов --------
    #[cfg(feature = "test-data")]
    let mut src = source_test::TestSource::new();

    #[cfg(feature = "adc")]
    let mut src = {
        // таймер для отсчёта периода 250 Гц (тип Delay выводится, не именуем его)
        let delay = dp.TIM2.delay_us(&mut rcc);
        source_adc::AdcSource::new(dp.ADC1, gpioa.pa1, delay, &mut rcc)
    };

    // -------- основной цикл --------
    loop {
        match src.next() {
            Some(s) => process_sample(&mut ana, s, &mut tx, &mut alarm, &mut st),
            None => break, // конец записи (бывает только в режиме test-data)
        }
    }

    // Сводка (достигается только когда запись кончилась — режим test-data).
    let avg = if st.n > 0 { st.cyc / st.n as u64 } else { 0 };
    let _ = writeln!(
        tx,
        "=== ГОТОВО: {} отсчётов, {} ударов ({} PVC), ~{} тактов/push ===\r",
        st.n, st.beats, st.pvc, avg
    );

    loop {
        cortex_m::asm::wfi();
    }
}
