//! Источник отсчётов с АЦП (режим adc) — для РЕАЛЬНОГО железа.
//!
//! Каждые 4 мс (250 Гц) снимает один отсчёт с аналогового входа PA1.
//! В Renode не запускается - плата
//!
//!
//! ВНИМАНИЕ: этот модуль скомпилируется на первом железном билде; API АЦП/таймера
//! HAL 0.23 использован по документации, но именно этот путь я не смог проверить
//! кросс-сборкой — если cargo ругнётся, правки будут здесь (тип пина/тайминга).
//! Тут я хз че делать только реальную плату

use embedded_hal::delay::DelayNs;
use stm32f4xx_hal::{
    adc::{
        config::{AdcConfig, SampleTime},
        Adc,
    },
    gpio::{Analog, Pin},
    pac::ADC1,
    prelude::*,
    rcc::Rcc,
};

/// Период дискретизации: 1_000_000 мкс/с / 250 Гц = 4000 мкс.
const PERIOD_US: u32 = 4000;

/// Обобщён по типу задержки (`D: DelayNs`), чтобы не завязываться на точный
/// тип `Delay<TIM2, ...>` — он выводится в main.
pub struct AdcSource<D: DelayNs> {
    adc: Adc<ADC1>,
    pin: Pin<'A', 1, Analog>,
    delay: D,
}

impl<D: DelayNs> AdcSource<D> {
    pub fn new<M>(adc1: ADC1, pa1: Pin<'A', 1, M>, delay: D, rcc: &mut Rcc) -> Self {
        let adc = Adc::new(adc1, true, AdcConfig::default(), rcc);
        let pin = pa1.into_analog();
        AdcSource { adc, pin, delay }
    }

    /// Ждёт следующий тик 250 Гц и возвращает свежий отсчёт. Никогда не None.
    pub fn next(&mut self) -> Option<i16> {
        self.delay.delay_us(PERIOD_US);
        // 12-битный код 0..4095 -> i16. DC/масштаб не важны: полосовой фильтр
        // ядра убирает постоянную составляющую, пороги адаптивные.
        let raw: u16 = self.adc.convert(&self.pin, SampleTime::Cycles_480);
        Some(raw as i16)
    }
}
