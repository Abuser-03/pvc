//! Источник отсчётов из зашитой записи .v1 (режим test-data).
//! Отдаёт отсчёты по одному, пока запись не кончится (тогда None).

use crate::ecg_data::SAMPLES;

pub struct TestSource {
    i: usize,
}

impl TestSource {
    pub fn new() -> Self {
        TestSource { i: 0 }
    }

    pub fn next(&mut self) -> Option<i16> {
        let v = *SAMPLES.get(self.i)?;
        self.i += 1;
        Some(v)
    }
}
