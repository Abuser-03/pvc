#!/usr/bin/env python3
# Превращает запись .v1 в Rust-файл с зашитым массивом отсчётов для прошивки.
#
#   python3 gen_data.py run5.v1            -> ecg-fw/src/ecg_data.rs (по умолчанию)
#   python3 gen_data.py mf1_II.v1 out.rs   -> произвольный путь
#
# Второй столбец .v1 = амплитуда ЭКГ (см. readme.txt). Кладём её как i16.

import sys
import os

def main():
    if len(sys.argv) < 2:
        print("использование: python3 gen_data.py <файл.v1> [выходной.rs]")
        sys.exit(1)
    src = sys.argv[1]
    dst = sys.argv[2] if len(sys.argv) > 2 else os.path.join(
        os.path.dirname(__file__), "ecg-fw", "src", "ecg_data.rs")

    vals = []
    with open(src) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            parts = line.split(",")
            if len(parts) < 2:
                continue
            try:
                y = int(round(float(parts[1])))
            except ValueError:
                continue
            # клип в диапазон i16 на всякий случай
            y = max(-32768, min(32767, y))
            vals.append(y)

    name = os.path.basename(src)
    with open(dst, "w") as f:
        f.write(f"// АВТОГЕНЕРАЦИЯ из {name} ({len(vals)} отсчётов, fs=250 Гц).\n")
        f.write("// Перегенерировать: python3 gen_data.py <файл.v1>\n")
        f.write("#![allow(clippy::all)]\n\n")
        f.write(f"pub static SAMPLES: &[i16] = &[\n")
        # по 16 значений в строке
        for i in range(0, len(vals), 16):
            chunk = ", ".join(str(v) for v in vals[i:i+16])
            f.write("    " + chunk + ",\n")
        f.write("];\n")

    print(f"записано {len(vals)} отсчётов -> {dst}")

if __name__ == "__main__":
    main()
