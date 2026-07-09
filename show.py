#!/usr/bin/env python3
# Просмотр результата детектора: рисует сигнал .v1 через ваш ZoomViewer и
# накладывает найденные удары из файла .markers (норма/PVC разным цветом).
#
# Должен лежать рядом с ZoomView.py, tables.py и файлами <name>.v1 / <name>.markers.
# Файл маркеров генерирует Rust-детектор:  ./ecg <name>.v1
#
#   python3 show.py <name>      (например: python3 show.py run5)

import sys
import os
import numpy as np
from ZoomView import ZoomViewer
from tables import table3

if __name__ == '__main__':
    fname = sys.argv[1] if len(sys.argv) > 1 else 'mf1_II'
    fname = fname[:-3] if fname.endswith('.v1') else fname
    print(fname)

    # --- сигнал ---
    x, y, _ = table3(fname + '.v1')
    x0 = x[0]                      # выравниваем по началу (аналог "x - 3922" в view.py)
    mean_y = np.mean(y)           # как в view.py: убираем постоянную составляющую
    xs = x - x0
    ys = y - mean_y

    viewer = ZoomViewer(500, [-1000, 1000])
    viewer.AddLine(xs, ys, "v1")

    # --- маркеры найденных ударов ---
    mpath = fname + '.markers'
    if os.path.exists(mpath):
        mx, my, mk = table3(mpath)          # x_raw, y_raw, kind (0=норма,1=PVC)
        mk = np.asarray(mk).astype(int)
        mxs = mx - x0
        mys = my - mean_y

        norm = mk == 0
        pvc = mk == 1
        if np.any(norm):
            viewer.AddMarkers(mxs[norm], mys[norm], 'o', "норма", 'limegreen')
        if np.any(pvc):
            viewer.AddMarkers(mxs[pvc], mys[pvc], 'v', "PVC", 'red')
        print(f"ударов: {len(mk)}  (норма {int(norm.sum())}, PVC {int(pvc.sum())})")
    else:
        print(f"нет файла {mpath} — сначала запусти детектор:  ./ecg {fname}.v1")

    viewer.Draw(fname)
