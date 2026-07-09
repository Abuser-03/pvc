import sys
import numpy as np
import os
from ZoomView import ZoomViewer
from tables import table2, table3, table4, table5

path = os.getcwd()

def check_path(filename, dp):
    f_path = os.path.join(path, filename + dp)
    return os.path.exists(f_path)
  
if __name__ == '__main__':
    if len(sys.argv) > 1:
        fname = sys.argv[1]
    else:    
        #fname = input("Введите название файла: ")
        fname = 'mf1_II'
    print( fname )    

    viewer = ZoomViewer( 200, [-1000, 1000])
    x, y, dummy = table3( fname + '.v1')
    x = x - 3922
    ys = y - np.mean(y)
    viewer.AddLine(x, ys, "v1")  
    #viewer.AddMarkers( x_plus, rr_plus, 's', "plus", 'magenta' ) 
    viewer.Draw(fname)