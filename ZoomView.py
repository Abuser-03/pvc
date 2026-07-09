# отображение набора кривых и маркеров
# с возможностью выбора участка просмотра и его увеличением

import sys
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
#from scipy.misc import electrocardiogram
#from scipy.signal import find_peaks
    
class ZoomViewer:
    def __init__(self, page=500, yrange=[-1000,1000]):
        self.fig = plt.figure()
        self.ax = self.fig.add_subplot(111)
        self.zoom = 1.0
        self.xstart = 0        # нижний предел по x
        #self.ymiddle = 0.5*(yrange[0]+yrange[1])
        self.ReInit(page, yrange)
        self.x = []
        self.xx = []
        self.yy = []
        self.labels = []
        #self.colors = []
        self.legends = []
        self.nlines = 0
        self.xm = []
        self.ym = []
        self.ms = []
        self.mlabels = []
        self.mcolors = []
        self.nmarkers = 0
        
    def ReInit(self, page, yrange):
        self.page = page
        self.yrange = yrange
        self.ymiddle = 0.5*(yrange[0]+yrange[1])
        self.xstop = self.xstart + self.page # верхний предел по x
        self.xstep = self.page//2     # шаг листания        
        
    def DoDraw(self):        
        self.ax.cla()
        self.ax.set_xlim([self.xstart, self.xstop])
        # увеличение с сохранением среднего
        y0 = self.ymiddle + self.zoom*(self.yrange[0]-self.ymiddle)
        y1 = self.ymiddle + self.zoom*(self.yrange[1]-self.ymiddle)
        self.ax.set_ylim([y0, y1]) 
        for i in range(self.nlines):    
            self.ax.plot( self.xx[i], self.yy[i], label=self.labels[i] ) #, xmaxs,ymaxs,'x')   
        for i in range(self.nmarkers):    
            self.ax.plot( self.xm[i], self.ym[i], self.ms[i], label=self.mlabels[i], color=self.mcolors[i] )        
        self.ax.grid(True)
        leg = self.ax.legend(fancybox=True, shadow=True, loc=4)
        self.fig.canvas.draw() 
        
    def key_event(self,e):
        if e.key == "z":
            self.page = self.page//2
        if e.key == "x":
            self.page = 2*self.page
        self.xstep = self.page//2
        if e.key == "up":
            self.zoom = 0.5*self.zoom
        if e.key == "down":
            self.zoom = 2.0*self.zoom
        if e.key == "right":
            self.xstart = self.xstart + self.xstep
        if e.key == "left":
            self.xstart = self.xstart - self.xstep
        self.xstop = self.xstart + self.page       
        self.DoDraw()
        
    def AddLine( self, cx, cy, label): 
        self.xx.append(cx) 
        self.yy.append(cy) 
        self.labels.append(label)
        self.nlines = self.nlines + 1
 
    def AddMarkers( self, cx, cy, m, label, colour ): 
        self.xm.append(cx) 
        self.ym.append(cy) 
        self.mlabels.append(label)
        self.mcolors.append(colour)
        self.ms.append(m)
        self.nmarkers = self.nmarkers + 1
        
    def Draw(self, filename):
        def gkey_event(e):
            self.key_event(e)
        self.fig.canvas.mpl_connect('key_press_event', gkey_event )
        self.fig.canvas.manager.set_window_title(filename)
        self.DoDraw()
        plt.title(label=filename, fontsize=25)
        plt.show()
           
def TestPoxView():
    #ecg0 = electrocardiogram()
    #ecg = ecg0[0:5000]    
    npoints = 100   #ecg.size
    x = 0.2*np.arange(npoints)
    y = x*x
    #peaks, _ = find_peaks(ecg, height=1.0)
    zviewer = ZoomViewer(20,[-100.0,100.0])
    zviewer.AddLine( x, y, "x*x" )    
    #zviewer.AddMarkers( peaks, ecg[peaks], '+', "peaks" )
    zviewer.Draw("test")
    
if __name__ == "__main__":
    TestPoxView()   
    
    
    