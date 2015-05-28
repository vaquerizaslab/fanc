'''
Created on May 20, 2015

@author: kkruse1
'''

from rpy2.robjects import pandas2ri as p2r
from rpy2.robjects.packages import importr
import numpy as np

def open_graphics_file(file_name):
    p2r.activate()
    gr = importr('grDevices')
    if file_name.endswith('.pdf'):
        gr.pdf(file_name)
    elif file_name.endswith('.png'):
        gr.png(file_name)
    elif file_name.endswith('.svg'):
        gr.svg(file_name)
    elif file_name.endswith('.jpg') or file_name.endswith('.jpeg'):
        gr.jpeg(file_name)
    else:
        raise ValueError("File ending not supported: " + file_name + " (try pdf, svg, png, or jpg)")

def close_graphics_file():
    p2r.activate()
    gr = importr('grDevices')
    gr.dev_off()

class HiCPlot(object):
    def __init__(self, hic, resolution, chrom=None, start=None, end=None,
                 zrange=[5,68], max_y=20, colors=["white","blue"]):
        self.data = hic
        self.resolution = resolution
        self.chrom = chrom
        self.start = start
        self.end = end
        self.zrange = np.array(zrange)
        self.max_y = max_y
        self.colors=np.array(colors)
    
    def show(self, output=None, showCoordinates=True):
        p2r.activate()
        sushi = importr('Sushi')
        grd = importr('grDevices')
        graphics = importr('graphics')
        
        df = self.data.as_data_frame(self.resolution, self.chrom, self.start, self.end)
        
        if output:
            open_graphics_file(output)
        
        if df.shape[0] > 0:
            dfr = p2r.py2ri(df)
            dfr.colnames = dfr.rownames # to correct for leading "X" in colnames
            sushi.plotHic(dfr,self.chrom,self.start,self.end,
                          palette=grd.colorRampPalette(self.colors),
                          zrange=self.zrange, max_y=self.max_y)
        else:
            #empty plot
            graphics.plot(0,type='n',axes=False,ann=False)
            
        if showCoordinates:
            sushi.labelgenome(self.chrom,chromstart=self.start,chromend=self.end,n=4,scale="Mb")
        
        if output:
            close_graphics_file()


class BedPlot(object):
    def __init__(self, bed, chrom=None, start=None, end=None,
                 plotType="region", showCoordinates=True):
        self.data = bed
        self.chrom = chrom
        self.start = start
        self.end = end
        self.type = plotType
        self.showCoordinates = showCoordinates
        
    def show(self, output=None):
        p2r.activate()
        sushi = importr('Sushi')
        graphics = importr('graphics')
        
        df = self.data.as_data_frame(self.chrom,self.start,self.end)
        
        if output:
            open_graphics_file(output)
        
        if df.shape[0] > 0:
            dfr = p2r.py2ri(df)
            sushi.plotBed(dfr,self.chrom,self.start,self.end, type=self.type)
        else:
            #empty plot
            graphics.plot(0,type='n',axes=False,ann=False)
            
            
        if self.showCoordinates:
            sushi.labelgenome(self.chrom,chromstart=self.start,chromend=self.end,n=4,scale="Mb")
        
        if output:
            close_graphics_file()


class BedpePlot(object):
    def __init__(self, bedpe, chrom=None, start=None, end=None,
                 plotType="loops", showCoordinates=True, heights=None):
        self.data = bedpe
        self.chrom = chrom
        self.start = start
        self.end = end
        self.type = plotType
        self.showCoordinates = showCoordinates
        self.heights = heights
        
    def show(self, output=None):
        p2r.activate()
        sushi = importr('Sushi')
        graphics = importr('graphics')
        
        df = self.data.as_data_frame(self.chrom,self.start,self.end)
        if not self.heights:
            self.heights = np.ones(df.shape[0])
            
        if output:
            open_graphics_file(output)
        
        if df.shape[0] > 0:
            dfr = p2r.py2ri(df)
            sushi.plotBedpe(dfr,self.chrom,self.start,self.end, plottype=self.type,heights=self.heights)
        else:
            #empty plot
            graphics.plot(0,type='n',axes=False,ann=False)
        
        if self.showCoordinates:
            sushi.labelgenome(self.chrom,chromstart=self.start,chromend=self.end,n=4,scale="Mb")
        
        if output:
            close_graphics_file()



class GenomicDataPlot(object):
    
    def __init__(self, chrom=None, start=None, end=None):
        self.start = start
        self.end = end 
        self.chrom = chrom
        
        
        self.panels = []
    
    def add_panel(self, panel):
        if not panel.chrom:
            panel.chrom = self.chrom
        if not panel.start:
            panel.start = self.start
        if not panel.end:
            panel.end = self.end
            
        self.panels.append(panel)
        
    def show(self, output=None):
        p2r.activate()
        graphics = importr('graphics')
        base = importr('base')
        
        if output:
            open_graphics_file(output)


        l = len(self.panels)
        graphics.layout(base.matrix(range(1,l+1), l, 1, byrow=True))
        graphics.par([3,4,1,1])
        
        for panel in self.panels:
            panel.show()

        
        if output:
            close_graphics_file()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        