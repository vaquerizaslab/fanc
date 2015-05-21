'''
Created on May 20, 2015

@author: kkruse1
'''

from rpy2.robjects import pandas2ri as p2r
from rpy2.robjects.packages import importr

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
    def __init__(self, hic, chrom=None, start=None, end=None):
        self.df = hic.as_data_frame()
        self.chrom = chrom
        self.start = start
        self.end = end
    
    def show(self, output=None):
        p2r.activate()
        sushi = importr('Sushi')
        
        dfr = p2r.py2ri(self.df)
        dfr.colnames = dfr.rownames # to correct for leading "X" in colnames
        
        if output:
            open_graphics_file(output)
        
        sushi.plotHic(dfr,self.chrom,self.start,self.end)
        
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
            
        self.plots.append(panel)
        
    def show(self, output):
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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        