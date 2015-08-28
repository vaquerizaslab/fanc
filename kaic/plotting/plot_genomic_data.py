'''
Created on May 20, 2015

@author: kkruse1
'''

from rpy2.robjects import pandas2ri as p2r
from rpy2.robjects.packages import importr
import numpy as np

def open_graphics_file(file_name, width=None, height=None):
    p2r.activate()
    gr = importr('grDevices')
    if file_name.endswith('.pdf'):
        if width is not None and height is not None:
            gr.pdf(file_name,width=width,height=height)
        elif width is not None:
            gr.pdf(file_name,width=width)
        elif height is not None:
            gr.pdf(file_name,height=height)
        else:
            gr.pdf(file_name)
    elif file_name.endswith('.png'):
        if width is not None and height is not None:
            gr.png(file_name,width=width,height=height)
        elif width is not None:
            gr.png(file_name,width=width)
        elif height is not None:
            gr.png(file_name,height=height)
        else:
            gr.png(file_name)
    elif file_name.endswith('.svg'):
        if width is not None and height is not None:
            gr.svg(file_name,width=width,height=height)
        elif width is not None:
            gr.svg(file_name,width=width)
        elif height is not None:
            gr.svg(file_name,height=height)
        else:
            gr.svg(file_name)
    elif file_name.endswith('.jpg') or file_name.endswith('.jpeg'):
        if width is not None and height is not None:
            gr.jpeg(file_name,width=width,height=height)
        elif width is not None:
            gr.jpeg(file_name,width=width)
        elif height is not None:
            gr.jpeg(file_name,height=height)
        else:
            gr.jpeg(file_name)
    else:
        raise ValueError("File ending not supported: " + file_name + " (try pdf, svg, png, or jpg)")

def close_graphics_file():
    p2r.activate()
    gr = importr('grDevices')
    gr.dev_off()
    



class BedAlignment(object):
    def __init__(self, bed1, bed2, chrom=None, start=None, end=None, n_bins=None, window_size=None):
        self.bed1 = bed1
        self.bed2 = bed2
        self.chrom = chrom
        self.start = start
        self.end = end
        self.n_bins = n_bins
        self.window_size = window_size
        
    def show(self, output=None):
        p2r.activate()
        genomation = importr('genomation')
        genomicRanges = importr('GenomicRanges')
        graphics = importr('graphics')
        
        bed1_df = self.bed1.as_data_frame(self.chrom,self.start,self.end)
        bed2_df = self.bed2.as_data_frame(self.chrom,self.start,self.end)
        
        # correct data frame
        if self.window_size is not None:
            mean = map(int, (bed2_df["start"]+bed2_df["end"])/2)
            bed2_df["start"] = [x - int(self.window_size/2) for x in mean]
            bed2_df["end"]   = [x + int(self.window_size/2) for x in mean]
        
        if output:
            open_graphics_file(output)
        
        if bed1_df.shape[0] > 0:
            # convert to R objects
            bed1_dfr = p2r.py2ri(bed1_df)
            bed2_dfr = p2r.py2ri(bed2_df)
            
            # get GenomicRange objects
            bed1_range = genomicRanges.makeGRangesFromDataFrame(bed1_dfr)
            bed2_range = genomicRanges.makeGRangesFromDataFrame(bed2_dfr)
            
            # get score matrix
            if self.n_bins is not None:
                sm = genomation.ScoreMatrixBin(target=bed1_range, windows=bed2_range, bin_num=self.n_bins)
            else:
                sm = genomation.ScoreMatrix(target=bed1_range, windows=bed2_range)
            
            genomation.heatMatrix(sm)
        else:
            #empty plot
            graphics.plot(0,type='n',axes=False,ann=False)

        if output:
            close_graphics_file()


class BedDistribution(object):
    def __init__(self, bed1, bed2, chrom=None, start=None, end=None, n_bins=None, window_size=None, rescale=True):
        self.bed1 = bed1
        if type(bed2) is not list:
            self.bed2 = [bed2]
        else:
            self.bed2 = bed2
        self.chrom = chrom
        self.start = start
        self.end = end
        self.n_bins = n_bins
        self.window_size = window_size
        self.rescale=rescale
        
    def show(self, output=None):
        p2r.activate()
        base = importr('base')
        genomation = importr('genomation')
        genomicRanges = importr('GenomicRanges')
        graphics = importr('graphics')
        grDevices = importr('grDevices')
        
        bed1_df = self.bed1.as_data_frame(self.chrom,self.start,self.end)
        
        if output:
            open_graphics_file(output)
        
        names = []
        if bed1_df.shape[0] > 0:
            bed1_dfr = p2r.py2ri(bed1_df)
            bed1_range = genomicRanges.makeGRangesFromDataFrame(bed1_dfr)
            
            l = base.vector(mode="list",length=len(self.bed2))
            i = 0
            for bed2 in self.bed2:
                bed2_df = bed2.as_data_frame(self.chrom,self.start,self.end)
                
                # correct data frame
                if self.window_size is not None:
                    mean = map(int, (bed2_df["start"]+bed2_df["end"])/2)
                    bed2_df["start"] = [x - int(self.window_size/2) for x in mean]
                    bed2_df["end"]   = [x + int(self.window_size/2) for x in mean]
                
                bed2_dfr = p2r.py2ri(bed2_df)
                bed2_range = genomicRanges.makeGRangesFromDataFrame(bed2_dfr)
                
                # get score matrix
                if self.n_bins is not None:
                    sm = genomation.ScoreMatrixBin(target=bed1_range, windows=bed2_range, bin_num=self.n_bins)
                else:
                    sm = genomation.ScoreMatrix(target=bed1_range, windows=bed2_range)
                    
                l[i] = sm
                names.append(bed2.name)
                i += 1
                
            l.names = names
            sml = genomation.ScoreMatrixList(l)
            
            # prepare plot
            graphics.layout(base.rbind(1,2), heights=base.c(4,1))
            
            graphics.par(mar=base.c(5,4,4,0), oma=base.c(0,0,0,0))
            colors = grDevices.rainbow(len(names))
            genomation.plotMeta(sml, meta_rescale=self.rescale, line_col=colors)
            if self.n_bins is not None:
                graphics.abline(v=self.n_bins/2)
            
            graphics.par(mar=base.c(0,0,0,0))
            graphics.plot_new()
            graphics.legend(x="center", legend=base.c(names), lty=base.c(1,1), lwd=base.c(2.5,2.5), col=colors, ncol=2, bty="n")
            
        else:
            #empty plot
            graphics.plot(0,type='n',axes=False,ann=False)

        if output:
            close_graphics_file()




class HiCPlotOld(object):
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
            
class HiCPlot(object):
    def __init__(self, hic, chrom=None, start=None, end=None,
                 zrange=[5,68], max_y=20, colors=["white","blue"]):
        self.hic = hic
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
        
        df = self.hic.as_data_frame(('%s:%d-%d' % (self.chrom, self.start, self.end),
                                     '%s:%d-%d' % (self.chrom, self.start, self.end)))
        
        if output:
            open_graphics_file(output)
        
        if df.shape[0] > 0:
            dfr = p2r.py2ri(df)
            dfr.colnames = dfr.rownames # to correct for leading "X" in colnames
            print dfr
            sushi.plotHic(dfr,self.chrom,self.start,self.end,
                          palette=grd.colorRampPalette(self.colors),
                          zrange=self.zrange, max_y=self.max_y)
            print "plotted sushi"
        else:
            print "no plot generated"
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

class GenePlot(object):
    def __init__(self, bed, chrom=None, start=None, end=None,
                 maxrows=70,plotgenetype="box", labeltext=True,
                 showCoordinates=True):
        self.data = bed
        self.chrom = chrom
        self.start = start
        self.end = end
        self.maxrows = maxrows
        self.plotgenetype = plotgenetype
        self.labeltext = labeltext
        self.showCoordinates = showCoordinates
        
    def show(self, output=None):
        p2r.activate()
        sushi = importr('Sushi')
        graphics = importr('graphics')
        
        df = self.data.as_data_frame(self.chrom,self.start,self.end)
        types = np.array(['exon'] * len(df))
        
        if output:
            open_graphics_file(output)
        
        if df.shape[0] > 0:
            dfr = p2r.py2ri(df)
            sushi.plotGenes(dfr,self.chrom,self.start,self.end,
                            maxrows=self.maxrows, plotgenetype=self.plotgenetype,
                            labeltext=self.labeltext, types=types)
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


# class HiCCorrelationPlot(object):
#     def __init__(self, hic1, hic2):
#         self.hic1 = hic1
#         self.hic2 = hic2
#         
#     def show(self, output=None):
#         p2r.activate()
#         graphics = importr('graphics')
#         base = importr('base')
#         
#         if output:
#             open_graphics_file(output)
# 
# 
#         l = len(self.panels)
#         graphics.layout(base.matrix(range(1,l+1), l, 1, byrow=True))
#         graphics.par([3,4,1,1])
#         
#         for panel in self.panels:
#             panel.show()
# 
#         
#         if output:
#             close_graphics_file()


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
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        