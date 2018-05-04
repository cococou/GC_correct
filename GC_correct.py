#!/usr/bin/env python3
import os,sys
import pysam
import math
import argparse
import numpy as np
from scipy.optimize import leastsq
import matplotlib.pyplot as plt
from pandas.core.frame import DataFrame
import re
import random

plt.switch_backend('agg')
__author__ = 'Marco Polo'

def PAR():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b","--bam",help="sorted index bamfile",required=True)
    parser.add_argument("-u",'--um',help="target region,bed file")
    parser.add_argument("-w",'--window',help="cnv correct bin size",default=1000000,type=int)
    parser.add_argument("-f",'--df',help="degree of freedom",default=4,type=int)
    parser.add_argument("-r",'--rm',help="remove region,bed file,such as um35-hs37d5",default='NA')
    parser.add_argument("-c",'--chrom',help="plot one chromsome,such as [chr2]",default="NA")
    args = parser.parse_args()
    return args

#def xopen_rm(par):
#    rm_region = None
#    if par.rm == 'NA':
#        return rm_region
#    names = ['chr','start','end']
#    dfty = {'chr':'str','start':'int','end':'int'}
#    rm_region = (pd.read_table(par.rm,sep="\t",names=names,header=False)
#                    
#                )
            

def xopen(par):
    # this module is used to get bed file by the window,if more than window's 20% will keep that region
    with open(par.um) as f:
        for line in f:
            regions = []
            line = line.strip().split(None)
            chrom,start,end = line[0:3]
            start = int(start);end = int(end)
            chu = (end-start)/par.window
            if chu < 0.2 :
                continue
            elif chu >= 0.2 and chu <= 1:
                regions.append((chrom,start,end))
            else:
                decimal,integer = math.modf(chu)
                integer = int(integer)
                for i in range(1,integer+1):
                    regions.append((chrom,start+par.window*(i-1),start+par.window*i-1))
                if decimal >= 0.2:
                    START = regions[len(regions)-1][1] + 1
                    END = end
                    regions.append((chrom,START,END))
            yield regions
            
                
def get_target(par,chrom,start,end,depths_average,GC_percentage,chroms,starts):
    # count indel and refskip site in depths
    samfile = pysam.AlignmentFile(par.bam, "rb")
    depths = 0;GC = 0;not_zero_pos_num = 0
    for pileupcolumn in samfile.pileup(chrom, start-1, end,truncate=True):
        pos = pileupcolumn.pos
        depths += pileupcolumn.n 
        not_zero_pos_num += 1
        #print('---',pos,depth)
        for pileupread in pileupcolumn.pileups:   
            if not pileupread.is_del and not pileupread.is_refskip:
                #depths += 1
                base = pileupread.alignment.query_sequence[pileupread.query_position]
                if base == "G" or base == "C":
                    GC += 1
                #print(base)
    if depths != 0:
        bin_depths_average = round(depths/(end-start+1),6)
        #bin_depths_average = round(depths/not_zero_pos_num,6)
        bin_GC_percentage = round((GC/depths)*100,4)
        depths_average.append(bin_depths_average)
        GC_percentage.append(bin_GC_percentage)
        chroms.append(chrom)
        starts.append(start)
        #real_win.append(end-start+1)
    return depths_average,GC_percentage,chroms,starts

def categore(depths_average,GC_percentage):
    categore_dic = {i:[] for i in range(0,101,2)} #0-100,step 2 
    for bin_dp_average,bin_GC_percentage in zip(depths_average,GC_percentage):
        bin_GC_percentage = round(bin_GC_percentage)
        if bin_GC_percentage % 2 == 1: bin_GC_percentage = bin_GC_percentage + 1
        categore_dic[bin_GC_percentage].append(bin_dp_average)
    xy = {key:float(np.median(value)) for key,value in categore_dic.items() if value}
    return xy

def polynomial(x,y,df):
    out = np.polyfit(x, y, df)
    outoj = np.poly1d(out)
    return outoj

def normalize(outoj,x):
    y = outoj(x)
    return y

def handle_interval(rm_region,chrom,start,end):
    pass

def main():
    par = PAR()
    #rm_region = xopen_rm(par)
    starts = []
    chroms =  []
    real_win = []
    depths_average = []
    GC_percentage = []
    window_warning = 0
    df = par.df
    for regions in xopen(par):
        for chrom,start,end in regions:              #window's bin       
            window_warning += 1
            depths_average,GC_percentage,chroms,starts = get_target(par,chrom,start,end,depths_average,GC_percentage,chroms,starts)
    if window_warning == 0:
        print("!!!!warning!!!!\n    your window size is too large or you bed region is too small\n    do noting and exit",file=sys.stderr)
        exit()

    xy = categore(depths_average,GC_percentage)
    #GC_categore,depths_average_median = zip(*sorted(xy.items(), key=lambda t: t[1]))
    GC_categore = sorted(xy.keys())
    depths_average_median = [xy[i] for i in GC_categore]
   
    outoj = polynomial(GC_categore,depths_average_median,df)
    depths_average_normalize = normalize(outoj,GC_percentage)
    #print(depths_average)
    #print(GC_percentage)
    #print(categore_dic)
    outfile = os.path.basename(par.bam).split('.')[0]+'.GC'
    outdf = DataFrame({'chrom':chroms,'start':starts,"real_depth_average":depths_average,"normalized_depth":depths_average_normalize})
    outdf.to_csv(outfile,index=False,header=True,sep="\t",columns = ['chrom','start','real_depth_average','normalized_depth'])

    norm = DataFrame({'GC_percentage':GC_percentage,"depths_average_normalize":depths_average_normalize})
    norm = norm.sort_values(['GC_percentage'],ascending=[True])

    fs = 8

    if par.chrom == "NA":
        pics = 3
    else:
        pics = 4
        cr = par.chrom
        y = outdf.query("chrom == @cr").normalized_depth
        x = outdf.query("chrom == @cr").start 
        ax = plt.subplot(pics, 1, pics)
        #ax.axis([0, 10*24, 0, ymax])
        ax.set_ylabel('normalised depth',fontsize=5)
        ax.set_xlabel(cr,fontsize=5)
        #ax.set_xticks([])
        #ax.set_ylim(bottom=0,top=np.mean(y)*1.5)
        ax.tick_params(labelsize = 3) 
        ax.plot(x,y,'bo',markersize=0.4)

    ax = plt.subplot(pics, 1, 1)
    
    ax.set_xlabel('polynomial GC',fontsize=fs)
    ax.set_ylabel('depth',fontsize=fs)
    ax.plot(GC_categore, depths_average_median, 'r--', label='Origin Line',linewidth=1)
    ax.plot(norm.GC_percentage, norm.depths_average_normalize, 'g-', label='Poly Fitting Line(deg={df})'.format(df=df),linewidth=1)
    xmax = np.max(GC_categore) + np.max(GC_categore) * 0.2
    ymax = np.max(depths_average_normalize + depths_average_normalize) * 0.8
    
    ax.tick_params(labelsize = fs)
    ax.axis([0, xmax, 0, ymax])
    ax.legend(fontsize=fs)
    ax.set_xlabel('GC percentage',fontsize=fs)

    ax = plt.subplot(pics, 1, 2)
    i = -10
    nn = -1
    for cr in outdf.chrom.unique():
        i += 10; nn += 1
        y = outdf.query("chrom == @cr").real_depth_average
        x = [random.uniform(i,i+9) for n in range(len(y))]
        if nn % 2 == 0: 
            col = 'ro'
        else:
            col = 'bo'
        ax.axis([0, 10*24, 0, ymax])
        ax.set_ylabel('real depth',fontsize=fs)
        ax.set_xlabel('chromosome',fontsize=fs)
        ax.set_xticks([])
        #ax.set_ylim(bottom=0,top=np.mean(y)*1.5)
        ax.tick_params(labelsize = fs)
        ax.plot(x,y,col,markersize=0.4)
    
    ax = plt.subplot(pics, 1, 3)
    i = -10 
    nn = -1
    for cr in outdf.chrom.unique():
        i += 10; nn += 1
        y = outdf.query("chrom == @cr").normalized_depth
        x = [random.uniform(i,i+9) for n in range(len(y))]
        if nn % 2 == 0:  
            col = 'ro'
        else:
            col = 'bo'
        ax.axis([0, 10*24, 0, ymax])
        ax.set_ylabel('normalised depth',fontsize=fs)
        ax.set_xlabel('chromosome',fontsize=fs)
        ax.set_xticks([])
        #ax.set_ylim(bottom=0,top=np.mean(y)*1.5)
        ax.tick_params(labelsize = fs)
        ax.plot(x,y,col,markersize=0.4)
        #ax.xticks(None)
            
    plt.savefig(os.path.basename(par.bam).split('.')[0]+'.SVG')
    plt.close()
    
    #for i,j in zip(depths_average,depths_average_normalize):
    #    print(i,j,sep="\t")
    #print("totoal {num} windows".format(num=window_warning),file=sys.stderr)
    #print(xy)
    print(outoj)
   
    #print(norm.GC_percentage)
    #print(depths_average_median)


if __name__ == "__main__":
    main()
