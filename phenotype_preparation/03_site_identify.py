import time
import argparse
import pysam
import os
import random

from numba import jit
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn import linear_model
import itertools as it
from scipy.stats import rankdata
from scipy.stats import norm
import sys
import pandas as pd

df=pd.read_csv("/data/zhangyuntian/gencode32_exons.txt",sep="\t")


def detect_exon(site,chr,dedusite):
    pos = df[(df['chr']==chr) & ((df['start']==site) | (df['end']==site))]
    if len(pos)==0:
        return False
    else:
        newsite=[]
        for i in range(len(pos)):
            newsite.append(pos.iloc[i,1])
            newsite.append(pos.iloc[i,2])
        tmp = list(set(newsite))
        for i in tmp:
            if i not in dedusite:
                tmp.remove(i)
        tmp.remove(site)
        return tmp

###divide each site into three classes: 5' & 3'; middle sites with common exon start ;middle sites without common exon start
def site_classify(clu,chr):
    site = []
    for item in clu:
        start, end = item
        site.append(start)
        site.append(end)
    dedusite = list(set(site))
    dedusite.sort()
    empty_anno = {}
    anno=""
    empty_anno_related_site = {}
    for site in dedusite:
        empty_anno_related_site[site] = []
        pos = detect_exon(site,chr,dedusite)
        if pos!=False:
            t=0
            direction = ""
            for i in pos:
                anosite = detect_exon(i,chr,dedusite)
                if [i in dedusite]:
                    anno = "middle_site"
                    for j in dedusite:
                        if [j in anosite] & j!=site:
                            empty_anno_related_site[site].append(j)
                            anno = "middle_site_coexon"
                else:
                    t=t+1
                    if i<site:
                        direction="5"
                    else:
                        direction="3"
            if t==len(pos):
                anno = direction+"prime_site"
        empty_anno[site] = anno
    five=[]
    three=[]
    for site in dedusite:
        if empty_anno[site]=="5prime_site":
            five.append(site)
        elif empty_anno[site]=="3prime_site":
            three.append(site)
    for site in dedusite:
        if empty_anno[site]=="5prime_site":
            empty_anno_related_site[site]=five
            empty_anno_related_site[site].remove(site)
        elif empty_anno[site]=="3prime_site":
            empty_anno_related_site[site]=three
            empty_anno_related_site[site].remove(site)
    return empty_anno,empty_anno_related_site


def site_identify(clu,chr):   ###to compute the splice site's inclusion and exclusion rate
    site = []
    for item in clu:
        start, end = item
        site.append(start)
        site.append(end)
    dedusite = list(set(site))
    dedusite.sort()
    site_class,site_related = site_classify(clu,chr)
    if len(dedusite) >= 3:
        include=[]
        exclude=[]
        for site_middle in dedusite:     ###start divide sites into three parts
            ###First need to check if the sites belong to which part
            ###5' or 3' site
            ###middle site with common exon start
            ###middle site without common exon start
            type_site = site_class[site_middle]
            related = site_related[site_middle]
            if type_site=="middle_site":
              for i in clu:
                start, end = i
                if (start==site_middle) | (end==site_middle):
                    include.append((site_middle,start,end))
                if (start<site_middle) & (site_middle<end):
                    exclude.append((site_middle,start,end))
            elif type_site=="middle_site_coexon":
              for i in clu:
                start, end = i
                if (start==site_middle) | (end==site_middle):
                    include.append((site_middle,start,end))
                if (start<site_middle) & (site_middle<end):
                    exclude.append((site_middle,start,end))
                if (start in site_related[site_middle]) | (end in site_related[site_middle]):
                    exclude.append((site_middle,start,end))
            else:
              for i in clu:
                start, end = i
                if (start==site_middle) | (end==site_middle):
                    include.append((site_middle,start,end))
                if (start<site_middle) & (site_middle<end):
                    exclude.append((site_middle,start,end))
                if (start in site_related[site_middle]) | (end in site_related[site_middle]):
                    exclude.append((site_middle,start,end))     
        return include,exclude
    else:
        return False,False
            


def read_site_pooled(cluster_file,run_dir,prefix):
    cluIntron = {}
    clu = 1
    fout = open(run_dir + "/" + prefix+"_include",'w')
    fout1 = open(run_dir + "/" + prefix+"_exclude",'w')
    pooled_cluster = cluster_file
    for ln in open(pooled_cluster):
        chr = ln.split()[0]
        for intron in ln.split()[1:]:
            start, end, count = intron.split(":")
            if clu not in cluIntron:
                cluIntron[clu] = {}
            cluIntron[clu][(start,end)]=count
        site_include,site_exclude = site_identify(cluIntron[clu],chr)
        if site_include==False:
            head = '%s ' % chr + '%d ' % clu 
            fout.write(head+'\n')
        elif len(site_include) > 0:
            Nev = 0
            head = '%s ' % chr + '%d ' % clu
            for event in site_include:
                A,B,C = event
                head += "%d:%d:%d " % (int(A),int(B),int(C))
                Nev += 1
            fout.write(head+'\n')
        if site_exclude==False:
            head = '%s ' % chr + '%d ' % clu 
            fout1.write(head+'\n')
        elif len(site_exclude) > 0:
            Nev = 0
            head = '%s ' % chr + '%d ' % clu
            for event in site_exclude:
                A,B,C = event
                head += "%d:%d:%d " % (int(A),int(B),int(C))
                Nev += 1
            fout1.write(head+'\n')
        clu += 1
        print(clu)



def main():
    parser = argparse.ArgumentParser(description='Set file path')
    parser.add_argument('--cluster_file', help='source cluster file')
    parser.add_argument('--run_dir', help='source run directory')
    parser.add_argument('--prefix', help='outfile prefix')

    args = parser.parse_args()
    sys.stderr.write("Start processing...\n")
    
    if args.cluster_file is None:
        sys.stderr.write("Error: no cluster file provided...\n")
        exit(0)

    read_site_pooled(args.cluster_file,args.run_dir,args.prefix)
    print("Site potential finished!")
 
    

if __name__ == "__main__":

    main()




