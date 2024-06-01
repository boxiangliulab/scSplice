import time
import argparse
import pysam
import os
import random
import gzip
from numba import jit
from sklearn.decomposition import PCA
from sklearn import preprocessing
from sklearn import linear_model
import itertools as it
from scipy.stats import rankdata
from scipy.stats import norm
import sys

def readin_junction(sample):
    clu_junc = {}
    head = True
    num = 0
    f = gzip.open(sample,'rb')
    for line in f.readlines():
        i = line.decode()
        if head == True:
            head = False
            continue
        chr_clu = i.split(' ')[0]
        chr, start, end, clu = chr_clu.split(':')
        clu_junc[(chr,start,end,clu)] = []
        num = len(i.split(' ')[1:])
        for j in i.split(' ')[1:]:
            clu_junc[(chr,start,end,clu)].append(int(j))
    return num, clu_junc
    


def compute_PSI(splice_site,sample,run_dir,prefix):
    PSI = {}
    splicesite = []
    num, clu_junc = readin_junction(sample)
    for ln in open(splice_site+"_include"):
        ln = ln.strip()
        chr_strand = ln.split(' ')[0]
        chr = chr_strand.split(':')[0]
        strand = chr_strand.split(':')[1]
        clu = ln.split(' ')[1]
        newclu = 'clu_%s' % clu + '_%s' % strand
        if len(ln.split(' ')) < 3:
            continue
        event = ln.split(' ')[2:]
        for e in event:
            A,B,C = e.split(':')
            splicesite.append((chr,newclu,A))
    for ln in open(splice_site+"_exclude"):
        ln = ln.strip()
        chr_strand = ln.split(' ')[0]
        chr = chr_strand.split(':')[0]
        strand = chr_strand.split(':')[1]
        clu = ln.split(' ')[1]
        newclu = 'clu_%s' % clu + '_%s' % strand
        if len(ln.split(' ')) < 3:
            continue
        event = ln.split(' ')[2:]
        for e in event:
            A,B,C = e.split(':')
            splicesite.append((chr,newclu,A))
    desplicesite=list(set(splicesite))
    desplicesite.sort()
    include_count={}
    exclude_count={}
    for item in desplicesite:
        include_count[item]=[]
        exclude_count[item]=[]
        for ind in range(0, int(num) , 1):
            include_count[item].append(0)
            exclude_count[item].append(0)
    for ln in open(splice_site+"_include"):
        ln = ln.strip()
        chr_strand = ln.split(' ')[0]
        chr = chr_strand.split(':')[0]
        strand = chr_strand.split(':')[1]
        clu = ln.split(' ')[1]
        newclu = 'clu_%s' % clu + '_%s' % strand
        print(newclu)
        if len(ln.split(' ')) < 3:
            continue
        event = ln.split(' ')[2:]
        for e in event:
            A,B,C = e.split(':')
            if (chr,newclu,A) in desplicesite:
                for ind in range(0, int(num) , 1):
                    inclusion = int(clu_junc[(chr,B,C,newclu)][ind]) 
                    include_count[(chr,newclu,A)][ind]=include_count[(chr,newclu,A)][ind]+inclusion
    for ln in open(splice_site+"_exclude"):
        ln = ln.strip()
        chr_strand = ln.split(' ')[0]
        chr = chr_strand.split(':')[0]
        strand = chr_strand.split(':')[1]
        clu = ln.split(' ')[1]
        newclu = 'clu_%s' % clu + '_%s' % strand
        print(newclu)
        if len(ln.split(' ')) < 3:
            continue
        event = ln.split(' ')[2:]
        for e in event:
            A,B,C = e.split(':')
            if (chr,newclu,A) in desplicesite:
                for ind in range(0, int(num) , 1):
                    print(ind)
                    exclusion = int(clu_junc[(chr,B,C,newclu)][ind]) 
                    exclude_count[(chr,newclu,A)][ind]=exclude_count[(chr,newclu,A)][ind]+exclusion
    for item in splicesite:
        PSI[item]=[]
        for ind in range(0,int(num),1):
            if exclude_count[item][ind]+include_count[item][ind] > 0:
               PSI[item].append(include_count[item][ind]/(exclude_count[item][ind]+include_count[item][ind]))
            if exclude_count[item][ind]+include_count[item][ind] == 0:
               PSI[item].append(-1)
    fout = open(run_dir + "/" + prefix + "_PSI",'w')
    for item in PSI:
        chr,clu,A = item
        out = '%s:%s:%s' % (chr,clu,A)
        for j in PSI[item]:
            out += ' %s' % str(j)
        fout.write(out + '\n')


def main():
    parser = argparse.ArgumentParser(description='Set file path')
    parser.add_argument('--sample', help='source sample junction read count file')
    parser.add_argument('--splice_site', help='source splice site include & exclude file')
    parser.add_argument('--run_dir', help='source run directory')
    parser.add_argument('--prefix', help='outfile prefix')

    args = parser.parse_args()
    sys.stderr.write("Start processing...\n")

    if args.splice_site is None or args.sample is None:
        sys.stderr.write("Error: no event file provided...\n")
        exit(0)
    compute_PSI(args.splice_site,args.sample,args.run_dir,args.prefix)


if __name__ == "__main__":

    main()
