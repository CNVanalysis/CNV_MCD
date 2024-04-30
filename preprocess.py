import os

import pysam
import numpy as np
from matplotlib import pyplot as plt
from rpy2 import robjects
from sklearn.mixture import BayesianGaussianMixture
from sklearn.neighbors import KernelDensity
from sklearn.cluster import KMeans
from sklearn.utils.validation import check_array
import pandas as pd
from sklearn import mixture
import itertools
from scipy import linalg
import matplotlib as mpl


def get_chrlist(filename):
    samfile = pysam.AlignmentFile(filename, "rb", ignore_truncation=True)
    List = samfile.references
    chrList = np.full(len(List), 0)
    for i in range(len(List)):
        chr = str(List[i]).strip('chr')
        if chr.isdigit():
            chrList[i] = int(chr)

    index = chrList > 0
    chrList = chrList[index]
    return chrList


def read_ref_file(filename, ref, num):
    # read reference file
    if os.path.exists(filename):
        print("Read reference file: " + str(filename))
        with open(filename, 'r') as f:
            line = f.readline()
            for line in f:
                linestr = line.strip()
                ref[num] += linestr
    else:
        print("Warning: can not open " + str(filename) + '\n')
    return ref


def get_RC(filename, chrList, ReadCount):
    samfile = pysam.AlignmentFile(filename, "rb", ignore_truncation=True)
    for line in samfile:
        if line.reference_name:
            chr = line.reference_name.strip('chr')
            if chr.isdigit():
                num = np.argwhere(chrList == int(chr))[0][0]
                posList = line.positions
                ReadCount[num][posList] += 1
    return ReadCount


def ReadDepth(ReadCount, binNum, ref, binSize):
    RD = np.full(binNum, 0.0)
    GC = np.full(binNum, 0)
    pos = np.arange(1, binNum + 1)
    print(str(binNum) + "binNum")
    for i in range(binNum):
        RD[i] = np.mean(ReadCount[i * binSize:(i + 1) * binSize])
        cur_ref = ref[i * binSize:(i + 1) * binSize]
        N_count = cur_ref.count('N') + cur_ref.count('n')
        if N_count == 0:
            gc_count = cur_ref.count('C') + cur_ref.count('c') + cur_ref.count('G') + cur_ref.count('g')
        else:
            RD[i] = -10000
            gc_count = 0
        GC[i] = int(round(gc_count / binSize, 3) * 1000)

    index = RD >= 0
    RD = RD[index]
    GC = GC[index]
    pos = pos[index]
    RD = gc_correct(RD, GC)
    return pos, RD


def gc_correct(RD, GC):
    bincount = np.bincount(GC)
    global_rd_ave = np.mean(RD)
    for i in range(len(RD)):
        if bincount[GC[i]] < 2:
            continue
        mean = np.mean(RD[GC == GC[i]])
        RD[i] = global_rd_ave * RD[i] / mean
    return RD


def plot(pos, data):
    plt.scatter(pos, data, s=5, c="red")
    plt.xlabel("pos")
    plt.ylabel("RD")
    plt.show()


def modeRD(RD):

    newRD = np.full(len(RD), 0)
    for i in range(len(RD)):
        newRD[i] = int(round(RD[i], 3) * 1000)
    count = np.bincount(newRD)
    countList = np.full(len(count) - 49, 0)
    for i in range(len(countList)):
        countList[i] = np.mean(count[i:i + 50])
    modemin = np.argmax(countList)
    modemax = modemin + 50
    mode = (modemax + modemin) / 2
    mode = mode / 1000
    return mode


def scaling_RD(RD, mode):
    for i in range(len(RD)):
        try:
            RD[i] = np.math.log(RD[i], mode) - 1
        except ValueError as e:
            print(str(i))
            print('RD[i]: ' + str(RD[i]) + " mode: " + str(mode))
            print(e)
        except ZeroDivisionError as e:
            print(str(i))
            print('RD[i]: ' + str(RD[i]) + " mode: " + str(mode))
            print(e)
    return RD


def preprocess(bam, binSize, chr_num):
    # chrList = get_chrlist(bam)
    # print(chrList)
    chrList = np.full(1, 0)
    chrList[0] = int(chr_num)
    print(chrList)

    chrNum = len(chrList)
    refList = [[] for i in range(chrNum)]

    for i in range(chrNum):
        refList = read_ref_file("../hg19/chr" + str(chr_num) + ".fa", refList,
                                i)
    chrLen = np.full(chrNum, 0)
    modeList = np.full(chrNum, 0.0)

    for i in range(chrNum):
        chrLen[i] = len(refList[i])
        print("the", i + 1, "base number of ref file：", chrLen[i])

    print("Read train bam file:", bam)

    ReadCount = np.full((chrNum, np.max(chrLen)), 0)
    ReadCount = get_RC(bam, chrList, ReadCount)

    for i in range(chrNum):
        binNum = int(chrLen[i] / binSize) + 1
        pos, RD = ReadDepth(ReadCount[0], binNum, refList[i], binSize)
        numbin = len(RD)
        for m in range(numbin):
            if np.isnan(RD[m]).any():
                if (m == numbin - 1):
                    RD[m] = RD[m - 1]
                else:
                    RD[m] = (RD[m - 1] + RD[m + 1]) / 2

        return pos, RD, RD

        # #CBS---------------------------------------------
        # print("segment count...")
        # v = robjects.FloatVector(scalRD)
        # col = round(chrLen[i] / 500000)
        # m = robjects.r['matrix'](v, ncol=col)
        # robjects.r.source("CBS_data.R")
        # robjects.r.CBS_data(m, os.path.abspath('..') + str("/seg"))
        #
        # # subprocess.call('Rscript CBS_data.R',shell=True)
        # num_col = int(numbin / col) + 1
        # seg_start, seg_end, seg_count, seg_len = Read_seg_file(num_col, numbin)
        # seg_count = np.array(seg_count)
        # seg_count, seg_start, seg_end = seg_RD(RD, pos, seg_start, seg_end, seg_count, binSize)
        #
        # all_RD = []
        # all_start = []
        # all_end = []
        # all_RD.extend(seg_count)
        # all_start.extend(seg_start)
        # all_end.extend(seg_end)
        # return all_start, all_end, all_RD
        # # CBS---------------------------------------------


def Read_seg_file(num_col, num_bin):
    """
    read segment file (Generated by DNAcopy.segment)
    seg file: col, chr, start, end, num_mark, seg_mean
    """
    seg_start = []
    seg_end = []
    seg_count = []
    seg_len = []
    with open("seg", 'r') as f:
        for line in f:
            linestrlist = line.strip().split('\t')
            start = (int(linestrlist[0]) - 1) * num_col + int(linestrlist[2]) - 1  # start从0开始
            end = (int(linestrlist[0]) - 1) * num_col + int(linestrlist[3]) - 1
            if start < num_bin:
                if end > num_bin:
                    end = num_bin - 1
                seg_start.append(start)
                seg_end.append(end)
                seg_count.append(float(linestrlist[5]))
                seg_len.append(int(linestrlist[4]))
    seg_start = np.array(seg_start)
    seg_end = np.array(seg_end)

    return seg_start, seg_end, seg_count, seg_len


def seg_RD(RD, binHead, seg_start, seg_end, seg_count, binSize):
    seg_RD = np.full(len(seg_count), 0.0)
    for i in range(len(seg_RD)):
        if seg_start[i] == seg_end[i]:
            seg_RD[i] = seg_RD[i - 1]
        else:
            seg_RD[i] = np.mean(RD[seg_start[i]:seg_end[i]])
        seg_start[i] = binHead[seg_start[i]] * binSize + 1
        if seg_end[i] == len(binHead):
            seg_end[i] = len(binHead) - 1
        seg_end[i] = binHead[seg_end[i]] * binSize + binSize
    return seg_RD, seg_start, seg_end


def read_seg_file(filename):
    segRD = []
    with open(filename, 'r') as f:
        line = f.readline()
        for line in f:
            linestr = line.strip()
            linestrlist = linestr.split('\t')
            segRD.append(float(linestrlist[1]))
    segRD = np.array(segRD)
    return segRD


def segment(pos, segrd):
    start = []
    end = []
    seg_rd = []
    i = 0
    j = 1
    while j < len(segrd):
        if j == len(segrd) - 1:
            start.append(int(pos[i]))
            end.append(int(pos[j - 1]))
            seg_rd.append(float(segrd[i]))
            j += 1
        else:
            if segrd[i] == segrd[j]:
                j += 1
            else:
                start.append(int(pos[i]))
                end.append(int(pos[j - 1]))
                seg_rd.append(float(segrd[i]))
                i = j
    return start, end, seg_rd
