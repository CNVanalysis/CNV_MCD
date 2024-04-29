import datetime
import numpy as np

from rpy2 import robjects
from rpy2.rinterface_lib.embedded import RRuntimeError

from matplotlib import pyplot as plt

from preprocess import preprocess, read_seg_file, modeRD, scaling_RD
from pyod.models.mcd import MCD
import pandas as pd
from denoising_cores import median_filter, detrend_score, remove_baseline_drift
from scipy.stats import chi2


def detect_cnv(train_bam_path, result_output_path, chr_num):
    starttime = datetime.datetime.now()

    binSize = 1000
    train_bam = train_bam_path
    outfile = result_output_path

    pos, RD, originRD = preprocess(train_bam, binSize, chr_num)

    np.savetxt('RD.txt', RD)
    np.savetxt('pos.txt', pos)

    RD = originRD

    try:
        for i in range(10):
            v = robjects.FloatVector(RD)
            m = robjects.r['matrix'](v, ncol=1)
            robjects.r.source("segment.R")
            robjects.r.segment(m)
            # segFile = "seg" + str(i) + ".txt"
            segFile = "seg.txt"
            RD = read_seg_file(segFile)
    except RRuntimeError as e:
        print(e)

    segRD = RD

    segRD = segRD + 1
    md = modeRD(segRD)
    print("mode:" + str(md))
    segRD = scaling_RD(segRD, md)

    train_pos = (pos - pos.min()) / (pos.max() - pos.min())
    train_seg = (segRD - segRD.min()) / (segRD.max() - segRD.min())
    segRD = train_seg

    X_train = np.c_[train_pos, train_seg]
    clf = MCD()
    clf.fit(X_train)
    clf_score = clf.decision_scores_
    clf_score = median_filter(clf_score)
    # clf_score = detrend_score(clf_score)
    # t1 = int(0.2 * 500)
    # clf_score = remove_baseline_drift(clf_score, t1)
    cnv_index = []
    n = 2  # choose the degrees of freedom, typically  between 1 and 3
    thres = chi2.ppf(0.01, n)
    index = np.where(clf_score > thres)
    cnv_index = np.union1d(cnv_index, index)
    cnv_index = cnv_index.astype('int64')

    # Visualize the location of outliers
    plt.scatter(pos, segRD,
                marker="o", c="b", s=5, label="RD")
    plt.scatter(pos[cnv_index], segRD[cnv_index], marker="o", c="r",
                s=5, label="OutlierRD")
    plt.title("MCD" + ",OutlierNum:" + str(len(cnv_index)))
    plt.legend()
    plt.show()

    cnv_start = pos[cnv_index]
    cnv_end = pos[cnv_index]
    cnv_RD = segRD[cnv_index]

    cnv_type = np.full(len(cnv_RD), 1)

    for i in range(len(cnv_index)):
        if segRD[cnv_index[i]] >= 0.05:
            cnv_type[i] = 2

    for i in range(len(cnv_index) - 1):
        if cnv_end[i] + 1 == cnv_start[i + 1] and cnv_type[i] == cnv_type[i + 1]:
            cnv_start[i + 1] = cnv_start[i]
            cnv_type[i] = 0

    index = cnv_type != 0
    cnv_type = cnv_type[index]
    cnv_start = cnv_start[index]
    cnv_end = cnv_end[index]
    cnv_start = (cnv_start - 1) * 1000 + 1
    cnv_end = cnv_end * 1000

    output = open(outfile, "w")
    for i in range(len(cnv_type)):
        if cnv_type[i] == 2:
            output.write("chr" + str(chr_num) + '\t' + str(int(cnv_start[i])) + '\t' + str(
                int(cnv_end[i])) + '\t' + str("gain") + '\n')
        else:
            output.write("chr" + str(chr_num) + '\t' + str(int(cnv_start[i])) + '\t' + str(
                int(cnv_end[i])) + '\t' + str("loss") + '\n')

    endtime = datetime.datetime.now()
    print("running time: " + str((endtime - starttime).seconds) + " seconds")


if __name__ == '__main__':
    chr_num = 21;
    train_bam_path = '/../test.bam'
    output_path = "/../test_result.txt"
    detect_cnv(train_bam_path, output_path, chr_num)
