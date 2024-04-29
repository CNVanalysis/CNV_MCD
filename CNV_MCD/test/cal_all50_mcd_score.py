# precision/sensitive

import sys
import numpy as np


def cal_score(input_file, output_file):
    result_start = []
    result_end = []
    result_type = []
    with open(input_file, 'r') as f:
        for line in f:
            linestr = line.strip()
            linestrlist = linestr.split('\t')
            result_start.append(int(linestrlist[1]))
            result_end.append(int(linestrlist[2]))
            result_type.append(linestrlist[4])

    truth_start = []
    truth_end = []
    truth_type = []
    with open("./ground_truth/SIM/GroundTruthCNV", 'r') as f:
        line = f.readline()
        for line in f:
            linestr = line.strip('\n')
            linestrlist = linestr.split('\t')
            truth_start.append(int(linestrlist[0]))
            truth_end.append(int(linestrlist[1]))
            if linestrlist[3] == 'gain':
                truth_type.append("gain")
            else:
                truth_type.append("loss")

    count = 0
    for i in range(len(result_type)):
        for j in range(len(truth_type)):
            if truth_start[j] <= result_start[i] <= truth_end[j] and truth_type[j] == result_type[i]:
                if result_end[i] <= truth_end[j]:
                    count += (result_end[i] - result_start[i] + 1)
                elif result_end[i] >= truth_end[j]:
                    count += (truth_end[j] - result_start[i] + 1)
                break
            elif truth_start[j] >= result_start[i] and truth_type[j] == result_type[i]:
                if truth_start[j] <= result_end[i] <= truth_end[j]:
                    count += (result_end[i] - truth_start[j] + 1)
                elif result_end[i] >= truth_end[j]:
                    count += (truth_end[j] - truth_start[j] + 1)
                break

    result_count = 0
    for i in range(len(result_start)):
        result_count += (result_end[i] - result_start[i] + 1)

    truth_count = 0
    for i in range(len(truth_start)):
        truth_count += (truth_end[i] - truth_start[i] + 1)

    print("count:", count)
    print("result：", result_count, "truth:", truth_count)
    output = open(output_file, "a")
    if (result_count == 0):
        output.write(str(0) + '\t' + str(count / truth_count) + '\n')
        print(str(0) + '\t' + str(count / truth_count) + '\n')
    else:
        output.write(str(count / result_count) + '\t' + str(count / truth_count) + '\n')
        print(str(count / result_count) + '\t' + str(count / truth_count) + '\n')


def cal_mean_score(scores_file, mean_score_file):
    precision = []
    sensitive = []
    result_type = []
    file = scores_file
    with open(file, 'r') as f:
        for line in f:
            linestr = line.strip()
            linestrlist = linestr.split('\t')
            if float(linestrlist[0]) > 0.1:
                print(float(linestrlist[0]), float(linestrlist[1]))
                precision.append(float(linestrlist[0]))
                sensitive.append(float(linestrlist[1]))

    precision = np.array(precision)
    sensitive = np.array(sensitive)
    print(len(precision), np.mean(precision), np.mean(sensitive))
    output = open(mean_score_file, "a")
    output.write(str(len(precision)) + '\t' + str(np.mean(precision)) + '\t' + str(np.mean(sensitive)))


if __name__ == '__main__':
    output_file = "/.. /scores.txt"
    mean_score_file = "/.. /mean_score.txt"
    for i in range(1, 51):
        input_file = "/.. /sim1x_" + str(i) + "_sorted.bam"
        cal_score(input_file, output_file)
    cal_mean_score(output_file, mean_score_file)
