import sys
import numpy as np
import re

if __name__ == '__main__':
    # init counter
    count_1 = 0
    count_2 = 0
    com_count = 0
    overlapratio = 0.5

    file1_start = []
    file1_end = []
    file1_type = []
    file2_start = []
    file2_end = []
    file2_type = []

    cnvnator = '/.. /cnvnator.txt'
    iftv = '/.. /iftv.txt'
    CNV_MCD = '/.. /CNV_MCD.txt'
    freec = '/.. /freec.txt'
    file1name = 'freec'
    file2name = 'CNV_MCD'
    file2 = CNV_MCD
    file1 = freec

    with open(file1, 'r') as f2:
        # line = f2.readline()
        for line in f2:
            linestr = line.strip('\n')
            linestrlist = linestr.split('\t')
            file1_start.append(int(linestrlist[1]))
            file1_end.append(int(linestrlist[2]))
            if linestrlist[3] == 'gain':
                file1_type.append("gain")
            else:
                file1_type.append("loss")

    # # sta iftv
    with open(file2, 'r') as f1:
        # line = f1.readline()
        for line in f1:
            linestr = line.strip('\n')
            linestrlist = linestr.split('\t')
            file2_start.append(int(linestrlist[1]))
            file2_end.append(int(linestrlist[2]))
            if linestrlist[3] == 'gain':
                file2_type.append("gain")
            else:
                file2_type.append("loss")

    # sta cnvnator
    # with open(file2, 'r') as f2:
    #     for line in f2:
    #         linestr = line.strip('\n')
    #         linestrlist = linestr.split('\t')
    #         if linestrlist[0] == 'duplication':
    #             file2_type.append("gain")
    #         else:
    #             file2_type.append("loss")
    #         file2_pos = linestrlist[1]
    #         start = re.findall(r'(?<=:)\d*',file2_pos)
    #         end = re.findall(r'(?<=-)\d*',file2_pos)
    #         file2_start.append(int(start[0]))
    #         file2_end.append(int(end[0]))

    count_1 += len(file1_type)  # dpcnv
    count_2 += len(file2_type)  # iftv
    for i in range(len(file1_type)):
        for j in range(len(file2_type)):
            length1 = file1_end[i] - file1_start[i] + 1
            length2 = file2_end[j] - file2_start[j] + 1

            if file2_start[j] <= file1_start[i] <= file2_end[j] and file2_type[j] == file1_type[i]:
                if file1_end[i] <= file2_end[j]:
                    com_length = file1_end[i] - file1_start[i] + 1
                    if com_length / min(length1, length2) >= overlapratio:
                        com_count += 1
                    # print(count)
                elif file1_end[i] >= file2_end[j]:
                    com_length = file2_end[j] - file1_start[i] + 1
                    if com_length / min(length1, length2) >= overlapratio:
                        com_count += 1
                    # print(count)
                break
            elif file2_start[j] >= file1_start[i] and file2_type[j] == file1_type[i]:
                if file2_start[j] <= file1_end[i] <= file2_end[j]:
                    com_length = file1_end[i] - file2_start[j] + 1
                    if com_length / min(length1, length2) >= overlapratio:
                        com_count += 1
                    # print(count)
                elif file1_end[i] >= file2_end[j]:
                    com_length = file2_end[j] - file2_start[j] + 1
                    if com_length / min(length1, length2) >= overlapratio:
                        com_count += 1
                    # print(count)
                break
    # print('count_iftv:',count_1,'\ncount_cnvnator:',count_2,'\ncommon_count:',com_count)
    # print('--------------------'+str(k+1)+'------------------')
print('entry_' + file1name + ': ' + str(count_1))
print('entry_' + file2name + ': ' + str(count_2))
print('entry_' + file1name + '&' + file2name + ': ' + str(com_count))

'''
with open('result_entry','a') as f3:
    f3.write('entry_' + file1name + ': '+ str(count_1) + '\n')
    f3.write('entry_' + file2name + ': '+ str(count_2) + '\n')
    f3.write('entry_' + file1name + '&' + file2name + ': '+ str(com_count) + '\n')
    f3.write('------------------------------------------------------\n')
'''
