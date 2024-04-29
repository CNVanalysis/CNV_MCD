from rpy2 import robjects
from rpy2.rinterface_lib.embedded import RRuntimeError

if __name__ == '__main__':
    robjects.r.source("chordDiagram.R")
    robjects.r.plotchord()

    # CNV_MCD = '/.. /'

    # for k in range(22):
    #     file = CNV_MCD + 'REF_' + str(k + 1) + '.txt'
    #     count = len(open(file, 'rU').readlines())
    #     print(count)
