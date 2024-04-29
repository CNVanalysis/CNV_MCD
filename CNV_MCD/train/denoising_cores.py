from scipy.signal import medfilt
import numpy as np
from matplotlib import pyplot as plt
from scipy import signal
import pywt
import scipy.io.wavfile as wav
import obspy
from obspy.signal.detrend import polynomial
from scipy.signal import medfilt, order_filter


def normalize(data):
    data = data.astype('float')
    mx = np.max(data, axis=0).astype(np.float64)
    mn = np.min(data, axis=0).astype(np.float64)
    # Workaround to solve the problem of ZeroDivisionError
    return np.true_divide(data - mn, mx - mn, out=np.zeros_like(data - mn), where=(mx - mn) != 0)


def detrend_score(anomaly_score):
    clf_score = normalize(anomaly_score)
    filter_score = signal.detrend(clf_score, axis=0, type='linear',
                                  bp=0, overwrite_data=False)

    plt.plot(filter_score, c="b")
    plt.plot(clf_score, c="r")
    plt.show()
    return filter_score


def polynomial_filter(anomaly_score):
    clf_score = normalize(anomaly_score)
    filter_score = polynomial(clf_score, order=3, plot=True)
    plt.plot(filter_score, c="b")
    plt.plot(clf_score, c="r")
    plt.show()


def median_filter(anomaly_score, window_size=500): #500
    clf_score1 = normalize(anomaly_score)
    t1 = int(0.2 * window_size)
    t2 = int(0.6 * window_size)
    clf_score2 = medfilt(clf_score1, t1 + 1)
    clf_score3 = medfilt(clf_score2, t2 + 1)
    filter_score = clf_score1 - clf_score3

    # plt.figure(figsize=(20, 4))
    # plt.plot(clf_score1)  
    # plt.show()
    # plt.figure(figsize=(20, 4))
    # plt.plot(clf_score3)  
    # plt.show()
    # plt.figure(figsize=(20, 4))
    # plt.plot(filter_score)  
    # plt.show()

    return filter_score


def remove_baseline_drift(signal, window_size):
    local_mean = np.zeros_like(signal)
    for i in range(len(signal)):
        start = max(1, i - window_size // 2)
        end = min(len(signal), i + window_size // 2 + 1)
        local_mean[i] = np.mean(signal[start:end])
    return signal - local_mean

