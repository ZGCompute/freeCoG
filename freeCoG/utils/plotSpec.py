import numpy as np
from matplotlib import pyplot as plt
from matplotlib.mlab import specgram

def plotSpec(dataChan):

    Fs = 512;
    y = dataChan;
    y *= np.hanning(len(y));y = y/np.mean(y);
    yy = np.concatenate((y, ([0] * 10 * len(y))))
    n = len(yy)  # length of the signal
    k = np.arange(n);
    T = n / Fs;
    frq = k / T;  # two sides frequency range
    frq = frq[range(n / 2)];  # one side frequency range

    Y = np.fft.fft(yy) / n;
    Y = Y[range(n / 2)] / max(Y[range(n / 2)]);    
    t = np.linspace(0,len(y)-1,len(y));
    

    # plot the signal and spectrogram
    fig = plt.figure(facecolor="white");
    plt.subplot(3, 1, 1);
    plt.plot(t* 1e3, y, 'r');
    plt.xlabel('Time (miliseconds)');
    plt.ylabel('Amplitude');
    plt.grid();
    
    plt.subplot(3, 1, 2);
    plt.plot(frq[0:600], abs(Y[0:600]), 'k');
    plt.xlabel('Freq (Hz)');
    plt.ylabel('|Y(freq)|');
    plt.grid();
    
    plt.subplot(3, 1, 3);
    plt.specgram(y, NFFT=512, Fs=Fs, noverlap=10);
    plt.ylim(0,255)
    plt.xlim(0,185)
    plt.show();