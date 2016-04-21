# get_CCEP_spec.py

# calculate and plot the event related spectral perturbations
# evoked by cortical and sub-cortical stimulation


import numpy as np
import scipy
import scipy.io
import os

def get_CCEP_spec(subj,stim_block):
    
    # load in epocked data
    data = scipy.io.loadmat(stim_block);
    trigg_chan = data.get('trigg_chan');
    data = data.get('avg_data');
    
    # define srate and signal to decompose
    Fs = 512;
    y = data[1,:];
    y *= np.hanning(len(y));
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
    subplot(3, 1, 1);
    plot(t * 1e3, y, 'r');
    xlabel('Time (miliseconds)');
    ylabel('Amplitude');
    grid();
    
    subplot(3, 1, 2);
    plot(frq[0:256], abs(Y[0:256]), 'k');
    xlabel('Freq (Hz)');
    ylabel('|Y(freq)|');
    grid();
    
    subplot(3, 1, 3);
    Pxx, freqs, bins, im = specgram(y, NFFT=512, Fs=Fs, noverlap=10);
    show();