#!/usr/bin/env python

"""
postProc.py:
Module of functions for post processing of retinal 2P data. Other post processing functions for general 
lab usage are in classFun.py
"""

__author__      = "Yannik Bauer"
__copyright__   = ""
__license__     = ""

# Import libs/mode
import sys
import numpy as np
import warnings
import scipy
import scipy.signal # Needs to be loaded explicitly
import matplotlib.pyplot as plt

def normalize(array, mode='z'):
    """
    normalize(array, mode='z'): array = [tPts, trial, cell/cluster]
    Note: to achieve correct scaling, either subtract min and divide by max in two steps;
    or do in one step by also subtracting min in the division term!
    """ 

    if mode is 'z': # z-normalization
        return (array - np.mean(array, axis=0)) / np.std(array, axis=0)

    elif mode is 'r': # range-normalization
        return (array - np.min(array, axis=0)) / (np.max(array, axis=0) - np.min(array, axis=0))
    
    elif mode is 'meanR': # mean range-normalization
        normed = array - np.min(np.mean(array, axis=1), axis=0)
        normed = normed / np.max(np.abs(np.mean(normed, axis=1)), axis=0)
        return normed

    elif mode is 'meanMax': # subtract baseline ((mean) of first 8 samples), then divide by max((abs(mean)) activity)
        normed = array - np.mean(array[0:8], axis=0)
        normed = normed / np.max(np.abs(np.mean(normed, axis=1)), axis=0)
        return normed

    elif mode is 'medMax': # subtract baseline ((median) of first 8 samples), then divide by max((abs(median)) activity)
        normed = array - np.median(array[0:8], axis=0)
        normed = normed / np.max(np.abs(np.median(normed, axis=1)), axis=0)
        return normed

    else:
        warnings.warn("Mode parameter not recognised")

        
def interpNewSRate(trace=None, newSRate=None, duration=None, kind='linear'):
    """
    Returns data resampled to new sampling rate for a given duration by interpolation.
    
    INPUT:
    ------
    trace : ndarray
        original trace.
    newSRate : int
        desired sampling rate.
    duration : int
        desired trace duration.
    kind : str, optional
        interpolation method.
                    
    OUTPUT:
    ------
    traceInterp : ndarray
        interpolated trace.
    """
    # TODO: create analogous function to bring two traces to same sampling rates: interpSameSRate
    
    # Get number of samples in input trace
    nSamplesOld = trace.shape[0] # e.g. nSamples in OGB1 ≈ 249        
    nSamplesNew = newSRate * duration
    
    # Set time vectors based on nSamples and desired duration
    tPointsOld = np.linspace(0, duration, nSamplesOld)
    tPointsNew = np.linspace(0, duration, nSamplesNew)

    # Interpolate trace
    try:
        interpolator = scipy.interpolate.interp1d(tPointsOld, trace[:,:], axis=0, kind=kind) # Interpolate multiple traces
    except:
        interpolator = scipy.interpolate.interp1d(tPointsOld, trace, axis=0, kind=kind) # Interpolate single trace
    traceInterp = interpolator(tPointsNew)
    
#     print("\tResampling data to %d Hz." % newSRate)
    
    return traceInterp


def qi(traces):
    """
    Takes [time, trial, roi] and returns Quality Index (QI).
    """    
    timeVarOfTrialMean = np.var(np.mean(traces, axis=1), axis=0) # cell responsiveness to stimulus
    trialMeanOftimeVar = np.mean(np.var(traces, axis=0), axis=0) # cell consistency / noise variance?
    qi = timeVarOfTrialMean / trialMeanOftimeVar
    
    # Set NaNs to 0
    qi[np.isnan(qi)] = 0
    
#     print('\tExtracting QI.')    
    return qi

def pearsonr2D(A,B):
    """Returns Pearson product-moment correlation coefficient for 2D arrays."""
    # Row-wise mean of input arrays & subtract from input arrays themselves
    A_mA = A - A.mean(1)[:,None]
    B_mB = B - B.mean(1)[:,None]

    # Sum of squares across rows
    ssA = (A_mA**2).sum(1);
    ssB = (B_mB**2).sum(1);

    # Get corr coeff
    r = np.dot(A_mA, B_mB.T) / np.sqrt(np.dot(ssA[:,None], ssB[None]))
    return r

def percentage(a, b):
    """Returns percentage of a/b."""
    return np.round(a / b * 100, 6)


def deconv(trace, fps=None, method=None, smooth=False, norm=True):
    """
    Returns deconvolved neural calcium traces to approximate spikes. 
    INPUT:
    ------
    trace : ndarray
        Array of calcium signals, either single trace [time] or multiple traces [time, trial].
    method : str
        Specify deconvolution method:
        gcamp6fKern : gcamp6 kernel extracted from own data (white noise stimulus)
        ogb1Kern : ogb1 kernel extracted from own data (white noise stimulus)
        artKern : artificial kernel, single- or double-term exponential (uses preset params)
        foopsi : "fast oopsi" (Vogelstein2010)
    smooth : bool
        Smooth trace and/or kernel (uses preset params).
    norm : bool
        Normalize deconvolved trace (range-normalization due to non-negativity of firing rates).

    OUTPUT:
    ------
    d : ndarray
        deconvolved traces
    trace : ndarray
        smoothed trace
        
    TODO: decide which smoothing sigma is reasonable
    """
    
    # Smooth trace with Gaussian filter
    if (smooth and fps is not None):
        sigma = 0.025 # Gaussian sigma (s) 25 ms seems to give best results <> (15 ms in Yaksi2006)
        sigmaN = fps * sigma # sigma in samples = fps * duration (s)
#         print('Smoothing trace with Gaussian kernel. Sigma =', sigma)
#         print(sigmaN)
        trace = scipy.ndimage.filters.gaussian_filter1d(trace, sigmaN, axis=0)
    
    ## Use specified deconvolution method
    # Use kernel extracted from own data
    if (method is 'gcamp6fKern' or method is 'ogb1Kern'):

        # Get kernel extracted from own data
        kern = getCaKern(kernel=method, fps=fps, smooth=True)
                
        # Deconvolve trace(s)
        # To avoid deconvolution clipping, pad trace with repeat of itself (for duration of kernel)
        tracePad = np.concatenate((trace, trace[-kern.shape[0]:-1]))  
        
        if len(trace.shape) == 1: # If single trace...          
            d, rmd = scipy.signal.deconvolve(tracePad, kern)
            
        elif len(trace.shape) > 1: # If multi trace, do for all traces
            d = np.zeros(trace.shape)
            for iTrace in range(trace.shape[1]):
                d[:,iTrace], _ = scipy.signal.deconvolve(tracePad[:,iTrace], kern)
                
        # Normalize deconvolved trace
        if norm:
            d = normalize(d, mode='r') # mode: range-norm due to non-negativity of firing rates
    
    # Fast oopsi non-negative deconvolution toolbox (Vogelstein2010)
    elif method == 'foopsi':

        sys.path.append('../../../code/py-oopsi-master')
        import oopsi
        
        if len(trace.shape) == 1: # If single trace...  
            d, _ = oopsi.fast(trace, dt=1/fps, iter_max=20)
            
        elif len(trace.shape) > 1: # If 2D trace [time,trial], do for all traces
            d = np.zeros(trace.shape)
            for iTrace in range(trace.shape[1]):
                d[:,iTrace], _ = oopsi.fast(trace[:,iTrace], dt=1/fps, iter_max=10)
        
    # Artificial deconvolution kernel
    # TODO: allow artKern1 vs artKern2
    elif method == 'artKern1':
        
        ## Create kernel
        # One-term exponential decay
        kernDur = 1 # Kernel duration (s) - arbitrary
        t = np.linspace(0, kernDur, fps*kernDur)
        tau = 0.205 # gcamp6f decay constant (s), see Chen2013: in vivo, V1, L2/3, 1 AP
        kern = np.exp(-t/tau) # Ca++-kernel
        kern -= min(kern)
        kern /= sum(kern); # normalize Ca++-kernel by AUC        
        
        ## Deconvolve
        # To avoid deconvolution clipping, pad trace with repeat of itself (for duration of kernel)
        tracePad = np.concatenate((trace, trace[-kern.shape[0]:-1]))  
        
        if len(trace.shape) == 1: # If single trace...          
            d, rmd = scipy.signal.deconvolve(tracePad, kern)
            
        elif len(trace.shape) > 1: # If multi trace, do for all traces
            d = np.zeros(trace.shape)
            for iTrace in range(trace.shape[1]):
                d[:,iTrace], _ = scipy.signal.deconvolve(tracePad[:,iTrace], kern)
        
        # Normalize deconvolved trace
        if norm:
            d = normalize(d, mode='r') # mode: range-norm due to non-negativity of firing rates

    elif method == 'artKern2':
        ## Create kernel        
        # Two-term exponential decay
        kernDur = 1 # Kernel duration (s) - arbitrary
        t = np.linspace(0, kernDur, fps*kernDur) # = np.arange(0,1,1/p['newSRate'])
        tau1 = 0.072 # gcamp6f rise constant (s), see Chen2013,  Dana2014: in vivo, V1, L2/3, 1 AP
        tau2 = 0.205 # gcamp6f decay constant (s)
        kern = 1*(1 - np.exp(-t/tau1)) + 1*np.exp(-t/tau2) # Ca++-kernel (weights to both terms set to 1 for now)
        kern -= min(kern)
        kern /= sum(kern); # normalize Ca++-kernel by AUC
        
        ## Deconvolve
        # To avoid deconvolution clipping, pad trace with repeat of itself (for duration of kernel)
        tracePad = np.concatenate((trace, trace[-kern.shape[0]:-1]))  
        
        if len(trace.shape) == 1: # If single trace...          
            d, rmd = scipy.signal.deconvolve(tracePad, kern)
            
        elif len(trace.shape) > 1: # If multi trace, do for all traces
            d = np.zeros(trace.shape)
            for iTrace in range(trace.shape[1]):
                d[:,iTrace], _ = scipy.signal.deconvolve(tracePad[:,iTrace], kern)
        
        # Normalize deconvolved trace
        if norm:
            d = normalize(d, mode='r') # mode: range-norm due to non-negativity of firing rates
    
    else:
        raise ValueError('No valid deconvolution method specified (ensure correct spelling)!')
    
    # Smooth trace with Gaussian filter
    if (smooth and fps is not None):
        #sigma = 0.03 # Gaussian sigma (s) (15 ms in Yaksi2006)
        #  sigmaN = fps * sigma # sigma in samples = fps * duration (s)
#      #   print('Smoothing trace with Gaussian kernel. Sigma =', sigma)
#       #  print(sigmaN)
        d = scipy.ndimage.filters.gaussian_filter1d(d, sigmaN, axis=0)

    return d, trace


def getCaKern(kernel='None', fps=None, smooth=True):
    """
    Get or make calcium kernel extracted from own data, and return processed kernel (smoothed, cut).
    INPUT:
    ------
    kernel : str
        options: {gcamp6fKern, ogb1Kern}
        
    OUTPUT:
    ------
    kern : array
        processed kernel
        
    TODO: 
     - encode loadDir, sRate
     - maybe fit decay constant
    """
    
    import scipy
    import warnings
    
    # Load calcium kernel
    # Extracted by Miro from gcamp6 data by taking cells w good responses to DN stim (n≈20), 
    # then thresholding ca-events, then taking average event, to output a kernel.
    try:
        if kernel == 'gcamp6fKern':
            caKern = scipy.io.loadmat("../../data/2P/proc/gcamp6f_kern.mat")
        elif kernel == 'ogb1Kern':
            caKern = scipy.io.loadmat("../../data/2P/proc/ogb1_kern.mat")
            
        kern = caKern['trace_norm_auc'][0] # Use kernel normalized by AUC
        ts = caKern['ts_trace'][0] # timestamps
        
    except:
        raise FileNotFoundError('No calcium kernel found. Aborting.'\
                                'TODO: in that case consider making it yourself: makeCaKernel().')

    # Resample kernel to sRate of trace (and optionally smooth)
    dur = ts[-1] - ts[0] # kernel duration
    if smooth:
        kern = interpNewSRate(kern, fps, dur, kind='cubic')
    elif not smooth:
        kern = interpNewSRate(kern, fps, dur, kind='linear')

    # Cut kernel: peak to 10th percentile
    # NOTE: full kernel produced ringing artefacts in deconvolution - reason unclear
    kern = kern[np.argmax(kern):-1]
    kern = kern[np.where(kern >= np.percentile(kern,10))]
    
    # Normalize kernel to AUC=1
    kern -= min(kern) # makes decon less noisy
    kern /= sum(kern)

    # plt.plot(kern)
    
    return kern