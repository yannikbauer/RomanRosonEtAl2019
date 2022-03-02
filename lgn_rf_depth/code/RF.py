import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
import numpy as np

import colors as colors
from data_joint_helpers import get_exp_monitor_info,\
                               get_exp_gratings_info


def get_RF_info(RF, vis_flag=False, mat_flag=False):
    """
    For one experiment
    Args:
        RF: receptive field map for each channel 
            array of the shape (1, number of channels)
        vis_flag: if True, RF maps and output figures will be generated, False by default
    
    Returns:
        stds: list with mean value for the RF map for each channel
        means: list with std for the RF map for each channel
        peaks of x value: list of positions of peaks in x values for the RF map for each channel
        peaks of y value: list of positions of peaks in y values for the RF map for each channel
        scores: list of scores given to every channel in the given experiment
        ch_shape: shape of every RF map 
    """
    means = []
    stds = []
    scores=[]
    peaks_y = []
    peaks_x = []
    if vis_flag:
        fig2, axs2 = plt.subplots(4,10, figsize=(20, 8))
        axs2 = axs2.ravel()
        j=0
    
    if mat_flag:
        chs = RF.reshape(-1)
    else:
        chs = RF

    for ch in chs:
        info_y = np.sum(ch,1) # collapse in y
        info_x = np.sum(ch,0) # collapse in x

        means.append(np.mean(ch))
        stds.append(np.std(ch))
        # calculate score for channel
        scores.append((np.amax(info_y) - np.amin(info_y))/np.amax(info_y))
        peaks_y.append(np.argmax(info_y))
        peaks_x.append(np.argmax(info_x))
        if vis_flag:
            axs2[j].imshow(ch)#, vmin=0, vmax=0.5)
            axs2[j].set_title(j)
            j+=1
        
    if vis_flag:
        axs2[33].plot(np.array(stds), 'r')
        axs2[33].set_xlabel('channel num')
        axs2[33].set_title('std')
        axs2[32].plot(np.array(means), 'g')
        axs2[32].set_xlabel('channel num')
        axs2[32].set_title('mean')
        
        axs2[39].plot(np.array(scores), 'c')
        axs2[39].set_xlabel('channel num')
        axs2[39].set_title('scores')

        axs2[38].plot(savgol_filter(np.array(peaks_x), 15, 3), savgol_filter(np.array(peaks_y), 15, 3))
        axs2[38].set_xlabel('x peaks')
        axs2[38].set_ylabel('y peaks')
        axs2[38].set_title('peaks')
        axs2[38].set_xlim((0,57))
        axs2[38].set_ylim((0,57))
        axs2[37].plot(savgol_filter(np.array(peaks_y), 15, 3), 'r')
        axs2[37].set_xlabel('channel num')
        axs2[37].set_title('peaks_y')
        axs2[37].set_ylim((0,53))
        axs2[36].plot(savgol_filter(np.array(peaks_x), 15, 3), 'g')
        axs2[36].set_xlabel('channel num')
        axs2[36].set_title('peaks_x')
        axs2[36].set_ylim((0,57))
        
        plt.show()
    return stds, means, peaks_x, peaks_y, scores, chs[0].shape


def channel_finder(scores, peaks_y, ground_truth=None, threshold=0,
                    name='This experiment', print_flag=False,
                    calculate_acc=False, accuracy_method='num_chs'):
    """
    given one experiment and assigned scores for each channel
    decide which channels have good quality RFs and label them
    ARGS:
        scores: (list)assigned to each channel
        peaks_y: (list)positions of peaks in y values
                  for the RF map for each channel
        ground_truth: (list) of channels
        threshold: assigned threshold between 0 and 1
        name: optional (string) of experiment name for debug printing if needed
        print_flag: (boolean) prints experiment name if given, scores and chosen channels
                     for debugging
        calculate_acc: (boolean) if set to True detected channels will be compared to given
                        ground truth, accuracy method can be chosen
        accuracy_method: (string) for accuracy calculation method
                         either num_chs
                         or false negative
                         or flase positive
                         or precision-recall

    RETURNS:
        channels_scores: (list)chosen chs according to scores combared to threshold
        channels_yslope: (list)chosen channels according to positive y_slope as indicator
                         of RF elevation shift with depth
    """
    if ground_truth is None:
        ground_truth = []

    # choose channels with score higher than threshold as threshold
    channels_scores = [i for i, s in enumerate(scores) if s > threshold]

    # choose channels with increasing y-position of peak activity
    dist = np.diff(np.array(savgol_filter(peaks_y, 15, 3)))
    channels_yslope = [int(x+1) for x in np.where(dist>0)[0]]

    if calculate_acc:
        acc_scores = calculate_accuracy(channels_scores, ground_truth, len(scores), accuracy_method)
        acc_yslope = calculate_accuracy(channels_yslope, ground_truth, len(scores), accuracy_method)

    if print_flag:
        print(name)
        print(colors.BLUE + 'Detected by score')
        print(channels_scores)
        print(colors.GREEN + 'Detected by slope')
        print(channels_yslope)
        if ground_truth:
            print(colors.RED + 'Manually labeled')
            print(ground_truth)
        print(colors.RESET)
        if calculate_acc:
            print('==========================')
            print(colors.BLUE + 'Detected by score')
            print(acc_scores)
            print(colors.GREEN + 'Detected by slope')
            print(acc_yslope)
            print(colors.RESET)
    return channels_scores, channels_yslope


def calculate_accuracy(chosen_chs, ground_truth_chs, num_chs, method='num_chs'):
    """
    ARGS:
        method: either num_chs
                or false negative
                or flase positive
                or precision-recall
    """
    cost = 0
    if method == 'num_chs':
        cost = abs(len(chosen_chs) - len(ground_truth_chs))/num_chs
    false_pos = [ch for ch in chosen_chs if ch not in ground_truth_chs]
    false_neg = [ch for ch in ground_truth_chs if ch not in chosen_chs]
    if method == 'false_pos':
        cost = len(false_pos)/num_chs
    if method == 'false-neg':
        cost = len(false_neg)/num_chs
    if method == 'precision-recall':
        p = (num_chs - len(false_pos)) / (num_chs)
        r = (num_chs - len(false_pos)) / (num_chs - len(false_pos) + len(false_neg))
        F1 = (2 * p * r) / (p + r)
        cost = 1 - F1
    return cost


def modify_chosen_chs(channels_scores, channels_yslope, units_chs):
    """
    For one experiment 
    ARGS:
        channels_scores -- (list)the chosen channels from the original 32 chs according to given score for variance of peak in y axis
        channels_yslope -- (list)the chosen channels from the original 32 chs according to the slope of peaks in y axis
        units_chs -- (List)the channels that correspond to recorded units
    Returns:
        chs -- (list)modified choice of channelss
    """
    detected_chs = channels_scores
    # correct indices for channels of recorded units
    units_chs = [int(ch - 1) for ch in units_chs]
    # reject outliers
    dist = np.diff(np.array(units_chs))
    if np.any(dist > 6):
        units_chs = units_chs[np.where(dist>6)[0][0]+1:]

    # reject outliers in detected channels
    if len(channels_scores) >7:
        dist = np.diff(savgol_filter(np.array(channels_scores), 7, 3))
    else:
        dist = np.diff(np.array(channels_scores))
    if np.any(dist > 3):
        if np.where(dist>3)[0][-1] > len(channels_scores)/2:
            detected_channels = channels_scores[:np.where(dist>3)[0][-1]+1]
        else:
            detected_channels = channels_scores[np.where(dist>3)[0][-1]+1:]
    else:
        detected_channels = channels_scores
    
    # fill in gaps and correct boundaries
    dist = np.diff(np.array(detected_channels))
    chs = []
    for k, x in enumerate(detected_channels):
        chs.append(int(x))
        if k < len(detected_channels) - 1 and detected_channels[k+1] - x > 1: # fill in gaps
            for l in range(1, detected_channels[k+1] - x):
                if x+l in channels_yslope or x+l in units_chs:
                    chs.append(x+l)
        elif k >= len(detected_channels): # correct boundaries
            for l in range(1,5):
                if x+l in channels_yslope or x+l in units_chs: 
                    chs.append(x+l)
    residual = [ch for ch in chs if ch not in units_chs]

    return chs


def position_to_ang(key, xpeaks, ypeaks, dims, monitor_data, chosen_chs, restrict_to_chosen_chs=True):
    """
    converts every RF position to angles in degrees relative to mouse eyes
    ARGs:
        name: a mouse name like "BL6_0289_05_(1)"
        xpeaks: positions of RFs on x-axis
        ypeka: positions of RFs on y-axis
        dims: dimensions of the RF map ex: (57,57)
        monitor_data: dict with monitor info {'distance', 'elevation offset'}
                      offset in case center of the screen is not aligned with mouse eyes
    returns:
        RF_h_ang: list of RF azimuth angles in degrees relative to mouse eyes
        RF_v_ang: list of RF elevation angles in degrees relative to mouse eyes
        x_ang: list of stimulus azimuth angles in degrees relative to mouse eyes
        y_ang: : list of stimulus angles in degrees relative to mouse eyes
    """
    print(colors.GREEN + "chosen channels for ")
    print(key)
    print(colors.RESET)
    print(chosen_chs)

    # get metadata
    mon_azi_angle, mon_elev = get_exp_monitor_info(key)
    # get gratings info
    grating_x, grating_y = get_exp_gratings_info(key)
    mon_elev_ang = np.degrees(np.arctan(float(mon_elev) / monitor_data['distance']))
    if monitor_data['elevation_offset'] is not None:
        mon_elev_ang += monitor_data['elevation_offset'] 
    x_ang = mon_azi_angle + grating_x
    y_ang = mon_elev_ang + grating_y

    x_range = x_ang[-1] - x_ang[0]
    y_range = y_ang[-1] - y_ang[0]

    RF_h_ang = [(float(x) * x_range / dims[0]) + x_ang[0] for x in xpeaks]
    RF_v_ang = [(float(y) * y_range / dims[1]) + y_ang[0] for y in ypeaks]
    if restrict_to_chosen_chs:
        RF_h_ang = [RF_h_ang[ch] for ch in chosen_chs]
        RF_v_ang = [RF_v_ang[ch] for ch in chosen_chs]

    return RF_h_ang, RF_v_ang, x_ang, y_ang


def visual_coverage(azimuth, elevation, bins=None, save_fig=False):
    """
    method for displaying histogram of azimuth and elevation angles
    ARGs:
        azimuth: (list) of angles
        elevation: (list) of angles
        bins: ([list of bins for azimuth], 
               [list of biins for elevation])
               if not set default will be used
               bins=(np.arange(min(azimuth) - 5, max(azimuth) + 5, 5).tolist(),
                     np.arange(min(elevation) - 5, max(elevation) + 5, 5).tolist())
        save_fig: (boolean) will save the histogram as fig in the directory from which script was started 
                    under the name of the function
                    otherwise, figure will be shown
    """
    if bins is None:
        bins=(np.arange(min(azimuth) - 5, max(azimuth) + 5, 5).tolist(),
              np.arange(min(elevation) - 5, max(elevation) + 5, 5).tolist())
    H = np.histogram2d(azimuth, elevation, bins)
    fig = plt.figure()
    im = plt.imshow(H[0].T, interpolation='none', origin='low', extent=[H[1][0], H[1][-1], H[2][0], H[2][-1]])
    plt.xlabel('Azimuth (deg)')
    plt.ylabel('Elevation (deg)')
    # plt.xticks([0,20,40,60,80,100])
    # plt.yticks([-20,0,20,40])
    cb = plt.colorbar(im, fraction=0.04, pad=0.04)#, ticks=[0,10,20,30])
    cb.ax.set_ylabel('number of channels')

    if save_fig:
        plt.savefig('visual_coverage.png')
    else:
        plt.show(block=True)