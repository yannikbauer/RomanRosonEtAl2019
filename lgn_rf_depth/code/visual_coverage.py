import os
import argparse

import scipy.io as sio
import matplotlib.pyplot as plt
import numpy as np

import datajoint as dj
from djd.util import mixin
import djd.util as util

from RF import get_RF_info, channel_finder, modify_chosen_chs,\
               calculate_accuracy, position_to_ang, visual_coverage
from data_joint_helpers import get_exp_units, translate_mouse_name
from json_tools import read_from_json
from data_joint_helpers import translate_mouse_name


def one_exp_visual_coverage_from_mat(file_name, name_version='new', cell_name=None,
                                     monitor_distance=25.0, monitor_elevation_offset=0.0,
                                     threshold=0, debug=False):
    """
    for one experiment
    ARGS:
        file_name: a "mat" file with RF map
                   example, "data/my_exp_RF_maps.mat"
        cell_name: name of the cell including the 2D RF map
                   example, 'IntRFByCh' or 'RF' cell
                   default is 'IntRFByCh'
    Returns:
        azimuth, elevation: angles in degrees
    """
    
    if cell_name is None:
        cell_name = 'IntRFByCh'
    RF_info = sio.loadmat(file_name)
    RF = RF_info[str(cell_name)]
    
    stds, means, peaks_x, peaks_y, scores, dims = get_RF_info(RF, vis_flag=debug, mat_flag=True)
    """
    calculate accuracy :
        If you already know where RF channels are, you can calculate accuracy
        for this purpose use the function like this
        chs, chy = channel_finder(scores, peaks_y,
                                  ground_truth=[17,18,19,20,21,22,23,24,25,26,27,28],
                                  calculate_acc=True, print_flag=True,
                                  accuracy_method='precision-recall', threshold=0.4)
        Other wise, you can just trust the function O_o
    Threshold:
        you can visualise RF maps from last function and check scores
        to get intuition about reasonable threshold
    """
    file_name= file_name.split('/')[-1]
    key = translate_mouse_name(file_name, name_version='old')
    channels_scores, channels_yslope = channel_finder(scores, peaks_y, threshold=threshold)
    units_chs = get_exp_units(key)
    chs = modify_chosen_chs(channels_scores, channels_yslope, units_chs)
    # You can calculate accuracy again if you want
    # calculate_accuracy(chs, ground_truth, num_chs, method='precision-recall')
    azimuth, elevation, _, _ = position_to_ang(key, peaks_x, peaks_y, dims,
                                    {'distance':monitor_distance,
                                     'elevation_offset':monitor_elevation_offset},
                                     chs, restrict_to_chosen_chs=True)
    return azimuth, elevation, chs


def one_exp_visual_coverage_from_key(key, RF, monitor_distance=25.0, monitor_elevation_offset=0.0, threshold=0, debug=False):
    """
    for one experiment
    ARGS:
        key : {'m','s','e'}
        RF: 2D RF map
    Returns:
        azimuth, elevation: angles in degrees
    """
    stds, means, peaks_x, peaks_y, scores, dims = get_RF_info(RF, vis_flag=debug)
    channels_scores, channels_yslope = channel_finder(scores, peaks_y, threshold=threshold)
    units_chs = get_exp_units(key)
    chs = modify_chosen_chs(channels_scores, channels_yslope, units_chs)
    # calculate_accuracy(chs, ground_truth, num_chs, method='precision-recall')
    azimuth, elevation, _, _ = position_to_ang(key, peaks_x, peaks_y, dims,
                                    {'distance':monitor_distance,
                                     'elevation_offset':monitor_elevation_offset},
                                     chs, restrict_to_chosen_chs=True)
    return azimuth, elevation, chs


"""
if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("database_name", help="database_name",
                        type=str)
    args = parser.parse_args()

    For mat data
    data_dir = 'sample_data/'
    # For many experiments:
    azimuth= []
    elevation = []
    for file_name in os.listdir(data_dir):
        print(file_name)
        info = translate_mouse_name(file_name)
        azi, elev = one_exp_visual_coverage_from_mat(data_dir + file_name)
        azimuth.extend(azi)
        elevation.extend(elev)
    visual_coverage(azimuth, elevation)
"""
    