import numpy as np

import datajoint as dj
from djd.util import mixin
import djd.util as util

from json_tools import read_from_json

import os


conf_path = os.path.realpath(__file__).rsplit('/',1)[0]
conf = read_from_json(conf_path + '/dj_local_conf.json')
dj.config.load(conf_path + '/dj_local_conf.json')
conn = dj.conn()
which_database = conf["DBNAME"]

schema = dj.schema(which_database, connection=conn, create_schema=True, create_tables=True)
print("Connected to database %r as %r" % (schema.database, conn.get_user()))
print("only needed tables for this script are imported")

if which_database == "dj_data":  # Tubingen
    schema.spawn_missing_classes()
    Mice = mixin(Mice)
    Series = mixin(Series)
    Experiments = mixin(Experiments)
    ClusterInfo = mixin(ClusterInfo)
elif which_database == "dj_d":  # Munich
    from djd.mouse import Mouse
    from djd.series import Series
    from djd.unit import Unit
    from djd.stimulus import Stimulus


def translate_mouse_name(name, name_version='old'):
    """
    ARGs:
        name: a mouse name like "BL6_0289_05_(1)"
    returns:
        a dictionary like 
        {'id':BL6_0289, 
        'exp':[1],
        'series':5}
    """
    name_split = name.split("_")
    if name_version == 'old':
        mouse_id = name_split[0] + '_' + name_split[1]
        series_num = int(name_split[2])
        if len(name_split[3].split(',')) > 1:
            exp_id = name_split[3][1:-1].split(',')
            exp_id = [int(s) for s in exp_id]
        else:
            exp_id = [int(name_split[3][1:-1])]
    elif name_version == 'new':
        mouse_id = name_split[0] + '_' + name_split[1] + '_' + name_split[2]
        series_num = int(name_split[3])
        if len(name_split[4].split(',')) > 1:
            exp_id = name_split[4][1:-1].split(',')
            exp_id = [int(s) for s in exp_id]
        else:
            exp_id = [int(name_split[4][1:-1])]
    return {'m': mouse_id, 'e': exp_id, 's': series_num}


def get_exp_units(key):
    """
    ARGS:
        key: {'m','s','e'}
    returns:
        units_chs : list of channels associated to
                    saved units for the specific experiment
    """
    global which_database
    if which_database == "dj_d":  # Munich
        unit_chs = np.unique((Unit.Properties() & key).fetch('u_maxchan'))
    elif which_database == "dj_data":  # Tubingen
        mouse_counter = (Mice() & {'mouse_id':key['m']}).fetch('mouse_counter')
        unit_chs = np.unique((ClusterInfo() & {'mouse_counter': int(mouse_counter),
                                               'series_num': int(key['s'])}).fetch('chan_num'))        
    return np.unique(unit_chs)


def get_exp_monitor_info(key):
    """
    ARGS:
        key: {'m','s','e'}
    returns:
        mon_angle : azimuth angle of monitor relative to mouse
        mon_elev : elevation angle of monitor relative to mouse
    """
    global which_database
    if which_database == "dj_d":  # Munich
        mon_angle = (Series.Experiment() & key).fetch('e_displayazim')
        mon_elev = (Series.Experiment() & key).fetch('e_displayelev')   
    elif which_database == "dj_data":  # Tubingen
        mouse_counter = (Mice() & {'mouse_id':key['m']}).fetch('mouse_counter')
        mon_angle = (Experiments() & {'mouse_counter': int(mouse_counter),
                                      'series_num': int(key['s']),
                                      'exp_num': int(key['e'][0])}).fetch('exp_monitorangle')
        mon_elev = (Experiments() & {'mouse_counter': int(mouse_counter),
                                     'series_num': int(key['s']),
                                     'exp_num': int(key['e'][0])}).fetch('exp_monitorelevation')
    return mon_angle, mon_elev


def get_exp_gratings_info(key):
    """
    ARGS:
        key: {'m','s','e'}
    returns:
        grating_x, grating_y : x & position of grating stimulus relative to screen
    """
    global which_database
    if which_database == "dj_d":  # Munich
        grating_x = (Stimulus.GratingCond() & key).fetch('grat_x_position')
        grating_y = (Stimulus.GratingCond() & key).fetch('grat_y_position')
    elif which_database == "dj_data":  # Tubingen
        mouse_counter = (Mice() & {'mouse_id':key['m']}).fetch('mouse_counter')
        grating_x = (GratingConditions() & {'mouse_counter': int(mouse_counter),
                                            'series_num': int(key['s']),
                                            'exp_num': int(key['e'][0])}).fetch('grat_x_position')
        grating_y = (GratingConditions() & {'mouse_counter': int(mouse_counter),
                                            'series_num': int(key['s']),
                                            'exp_num': int(key['e'][0])}).fetch('grat_y_position')
    grating_x = np.unique([float(x) for x in grating_x])
    grating_y = np.unique([float(y) for y in grating_y])
    return grating_x, grating_y