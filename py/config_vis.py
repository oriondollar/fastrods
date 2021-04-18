import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from hardrod_model import *
from util import *

# read config file
config = {'n_dim': None,
          'rod_length': None,
          'aspect_ratio': None,
          'box_length': None,
          'mu': None,
          'n_cycle': None}
with open('../config.dat', 'r') as f: ### PATH TO CONFIG HERE
    for line in f:
        key, val = line.split(':')
        if key in config.keys():
            config[key] = int(val)

n_rods = 10
hrm = HardRodModel(n_dim=config['n_dim'], rod_length=config['rod_length'],
                   aspect_ratio=config['aspect_ratio'], n_rods=n_rods,
                   box_length=config['box_length'], params=config, restrict_orientations=True)

traj = pd.read_csv('../traj.dat') ### PATH TO TRAJ HERE
final_config = plot_cycle(traj.cycle.unique()[-1], traj, hrm, style='colors')
final_config.show()
