import os
import numpy as np
import pandas as pd
from tqdm.auto import trange
import matplotlib.pyplot as plt

traj = pd.read_csv('../traj.dat') ### PATH TO TRAJ HERE

n_120 = []
n_60 = []
n_0 = []

for cycle in traj.cycle.unique():
    val_counts = traj[traj.cycle == cycle].orientation.value_counts()
    count_120 = val_counts[120]
    count_60 = val_counts[60]
    count_0 = val_counts[0]
    count_all = count_120 + count_60 + count_0
    n_120.append(count_120 / count_all)
    n_60.append(count_60 / count_all)
    n_0.append(count_0 / count_all)

plt.plot(traj.cycle.unique(), n_120, label='120')
plt.plot(traj.cycle.unique(), n_60, label='60')
plt.plot(traj.cycle.unique(), n_0, label='0')
plt.xticks(rotation=45)
plt.legend(loc=(1.05, 0.35))
plt.show()
