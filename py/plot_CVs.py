import os
import numpy as np
import pandas as pd
from tqdm.auto import trange
import matplotlib.pyplot as plt

colvar = pd.read_csv('../COLVAR.dat') ### PATH TO COLVAR HERE
fig, ax = plt.subplots(1, 2, figsize=(12,6))

ax[0].plot(colvar.cycle, colvar.S)
ax[0].set_xlabel('Cycle')
ax[0].set_title('S', fontweight='bold')
ax[1].plot(colvar.cycle, colvar.density)
ax[1].set_xlabel('Cycle')
ax[1].set_title('Rod Density', fontweight='bold')
plt.show()
