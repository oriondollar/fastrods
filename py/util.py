import os
import math, sys, time, shutil, imageio, copy
import numpy as np
import itertools
from itertools import permutations
from time import perf_counter
import matplotlib.pyplot as plt

def get_rand_loc(box_length, n_dim=2):
    """
    Returns a random x/y/z coordinate
    """
    x = np.random.rand() * box_length
    y = np.random.rand() * box_length
    z = np.random.rand() * box_length
    if n_dim == 2:
        loc = np.array([x,y])
    elif n_dim == 3:
        loc = np.array([x,y,z])
    return loc

def get_rand_orientation(n_dim=2, restricted=False):
    """
    Returns a random vector rotated between 0-180 degrees
    """
    if n_dim == 2:
        if restricted:
            rot = np.random.choice([0,60,120], size=1)
        else:
            rot = np.random.random() * 180
    elif n_dim == 3:
        if restricted:
            rot = [get_rand_orientation(n_dim=2, restricted=restricted),
                   get_rand_orientation(n_dim=2, restricted=restricted)]
        else:
            rot = [get_rand_orientation(n_dim=2), get_rand_orientation(n_dim=2)]
    return rot

def get_grid_coordinate(x, y, L, bins):
    x = np.round(x, 3)
    y = np.round(y, 3)
    if x > L:
        x -= L
    elif x < 0:
        x = L + x
    if y > L:
        y -= L
    elif y < 0:
        y = L + y

    n_bins = len(bins)
    found_x = False
    found_y = False
    x_bin = n_bins-1
    y_bin = n_bins-1
    for i in range(n_bins-1):
        bin = bins[i]
        if x < bin and not found_x:
            x_bin = i
            found_x = True
        if y < bin and not found_y:
            y_bin = i
            found_y = True
        if found_x and found_y:
            break
    grid_coord = x_bin + y_bin * n_bins
    return grid_coord

def get_grid_coordinates(rod, dx, L, grid_bins):
    grid_coordinates = []
    x1, y1, x2, y2 = rod.get_endpoints()
    grid_coordinates.append(get_grid_coordinate(x1, y1, L, grid_bins))
    v = np.array([x2, y2]) - np.array([x1, y1])
    for d in np.arange(dx, 1+dx, dx):
        x, y = np.array([x1, y1]) + d*v
        grid_coordinates.append(get_grid_coordinate(x, y, L, grid_bins))
    return grid_coordinates

def get_coord_from_grid_id(grid_id, n_bins, grid_spacing):
    x = grid_id % n_bins
    y = int(grid_id / n_bins)
    x = x * grid_spacing + grid_spacing / 2
    y = y * grid_spacing + grid_spacing / 2
    return x, y

def get_grid_neighbors(x, y, n_bins, nn_window):
    grid_neighbors = []
    x_neighbors = [x]
    y_neighbors = [y]
    for i in range(1,nn_window+1):
        x_up = x + i
        x_down = x - i
        y_up = y + i
        y_down = y - i
        if x_up > n_bins - 1:
            x_up = x_up - n_bins
        if x_down < 0:
            x_down = n_bins + x_down
        if y_up > n_bins - 1:
            y_up = y_up - n_bins
        elif y_down < 0:
            y_down = n_bins + y_down
        x_neighbors.append(x_up)
        x_neighbors.append(x_down)
        y_neighbors.append(y_up)
        y_neighbors.append(y_down)
    for x_neighbor in x_neighbors:
        for y_neighbor in y_neighbors:
            grid_neighbor = x_neighbor + y_neighbor * n_bins
            grid_neighbors.append(grid_neighbor)
    return sorted(grid_neighbors)

def get_nearest_neighbors(rod_id, rod_dict, L, n_dim, nn_cutoff):
    nn_list = []
    rod = rod_dict[rod_id]
    rod_ids = set(rod_dict.keys())
    rod_ids.remove(rod_id)
    for id in rod_ids:
        dist = calc_pbc_dist(rod.loc, rod_dict[id].loc, L, n_dim)
        if dist <= nn_cutoff:
            nn_list.append(id)
        else:
            pass
    return nn_list

def check_overlap_spheres(rod1, rod2):
    overlap = False
    L = rod1.box_length
    d = rod1.sphere_diameter
    for i in range(rod1.sphere_locs.shape[0]):
        for j in range(rod2.sphere_locs.shape[0]):
            loc1 = rod1.sphere_locs[i,:]
            loc2 = rod2.sphere_locs[j,:]
            dist = calc_pbc_dist(loc1, loc2, L)
            if dist < d:
                overlap = True
                return overlap
            else:
                pass
    return overlap

def check_overlap_grids(rod1, rod2):
    if rod1.id == rod2.id:
        return False
    grids1 = rod1.grids_occupied
    grids2 = rod2.grids_occupied
    overlap = False
    for grid_space in grids1:
        if grid_space in grids2:
            overlap = True
            break
    return overlap

def check_overlap(rod1, rod2):
    if rod1.id == rod2.id:
        return False
    move_rod = False
    r = rod1.loc - rod2.loc
    for d in range(r.shape[0]):
        f = round(r[d]/rod1.L)
        if round(f) != 0:
            r[d] = r[d] - rod1.L*f
            move_rod = True
    if move_rod:
        rod2 = copy.copy(rod2)
        rod2.loc = rod1.loc - r
        rod2.get_axes()
        rod2.get_vertices()
    projections = np.zeros((2,4))
    axes = [rod1.long_axis, rod1.short_axis, rod2.long_axis, rod2.short_axis]
    rod_vertices = np.stack([rod1.vertices, rod2.vertices])
    overlap = True
    for i in range(4):
        axis = axes[i]
        projections = rod_vertices@axis

        min_proj_1, min_proj_2 = np.amin(projections, axis=1)
        max_proj_1, max_proj_2 = np.amax(projections, axis=1)

        if (min_proj_1 > max_proj_2) or (max_proj_1 < min_proj_2):
            overlap = False
            return overlap
        else:
            pass
    return overlap

def rotate_vector(v, rot):
    rot = -rot
    rad = rot * math.pi / 180
    cos_rad = math.cos(rad)
    sin_rad = math.sin(rad)
    x_new = v[0]*cos_rad - v[1]*sin_rad
    y_new = v[0]*sin_rad - v[1]*cos_rad
    v = np.array([x_new, y_new])
    return v

def calc_pbc_dist(loc1, loc2, L, n_dim=2):
    if n_dim == 2:
        h = np.array([[L, 0], [0, L]])
    elif n_dim == 3:
        h = np.array([[L, 0, 0], [0, L, 0], [0, 0, L]])
    h_inv = np.linalg.inv(h)
    r = loc1 - loc2
    s = h_inv@r
    move_rod = False
    for f in s:
        if round(f) != 0:
            move_rod = True
    s -= np.round(s)
    r = h@s
    dist = np.sqrt(r@r)
    return dist, move_rod

###############################################################################
####################### MONTE CARLO FUNCTIONS #################################
###############################################################################

def select_weighted_config(config_list, weights, weights_sum, k):
    w_rand = np.random.rand() * weights_sum
    running_w = 0
    idx_sel = k - 1
    for i in range(k):
        running_w += weights[i]
        if running_w > w_rand:
            idx_sel = i
            break
        else:
            pass
    config = config_list[idx_sel]
    return config

###############################################################################
##################### VISUALIZATION FUNCTIONS #################################
###############################################################################

def plot_cycle(cycle, traj, hrm, style='colors'):
    traj = traj[traj.cycle == cycle]

    hrm.rod_dict = {}
    for i, row in traj.iterrows():
        x = row.x
        y = row.y
        orientation = row.orientation
        hrm.rod_dict[i] = Rod(i, np.array([x,y]), orientation, hrm.rod_length,
                              hrm.rod_width, hrm.box_length, hrm.grid_bins)

    config_plot = hrm.plot_config(style=style)
    config_plot.title('Cycle {}'.format(cycle), fontweight='bold')
    return config_plot
