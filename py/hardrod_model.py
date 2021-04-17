import os
import copy, math, sys, time, shutil, imageio
import numpy as np
import seaborn as sns
import scipy.stats
from tqdm.auto import tqdm
from tqdm.auto import trange
from time import perf_counter
from scipy.interpolate import make_interp_spline
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import matplotlib.animation as animation

from util import *

class Rod:
    def __init__(self, id, loc, orientation, length, width, L, bins):
        self.id = id
        self.loc = loc
        self.orientation = orientation
        self.length = length
        self.width = width
        self.length_by_2 = length / 2
        self.width_by_2 = width / 2
        self.L = L
        self.bins = bins
        self.update()

    def get_endpoints(self):
        v = self.long_axis
        x1 = self.loc[0] - (v[0] / 2 * self.length)
        x2 = self.loc[0] + (v[0] / 2 * self.length)
        y1 = self.loc[1] - (v[1] / 2 * self.length)
        y2 = self.loc[1] + (v[1] / 2 * self.length)
        return x1, y1, x2, y2

    def get_grid_id(self, L, bins):
        self.grid_id = get_grid_coordinate(self.loc[0], self.loc[1], L, bins)

    def get_axes(self):
        self.long_axis = rotate_vector(np.array([0, 1]), self.orientation)
        self.short_axis = rotate_vector(np.array([0, 1]), self.orientation-90)
        self.rotation_matrix = np.zeros((2,2))
        self.rotation_matrix[0,:] = self.long_axis
        self.rotation_matrix[1,:] = self.short_axis

    def get_vertices(self):
        self.vertical_vertices = np.zeros((4,2))
        self.vertical_vertices[0,:] = np.array([-self.length_by_2, self.width_by_2])
        self.vertical_vertices[1,:] = np.array([self.length_by_2, self.width_by_2])
        self.vertical_vertices[2,:] = np.array([self.length_by_2, -self.width_by_2])
        self.vertical_vertices[3,:] = np.array([-self.length_by_2, -self.width_by_2])
        self.vertices = np.zeros((4,2))
        for i in range(4):
            self.vertices[i,:] = self.rotation_matrix@self.vertical_vertices[i,:]+self.loc

    def update(self):
        self.get_grid_id(self.L, self.bins)
        self.get_axes()
        self.get_vertices()

    def __copy__(self):
        return Rod(self.id, self.loc, self.orientation, self.length,
                   self.width, self.L, self.bins)


class GridSpace:
    def __init__(self, grid_neighbors, rod_neighbors):
        self.grid_neighbors = grid_neighbors
        self.rod_neighbors = rod_neighbors

class HardRodModel:
    def __init__(self, n_dim, rod_length=1, aspect_ratio=10, box_length=100,
                 mc_alg='grand_canonical', n_rods=24, cutoff_ratio=1.5,
                 temp=300, params={}, restrict_orientations=False):
        self.n_dim = n_dim
        self.n_rods = n_rods
        self.rod_length = rod_length
        self.rod_width = rod_length / aspect_ratio
        self.aspect_ratio = aspect_ratio
        self.nn_cutoff = self.rod_length * cutoff_ratio
        self.box_length = box_length
        self.n_bins = int(self.box_length / (self.nn_cutoff / 2))
        self.grid_spacing = self.box_length / self.n_bins
        self.grid_bins = np.linspace(self.grid_spacing, self.box_length, self.n_bins)
        # self.grid_spacing = self.grid_bins[1] - self.grid_bins[0]
        self.nn_window = 1
        self.grid = {}
        for i in range(self.n_bins**2):
            x = i % self.n_bins
            y = int(i / self.n_bins)
            grid_neighbors = get_grid_neighbors(x, y, self.n_bins, self.nn_window)
            self.grid[i] = GridSpace(grid_neighbors, [])
        self.mc_alg = mc_alg
        self.cutoff_ratio = cutoff_ratio
        self.temp = temp
        self.kb = 1
        self.beta = 1 / (self.temp * self.kb)
        self.params = params
        self.restrict_orientations = restrict_orientations

        self.model_init(n_rods, n_dim)
        self.rod_history = {0: copy.deepcopy(self.rod_dict)}
        self.history = {'n_rods': [],
                        'swap_success': [],
                        'insert_success': [],
                        'delete_success': [],
                        'rotate_success': [],
                        'translate_success': []}
        self.swap_attempts = 0
        self.swap_successes = 0
        self.insert_attempts = 0
        self.insert_successes = 0
        self.delete_attempts = 0
        self.delete_successes = 0
        self.rotate_attempts = 0
        self.rotate_successes = 0
        self.translate_attempts = 0
        self.translate_successes = 0

    def model_init(self, n_rods, n_dim):
        """
        Randomly initialize rods on 2d surface
        """
        self.rod_dict = {}
        self.available_rod_ids = []
        for i in range(n_rods):
            loc = get_rand_loc(self.box_length, self.n_dim)
            orientation = get_rand_orientation(self.n_dim,
                          restricted=self.restrict_orientations)
            self.rod_dict[i] = Rod(i, loc, orientation, self.rod_length,
                                   self.rod_width, self.box_length, self.grid_bins)
            ### Add rod_id to nearest grid neighbors
            for grid_id in self.grid[self.rod_dict[i].grid_id].grid_neighbors:
                self.grid[grid_id].rod_neighbors.append(i)
        self.available_rod_ids.append(n_rods)

        # re-draw overlapping rods if necessary
        self.remove_overlaps()

    def remove_overlap(self, rod, n_iter=1000):
        for i in range(n_iter):
            no_overlaps = True
            for neighbor_id in self.grid[rod.grid_id].rod_neighbors:
                overlap = check_overlap(rod, self.rod_dict[neighbor_id])
                if overlap:
                    no_overlaps = False
                else:
                    pass
            if no_overlaps:
                break
            else: ### move rod if overlap
                ### remove current rod from neighbor_lists
                for grid_id in self.grid[rod.grid_id].grid_neighbors:
                    rod_neighbors = self.grid[grid_id].rod_neighbors
                    rod_neighbors.pop(rod_neighbors.index(rod.id))
                    self.grid[grid_id].rod_neighbors = rod_neighbors
                ### move and re-orient
                rod.loc = get_rand_loc(self.box_length, self.n_dim)
                rod.orientation = get_rand_orientation(self.n_dim,
                                     restricted=self.restrict_orientations)
                ### add new rod to neighbor lists
                rod.update()
                for grid_id in self.grid[rod.grid_id].grid_neighbors:
                    rod_neighbors = self.grid[grid_id].rod_neighbors
                    rod_neighbors.append(rod.id)
                    self.grid[grid_id].rod_neighbors = rod_neighbors
                self.rod_dict[rod.id] = rod

        # print('Rod {} -> {}'.format(rod['id'], i))
        if i == n_iter - 1:
            print('Could not relocate rod {}'.format(rod.id))

    def remove_overlaps(self):
        for k, rod in tqdm(self.rod_dict.items()):
            self.remove_overlap(rod)

    def mc_swap(self, rod_id, k):
        """
        Attempts to move a single rod to a new location and rotate
        """
        # Move rod to new location and re-calculate nearest neighbors
        old_rod = copy.copy(self.rod_dict[rod_id])
        og_orientation = old_rod.orientation
        new_rod = copy.copy(old_rod)
        self.rod_dict[rod_id] = new_rod
        new_rod.loc = get_rand_loc(self.box_length, self.n_dim)
        new_rod.get_grid_id(self.box_length, self.grid_bins)
        self.rod_dict[rod_id] = old_rod

        w_olds = [1]
        w_news = []
        new_orientations = []

        # generate new weights
        for _ in range(k):
            new_rod.orientation = get_rand_orientation(self.n_dim,
                                  restricted=self.restrict_orientations)
            new_rod.get_axes()
            new_rod.get_vertices()
            new_orientations.append(new_rod.orientation)
            no_overlaps = True
            for neighbor_id in self.grid[new_rod.grid_id].rod_neighbors:
                overlap = check_overlap(new_rod, self.rod_dict[neighbor_id])
                if overlap:
                    no_overlaps = False
                    break
                else:
                    pass
            if no_overlaps:
                w_news.append(1)
            else:
                w_news.append(0)

        # generate old weights
        self.rod_dict[rod_id] = old_rod
        for _ in range(k-1):
            old_rod.orientation = get_rand_orientation(self.n_dim,
                                  restricted=self.restrict_orientations)
            old_rod.get_axes()
            old_rod.get_vertices()
            no_overlaps = True
            for neighbor_id in self.grid[old_rod.grid_id].rod_neighbors:
                overlap = check_overlap(old_rod, self.rod_dict[neighbor_id])
                if overlap:
                    no_overlaps = False
                    break
                else:
                    pass
            if no_overlaps:
                w_olds.append(1)
            else:
                w_olds.append(0)

        w_old = sum(w_olds)
        w_new = sum(w_news)

        # select new orientation to keep
        new_rod.orientation = select_weighted_config(new_orientations, w_news,
                                                     w_new, k)

        acc = w_new / w_old
        if np.random.rand() < acc:
            ### remove current rod from neighbor_lists
            for grid_id in self.grid[self.rod_dict[rod_id].grid_id].grid_neighbors:
                rod_neighbors = self.grid[grid_id].rod_neighbors
                rod_neighbors.pop(rod_neighbors.index(rod_id))
                self.grid[grid_id].rod_neighbors = rod_neighbors
            ### add new rod to neighbor lists
            new_rod.get_axes()
            new_rod.get_vertices()
            self.rod_dict[rod_id] = new_rod
            for grid_id in self.grid[self.rod_dict[rod_id].grid_id].grid_neighbors:
                rod_neighbors = self.grid[grid_id].rod_neighbors
                rod_neighbors.append(rod_id)
                self.grid[grid_id].rod_neighbors = rod_neighbors
            self.swap_successes += 1
        else:
            self.rod_dict[rod_id].orientation = og_orientation
            self.rod_dict[rod_id].get_axes()
            self.rod_dict[rod_id].get_vertices()

        self.swap_attempts += 1

    def mc_insert(self, k):
        # select rod id and pick random location
        rod_id = np.random.choice(self.available_rod_ids, size=1)[0]
        loc = get_rand_loc(self.box_length, self.n_dim)
        orientation = get_rand_orientation(self.n_dim,
                      restricted=self.restrict_orientations)

        new_rod = Rod(rod_id, loc, orientation, self.rod_length, self.rod_width,
                      self.box_length, self.grid_bins)

        w_news = []
        new_orientations = []
        for i in range(k):
            new_orientations.append(new_rod.orientation)
            no_overlaps = True
            for neighbor_id in self.grid[new_rod.grid_id].rod_neighbors:
                overlap = check_overlap(new_rod, self.rod_dict[neighbor_id])
                if overlap:
                    no_overlaps = False
                    break
                else:
                    pass
            if no_overlaps:
                w_news.append(1)
            else:
                w_news.append(0)

            if i == k - 1:
                pass
            else:
                new_rod.orientation = get_rand_orientation(self.n_dim,
                                      restricted=self.restrict_orientations)
                new_rod.get_axes()
                new_rod.get_vertices()

        w_new = sum(w_news)

        # select new orientation to keep
        new_rod.orientation = select_weighted_config(new_orientations, w_news,
                                                     w_new, k)

        # calculate acceptance
        V = self.box_length ** 2
        N = self.n_rods
        mu = self.params['mu']
        acc = (V*np.exp(self.beta*mu)/(N+1))*(w_new / k)
        if np.random.rand() < acc:
            ### add new rod to dictionary
            new_rod.get_axes()
            new_rod.get_vertices()
            self.rod_dict[rod_id] = new_rod
            self.n_rods += 1
            ### remove rod id from available and add to list if no more rod ids left
            self.available_rod_ids.pop(self.available_rod_ids.index(rod_id))
            if len(self.available_rod_ids) == 0:
                self.available_rod_ids.append(len(self.rod_dict.keys()))
            ### add new rod to neighbor lists
            for grid_id in self.grid[self.rod_dict[rod_id].grid_id].grid_neighbors:
                rod_neighbors = self.grid[grid_id].rod_neighbors
                rod_neighbors.append(rod_id)
                self.grid[grid_id].rod_neighbors = rod_neighbors
            self.insert_successes += 1
        else:
            del new_rod

        self.insert_attempts += 1


    def mc_delete(self, rod_id):
        rod = self.rod_dict[rod_id]

        # calculate acceptance
        N = self.n_rods
        V = self.box_length ** 2
        mu = self.params['mu']
        acc = (N/(V*np.exp(self.beta*mu)))
        if np.random.rand() < acc:
            ### delete rod from existence
            for grid_id in self.grid[rod.grid_id].grid_neighbors:
                rod_neighbors = self.grid[grid_id].rod_neighbors
                rod_neighbors.pop(rod_neighbors.index(rod_id))
                self.grid[grid_id].rod_neighbors = rod_neighbors
            self.available_rod_ids.append(rod_id)
            del self.rod_dict[rod_id]
            del rod
            self.n_rods -= 1
            self.delete_successes += 1
        else:
            pass

        self.delete_attempts += 1

    def mc_translate(self, rod_id):
        rod = self.rod_dict[rod_id]
        og_loc = rod.loc
        og_grid_id = rod.grid_id

        max_r = 1
        max_theta = 2*math.pi
        r = np.random.rand() * max_r
        theta = np.random.rand() * max_theta
        x = r * math.sin(theta)
        y = r * math.cos(theta)
        v = np.array([x,y])
        new_loc = og_loc + v
        if new_loc[0] > self.box_length:
            new_loc[0] -= self.box_length
        if new_loc[1] > self.box_length:
            new_loc[1] -= self.box_length
        if new_loc[0] < 0:
            new_loc[0] = self.box_length + new_loc[0]
        if new_loc[1] < 0:
            new_loc[1] = self.box_length + new_loc[1]
        rod.loc = new_loc
        rod.update()
        new_grid_id = rod.grid_id

        no_overlaps = True
        for neighbor_id in self.grid[new_grid_id].rod_neighbors:
            overlap = check_overlap(rod, self.rod_dict[neighbor_id])
            if overlap:
                no_overlaps = False
                break
            else:
                pass

        if no_overlaps:
            if og_grid_id == new_grid_id:
                self.rod_dict[rod_id] = rod
            else:
                ### remove current rod from neighbor lists
                for grid_id in self.grid[og_grid_id].grid_neighbors:
                    rod_neighbors = self.grid[grid_id].rod_neighbors
                    rod_neighbors.pop(rod_neighbors.index(rod_id))
                    self.grid[grid_id].rod_neighbors = rod_neighbors
                ### add new rod to neighbor lists
                self.rod_dict[rod_id] = rod
                for grid_id in self.grid[new_grid_id].grid_neighbors:
                    rod_neighbors = self.grid[grid_id].rod_neighbors
                    rod_neighbors.append(rod_id)
                    self.grid[grid_id].rod_neighbors = rod_neighbors
            self.translate_successes += 1
        else:
            rod.loc = og_loc
            rod.grid_id = og_grid_id
            rod.get_axes()
            rod.get_vertices()
            self.rod_dict[rod_id] = rod

        self.translate_attempts += 1


    def mc_rotate(self, rod_id, k):
        # rotate_stats = {'retrieve_rod_from_storage': None,
        #                 'copy_rod': None,
        #                 'generate_rand_orientation': None,
        #                 'recalc_axes': None,
        #                 'recalc_vertices': None,
        #                 'retrieve_grid_from_storage': None,
        #                 'full_nearest_neighbor_loop': None,
        #                 'number_of_neighbors': None,
        #                 'avg_neighbor_retrieval_time': None,
        #                 'avg_check_overlap_time': None,
        #                 'sum_and_select_config': None,
        #                 'acceptance_criteria': None,
        #                 'recalc_axes_final': None,
        #                 'recalc_vertices_final': None,
        #                 'update_rod_storage': None}
        # start = perf_counter()
        old_rod = self.rod_dict[rod_id]
        # stop = perf_counter()
        og_orientation = old_rod.orientation
        # rotate_stats['retrieve_rod_from_storage'] = stop - start
        # start = perf_counter()
        new_rod = copy.copy(old_rod)
        # stop = perf_counter()
        # rotate_stats['copy_rod'] = stop - start

        # new random orientation
        # start = perf_counter()
        new_rod.orientation = get_rand_orientation(self.n_dim,
                              restricted=self.restrict_orientations)
        # stop = perf_counter()
        # rotate_stats['generate_rand_orientation'] = stop - start
        # start = perf_counter()
        new_rod.get_axes()
        # stop = perf_counter()
        # rotate_stats['recalc_axes'] = stop - start
        # start = perf_counter()
        new_rod.get_vertices()
        # stop = perf_counter()
        # rotate_stats['recalc_vertices'] = stop - start

        w_news = []
        new_orientations = []
        # start_loop = perf_counter()
        for i in range(k):
            new_orientations.append(new_rod.orientation)
            no_overlaps = True
            # start = perf_counter()
            rod_neighbors = self.grid[new_rod.grid_id].rod_neighbors
            # stop = perf_counter()
            # rotate_stats['retrieve_grid_from_storage'] = stop - start
            # rotate_stats['number_of_neighbors'] = len(rod_neighbors)
            # neighbor_retrieval_times = []
            # check_overlap_times = []
            for neighbor_id in rod_neighbors:
                # start = perf_counter()
                neighbor_rod = self.rod_dict[neighbor_id]
                # stop = perf_counter()
                # neighbor_retrieval_times.append(stop - start)
                # start = perf_counter()
                overlap = check_overlap(new_rod, neighbor_rod)
                # stop = perf_counter()
                # check_overlap_times.append(stop - start)
                if overlap:
                    no_overlaps = False
                    break
                else:
                    pass
            if no_overlaps:
                w_news.append(1)
            else:
                w_news.append(0)

            if i == k - 1:
                pass
            else:
                new_rod.orientation = get_rand_orientation(self.n_dim,
                                      restricted=self.restrict_orientations)
                new_rod.get_axes()
                new_rod.get_vertices()
        # stop_loop = perf_counter()
        # rotate_stats['full_nearest_neighbor_loop'] = stop_loop - start_loop
        # rotate_stats['avg_neighbor_retrieval_time'] = np.mean(neighbor_retrieval_times)
        # rotate_stats['avg_check_overlap_time'] = np.mean(check_overlap_times)

        # start = perf_counter()
        w_new = sum(w_news)

        # select new orientation to keep
        new_rod.orientation = select_weighted_config(new_orientations, w_news,
                                                     w_new, k)
        # stop = perf_counter()
        # rotate_stats['sum_and_select_config'] = stop - start

        w_olds = [1]
        for i in range(k-1):
            old_rod.orientation = get_rand_orientation(self.n_dim,
                                  restricted=self.restrict_orientations)
            old_rod.get_axes()
            old_rod.get_vertices()
            no_overlaps = True
            for neighbor_id in self.grid[old_rod.grid_id].rod_neighbors:
                overlap = check_overlap(old_rod, self.rod_dict[neighbor_id])
                if overlap:
                    no_overlaps = False
                else:
                    pass
            if no_overlaps:
                w_olds.append(1)
            else:
                w_olds.append(0)

        w_old = sum(w_olds)

        # start = perf_counter()
        acc = w_new / w_old
        # stop = perf_counter()
        # rotate_stats['acceptance_criteria'] = stop - start
        if np.random.rand() < acc:
            # start = perf_counter()
            new_rod.get_axes()
            # stop = perf_counter()
            # recalc_axes_final = stop - start
            # start = perf_counter()
            new_rod.get_vertices()
            # stop = perf_counter()
            # recalc_vertices_final = stop - start
            # start = perf_counter()
            self.rod_dict[rod_id] = new_rod
            # stop = perf_counter()
            # update_rod_storage = stop - start
            self.rotate_successes += 1
        else:
            old_rod.orientation = og_orientation
            # start = perf_counter()
            old_rod.get_axes()
            # stop = perf_counter()
            # recalc_axes_final = stop - start
            # start = perf_counter()
            old_rod.get_vertices()
            # stop = perf_counter()
            # recalc_vertices_final = stop - start
            # start = perf_counter()
            self.rod_dict[rod_id] = old_rod
            # stop = perf_counter()
            # update_rod_storage = stop - start
        # rotate_stats['recalc_axes_final'] = recalc_axes_final
        # rotate_stats['recalc_vertices_final'] = recalc_vertices_final
        # rotate_stats['update_rod_storage'] = update_rod_storage

        self.rotate_attempts += 1
        # return rotate_stats

    def monte_carlo(self, n_iter, k_rosenbluth=10):
        pass

    def plot_config(self, figsize=(8,8), rod_dict=None, grid=False, highlight_ids=[],
                    style='mica'):
        """
        Visualize the model configuration
        """
        if rod_dict is None:
            rod_dict = self.rod_dict

        fig, ax = plt.subplots(figsize=figsize)
        if style == 'mica':
            ax.set_facecolor("#6F1A06")
        else:
            pass
        ax.set_xticks([])
        ax.set_yticks([])
        for k, rod in rod_dict.items():
            x1, y1, x2, y2 = rod.get_endpoints()
            m = (y2 - y1) / (x2 - x1)
            b = y1 - m * x1
            L = self.box_length
            color = "#AC7C1A"
            if k in highlight_ids:
                if style == 'mica':
                    color = 'white'
                else:
                    color = 'black'
            else:
                if style == 'mica':
                    color = "#AC7C1A"
                else:
                    if rod.orientation < 60:
                        color = '#264653'
                    elif rod.orientation >=60 and rod.orientation < 120:
                        color = '#2a9d8f'
                    elif rod.orientation >= 120 and rod.orientation <=180:
                        color = '#e76f51'
            if x1 > L:
                x_pbc = x1 - L
                print('x1 > L', x1, x2, L, x_pbc, k)
                y_cross = m * L + b
                ax.plot([x_pbc, 0], [y1, y_cross], c=color)
            elif x1 < 0:
                x_pbc = L + x1
                y_cross = b
                ax.plot([x_pbc, L], [y1, y_cross], c=color)
            if y1 > L:
                y_pbc = y1 - L
                x_cross = (L - b) / m
                ax.plot([x1, x_cross], [y_pbc, 0], c=color)
            elif y1 < 0:
                y_pbc = L + y1
                x_cross = -b / m
                ax.plot([x1, x_cross], [y_pbc, L], c=color)
            if x2 > L:
                x_pbc = x2 - L
                y_cross = m * L + b
                ax.plot([0, x_pbc], [y_cross, y2], c=color)
            elif x2 < 0:
                x_pbc = L + x2
                print('x2 < 0', x1, x2, 0, x_pbc, k)
                y_cross = b
                ax.plot([L, x_pbc], [y_cross, y2], c=color)
            if y2 > L:
                y_pbc = y2 - L
                x_cross = (L - b) / m
                ax.plot([x_cross, x2], [0, y_pbc], c=color)
            elif y2 < 0:
                y_pbc = L + y2
                x_cross = -b / m
                ax.plot([x_cross, x2], [L, y_pbc], c=color)

            ## ignore rotation for now
            ax.plot([x1, x2], [y1, y2], c=color)
        ax.set_xlim([0,self.box_length])
        ax.set_ylim([0,self.box_length])
        return plt
