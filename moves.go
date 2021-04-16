package main

import (
	"math"
	"math/rand"
)

func Swap(rod *Rod, grid []*GridSpace, config *Config, rods []*Rod) {
	// copy rod structs to use for rosenbluth trials
	og_rod := RodDeepCopy(rod)
	og_grid_id := og_rod.grid_id
	new_rod := RodDeepCopy(rod)

	// move rod to new location
	GetRandLoc(config.n_dim, config.box_length, new_rod)
	GetGridID(config.box_length, config.n_bins, &config.grid_bins, new_rod)
	new_grid_id := new_rod.grid_id

	// run k rosenbluth trials to generate swap weights (new location - new orientations)
	GetRandOrientation(config.restrict_orientations, new_rod)
	GetAxes(new_rod)
	GetVertices(config.n_dim, config.n_vertices, new_rod)
	new_weights := make([]float64, config.k)
	var new_weight_sum float64 = 0
	new_orientations := make([]float64, config.k)
	for i := 0; i < config.k; i++ {
		new_orientations[i] = new_rod.orientation
		no_overlaps := CheckNeighborOverlaps(new_rod, grid, rods, config)
		if no_overlaps {
			new_weights[i] = 1
			new_weight_sum += 1
		} else {
			new_weights[i] = 0
		}

		if i != (config.k - 1) {
			GetRandOrientation(config.restrict_orientations, new_rod)
			GetAxes(new_rod)
			GetVertices(config.n_dim, config.n_vertices, new_rod)
		}
	}

	// select new configuration
	new_rod.orientation = SelectWeightedConfig(new_orientations, new_weights, new_weight_sum, config.k)

	// run k-1 rosenbluth trials to generate weights for original configuration
	var og_weight_sum float64 = 1
	for i := 0; i < (config.k - 1); i++ {
		GetRandOrientation(config.restrict_orientations, og_rod)
		GetAxes(og_rod)
		GetVertices(config.n_dim, config.n_vertices, og_rod)
		no_overlaps := CheckNeighborOverlaps(og_rod, grid, rods, config)
		if no_overlaps {
			og_weight_sum += 1
		}
	}

	// calculate acceptance probability for swap
	acc := new_weight_sum / og_weight_sum
	if rand.Float64() < acc {
		// move accepted
		GetAxes(new_rod)
		GetVertices(config.n_dim, config.n_vertices, new_rod)
		if og_grid_id == new_rod.grid_id {
			rods[rod.id] = new_rod
		} else {
			new_rod.grid_id = og_grid_id // temporarily store old grid_id in rod to update old neighbor lists
			RemFromNeighborLists(new_rod, grid)
			new_rod.grid_id = new_grid_id // restore new grid_id to add to new neighbor lists
			AddToNeighborLists(new_rod, grid)
			rods[rod.id] = new_rod
		}
		config.swap_successes++
	}
	config.swap_attempts++
}

func Translate(rod *Rod, grid []*GridSpace, config *Config, rods []*Rod) {
	// store original rod location and grid id
	og_loc := rod.loc
	og_grid_id := rod.grid_id

	// get new location within 1 distance away from current rod COM
	var max_r float64 = 1
	max_theta := 2 * math.Pi
	r := rand.Float64() * max_r
	theta := rand.Float64() * max_theta
	v := [2]float64{r * math.Sin(theta), r * math.Cos(theta)}
	new_loc := make([]float64, config.n_dim)
	new_loc[0] = og_loc[0] + v[0]
	new_loc[1] = og_loc[1] + v[1]
	if new_loc[0] > config.box_length {
		new_loc[0] -= config.box_length
	}
	if new_loc[1] > config.box_length {
		new_loc[1] -= config.box_length
	}
	if new_loc[0] < 0 {
		new_loc[0] += config.box_length
	}
	if new_loc[1] < 0 {
		new_loc[1] += config.box_length
	}
	rod.loc = new_loc
	GetGridID(config.box_length, config.n_bins, &config.grid_bins, rod)
	GetAxes(rod)
	GetVertices(config.n_dim, config.n_vertices, rod)
	new_grid_id := rod.grid_id

	// check new rod location for overlaps
	no_overlaps := CheckNeighborOverlaps(rod, grid, rods, config)

	// automatically accept move if there are no overlaps
	if no_overlaps {
		// move accepted
		if og_grid_id == rod.grid_id {
			rods[rod.id] = rod
		} else {
			rod.grid_id = og_grid_id // temporarily store old grid_id in rod to update old neighbor lists
			RemFromNeighborLists(rod, grid)
			rod.grid_id = new_grid_id // restore new grid_id to add to new neighbor lists
			AddToNeighborLists(rod, grid)
		}
		config.translation_successes++
	} else {
		// move rejected
		rod.loc = og_loc
		rod.grid_id = og_grid_id
		GetAxes(rod)
		GetVertices(config.n_dim, config.n_vertices, rod)
	}
	config.translation_attempts++
}

func Rotate(rod *Rod, grid []*GridSpace, config *Config, rods []*Rod) {
	// copy rod structs to use for rosenbluth trials
	og_rod := RodDeepCopy(rod)
	new_rod := RodDeepCopy(rod)

	// set rosenbluth k to 1 for rotation
	k := 1

	// run k rosenbluth trials to generate rotation weights (same location - new orientations)
	GetRandOrientation(config.restrict_orientations, new_rod)
	GetAxes(new_rod)
	GetVertices(config.n_dim, config.n_vertices, new_rod)
	new_weights := make([]float64, k)
	var new_weight_sum float64 = 0
	new_orientations := make([]float64, k)
	for i := 0; i < k; i++ {
		new_orientations[i] = new_rod.orientation
		no_overlaps := CheckNeighborOverlaps(new_rod, grid, rods, config)
		if no_overlaps {
			new_weights[i] = 1
			new_weight_sum += 1
		} else {
			new_weights[i] = 0
		}

		if i != (k - 1) {
			GetRandOrientation(config.restrict_orientations, new_rod)
			GetAxes(new_rod)
			GetVertices(config.n_dim, config.n_vertices, new_rod)
		}
	}

	// select new configuration
	new_rod.orientation = SelectWeightedConfig(new_orientations, new_weights, new_weight_sum, k)

	// run k-1 rosenbluth trials to generate weights for original configuration
	var og_weight_sum float64 = 1
	for i := 0; i < (k - 1); i++ {
		GetRandOrientation(config.restrict_orientations, og_rod)
		GetAxes(og_rod)
		GetVertices(config.n_dim, config.n_vertices, og_rod)
		no_overlaps := CheckNeighborOverlaps(og_rod, grid, rods, config)
		if no_overlaps {
			og_weight_sum += 1
		}
	}

	// calculate acceptance probability for rotation
	acc := new_weight_sum / og_weight_sum
	if rand.Float64() < acc {
		// move accepted
		GetAxes(new_rod)
		GetVertices(config.n_dim, config.n_vertices, new_rod)
		rods[rod.id] = new_rod
		config.rotation_successes++
	}
	config.rotation_attempts++
}

func Insert(grid []*GridSpace, config *Config, rods *[]*Rod) {
	// find rod id to insert and either create new rod or refresh rod state
	fromScratch := false
	var rod *Rod
	if len(config.avail_rod_ids) > 0 {
		rand_id := rand.Intn(len(config.avail_rod_ids))
		rod_id := config.avail_rod_ids[rand_id]
		rod = (*rods)[rod_id]
		RodRefresh(config, rod)
	} else {
		fromScratch = true
		rod_id := config.next_unused_rod_id
		rod = &Rod{}
		rod.id = rod_id
		RodInit(config, rod)
	}

	// set rosenbluth k to 1 for insertion
	k := 1

	// run k rosenbluth trials to generate insertion weights (same location - new orientation)
	new_weights := make([]float64, k)
	var new_weight_sum float64 = 0
	new_orientations := make([]float64, k)
	for i := 0; i < k; i++ {
		new_orientations[i] = rod.orientation
		no_overlaps := CheckNeighborOverlaps(rod, grid, *rods, config)
		if no_overlaps {
			new_weights[i] = 1
			new_weight_sum += 1
		} else {
			new_weights[i] = 0
		}

		if i != (k - 1) {
			GetRandOrientation(config.restrict_orientations, rod)
			GetAxes(rod)
			GetVertices(config.n_dim, config.n_vertices, rod)
		}
	}

	// select new configuration
	rod.orientation = SelectWeightedConfig(new_orientations, new_weights, new_weight_sum, k)

	// calculate acceptance probability for insertion
	N := float64(config.n_rods)
	acc := (config.V * math.Exp(config.beta*config.mu) / (N + 1)) * (new_weight_sum / float64(k))
	if rand.Float64() < acc {
		// move accepted
		rod.exists = true
		GetAxes(rod)
		GetVertices(config.n_dim, config.n_vertices, rod)

		// add rod to neighbor lists
		AddToNeighborLists(rod, grid)

		// update avail ids
		if !fromScratch {
			cur_rod_idx := 0
			for i := 0; i < len(config.avail_rod_ids); i++ {
				if config.avail_rod_ids[i] == rod.id {
					cur_rod_idx = i
					break
				}
			}
			config.avail_rod_ids = append(config.avail_rod_ids[:cur_rod_idx], config.avail_rod_ids[cur_rod_idx+1:]...)
		} else if fromScratch {
			config.next_unused_rod_id++
			(*rods) = append((*rods), rod)
		}
		config.n_rods++
		config.insertion_successes++
	} else {
		// move rejected
		rod.exists = false
	}
	config.insertion_attempts++
}

func Delete(rod *Rod, grid []*GridSpace, config *Config, rods []*Rod) {
	// calculate acceptance probability for deletion
	N := float64(config.n_rods)
	acc := (N / (config.V * math.Exp(config.beta*config.mu)))
	if rand.Float64() < acc {
		// if delete then remove rod from nearest neighbors lists
		RemFromNeighborLists(rod, grid)
		// and remove from rod list
		rod.exists = false
		config.avail_rod_ids = append(config.avail_rod_ids, rod.id)
		config.n_rods--
		config.deletion_successes++
	}
	config.deletion_attempts++
}
