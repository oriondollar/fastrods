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

	// calculate initial surface energy
	og_surface_energy := CalcSurfaceEnergy(og_rod, config)

	// move rod to new location
	GetRandLoc(config, new_rod)
	GetGridID(config.box_dims, config.n_bins, config.grid_bins, new_rod)
	new_grid_id := new_rod.grid_id

	// run k rosenbluth trials to generate swap weights (new location - new orientations)
	GetRandOrientation(config, new_rod)
	GetAxes(new_rod)
	GetVertices(config.n_dim, config.n_vertices, new_rod)
	new_weights := make([]float64, config.k)
	var new_weight_sum float64 = 0
	new_orientations := make([]float64, config.k)
	new_surface_energies := make([]float64, config.k)
	for i := 0; i < config.k; i++ {
		new_orientations[i] = new_rod.orientation
		surface_energy := CalcSurfaceEnergy(new_rod, config)
		no_overlaps := CheckNeighborOverlaps(new_rod, grid, rods, config)
		if no_overlaps {
			p := math.Exp(-config.beta * surface_energy)
			new_weights[i] = p
			new_weight_sum += p
			new_surface_energies[i] = surface_energy
		} else {
			new_weights[i] = 0
			new_surface_energies[i] = surface_energy + 10e8
		}

		if i != (config.k - 1) {
			GetRandOrientation(config, new_rod)
			GetAxes(new_rod)
			GetVertices(config.n_dim, config.n_vertices, new_rod)
		}
	}

	// select new configuration
	var new_surface_energy float64
	new_rod.orientation, new_surface_energy = SelectWeightedConfig(new_orientations, new_surface_energies, new_weights, new_weight_sum, config.k)

	// run k-1 rosenbluth trials to generate weights for original configuration
	var og_weight_sum float64 = math.Exp(-config.beta * og_surface_energy)
	for i := 0; i < (config.k - 1); i++ {
		GetRandOrientation(config, og_rod)
		GetAxes(og_rod)
		GetVertices(config.n_dim, config.n_vertices, og_rod)
		no_overlaps := CheckNeighborOverlaps(og_rod, grid, rods, config)
		if no_overlaps {
			surface_energy := CalcSurfaceEnergy(og_rod, config)
			og_weight_sum += math.Exp(-config.beta * surface_energy)
		}
	}

	// calculate acceptance probability for swap
	var acc float64
	if config.k == 1 {
		acc = math.Exp(-config.beta * (new_surface_energy - og_surface_energy))
	} else {
		acc = (og_weight_sum / new_weight_sum) * math.Exp(-config.beta*(new_surface_energy-og_surface_energy))
	}
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
		config.potential_energy += (new_surface_energy - og_surface_energy)
		config.swap_successes++
	}
	config.swap_attempts++
}

func Translate(rod *Rod, grid []*GridSpace, config *Config, rods []*Rod) {
	// store original rod location and grid id
	og_loc := rod.loc
	og_grid_id := rod.grid_id

	// get new location within 1 distance away from current rod COM
	new_loc := make([]float64, config.n_dim)
	if !config.restrict_translations {
		var max_r float64 = 1
		max_theta := 2 * math.Pi
		r := rand.Float64() * max_r
		theta := rand.Float64() * max_theta
		v := [2]float64{r * math.Sin(theta), r * math.Cos(theta)}
		new_loc[0] = og_loc[0] + v[0]
		new_loc[1] = og_loc[1] + v[1]
		if new_loc[0] >= config.box_dims[0] {
			new_loc[0] -= config.box_dims[0]
		}
		if new_loc[1] >= config.box_dims[1] {
			new_loc[1] -= config.box_dims[1]
		}
		if new_loc[0] < 0 {
			new_loc[0] += config.box_dims[0]
		}
		if new_loc[1] < 0 {
			new_loc[1] += config.box_dims[0]
		}
	} else {
		rand_move_idx := rand.Intn(len(config.lattice_moves))
		rand_move := config.lattice_moves[rand_move_idx]
		new_lattice_x := rod.lattice_x + rand_move[0]
		new_lattice_y := rod.lattice_y + rand_move[1]
		if new_lattice_x >= config.lattice_x {
			new_lattice_x -= config.lattice_x
		}
		if new_lattice_y >= config.lattice_y {
			new_lattice_y -= config.lattice_y
		}
		if new_lattice_x < 0 {
			new_lattice_x += config.lattice_x
		}
		if new_lattice_y < 0 {
			new_lattice_y += config.lattice_y
		}
		rod.lattice_x = new_lattice_x
		rod.lattice_y = new_lattice_y
		new_loc[0] = config.lattice_grid[rod.lattice_x][rod.lattice_y][0]
		new_loc[1] = config.lattice_grid[rod.lattice_x][rod.lattice_y][1]
		if new_loc[0] > config.box_dims[0] {
			new_loc[0] -= config.box_dims[0]
		}
		if new_loc[1] > config.box_dims[1] {
			new_loc[1] -= config.box_dims[1]
		}
		if new_loc[0] < 0 {
			new_loc[0] += config.box_dims[0]
		}
		if new_loc[1] < 0 {
			new_loc[1] += config.box_dims[0]
		}
	}

	rod.loc = new_loc
	GetGridID(config.box_dims, config.n_bins, config.grid_bins, rod)
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
	og_rod_orientation := og_rod.orientation
	new_rod := RodDeepCopy(rod)

	// calculate initial surface energy
	og_surface_energy := CalcSurfaceEnergy(og_rod, config)
	// fmt.Printf("Original Surface Energy: %f\n", og_surface_energy)

	// run k rosenbluth trials to generate rotation weights (same location - new orientations)
	GetRandOrientation(config, new_rod)
	GetAxes(new_rod)
	GetVertices(config.n_dim, config.n_vertices, new_rod)
	new_weights := make([]float64, config.k)
	var new_weight_sum float64 = 0
	new_orientations := make([]float64, config.k)
	new_surface_energies := make([]float64, config.k)
	for i := 0; i < config.k; i++ {
		new_orientations[i] = new_rod.orientation
		surface_energy := CalcSurfaceEnergy(new_rod, config)
		no_overlaps := CheckNeighborOverlaps(new_rod, grid, rods, config)
		if no_overlaps {
			p := math.Exp(-config.beta * surface_energy)
			// fmt.Printf("New Surface Energy: %f\n", surface_energy)
			new_weights[i] = p
			new_weight_sum += p
			new_surface_energies[i] = surface_energy
		} else {
			// fmt.Printf("New Surface Energy: %f\n", surface_energy+10e8)
			new_weights[i] = 0
			new_surface_energies[i] = surface_energy + 10e8
		}

		if i != (config.k - 1) {
			GetRandOrientation(config, new_rod)
			GetAxes(new_rod)
			GetVertices(config.n_dim, config.n_vertices, new_rod)
		}
	}

	// select new configuration
	var new_surface_energy float64
	new_rod.orientation, new_surface_energy = SelectWeightedConfig(new_orientations, new_surface_energies, new_weights, new_weight_sum, config.k)
	new_rod_orientation := new_rod.orientation

	// run k-1 rosenbluth trials to generate weights for original configuration
	var og_weight_sum float64 = math.Exp(-config.beta * og_surface_energy)
	for i := 0; i < (config.k - 1); i++ {
		GetRandOrientation(config, og_rod)
		GetAxes(og_rod)
		GetVertices(config.n_dim, config.n_vertices, og_rod)
		no_overlaps := CheckNeighborOverlaps(og_rod, grid, rods, config)
		if no_overlaps {
			surface_energy := CalcSurfaceEnergy(og_rod, config)
			og_weight_sum += math.Exp(-config.beta * surface_energy)
		}
	}

	// calculate acceptance probability for rotation
	var acc float64
	if config.k == 1 {
		acc = math.Exp(-config.beta * (new_surface_energy - og_surface_energy))
	} else {
		acc = (og_weight_sum / new_weight_sum) * math.Exp(-config.beta*(new_surface_energy-og_surface_energy))
	}
	// fmt.Printf("acceptance probability: %f\n", acc)
	if rand.Float64() < acc {
		// move accepted
		GetAxes(new_rod)
		GetVertices(config.n_dim, config.n_vertices, new_rod)
		rods[rod.id] = new_rod
		config.potential_energy += (new_surface_energy - og_surface_energy)
		if og_rod_orientation != new_rod_orientation {
			config.rotation_successes++
		}
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
	new_surface_energies := make([]float64, k)
	for i := 0; i < k; i++ {
		new_orientations[i] = rod.orientation
		surface_energy := CalcSurfaceEnergy(rod, config)
		no_overlaps := CheckNeighborOverlaps(rod, grid, *rods, config)
		if no_overlaps {
			p := math.Exp(-config.beta * surface_energy)
			new_weights[i] = p
			new_weight_sum += p
			new_surface_energies[i] = surface_energy
		} else {
			new_weights[i] = 0
			new_surface_energies[i] = surface_energy + 10e8
		}

		if i != (k - 1) {
			GetRandOrientation(config, rod)
			GetAxes(rod)
			GetVertices(config.n_dim, config.n_vertices, rod)
		}
	}

	// select new configuration
	var new_surface_energy float64
	rod.orientation, new_surface_energy = SelectWeightedConfig(new_orientations, new_surface_energies, new_weights, new_weight_sum, k)

	// calculate acceptance probability for insertion
	N := float64(config.n_rods)
	acc := (config.V * math.Exp(config.beta*(config.mu-new_surface_energy)) / (N + 1))
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
		config.potential_energy += new_surface_energy
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
	surface_energy := CalcSurfaceEnergy(rod, config)
	acc := (N / (config.V * math.Exp(config.beta*(config.mu-surface_energy))))
	if rand.Float64() < acc {
		// if delete then remove rod from nearest neighbors lists
		RemFromNeighborLists(rod, grid)
		// and remove from rod list
		rod.exists = false
		config.avail_rod_ids = append(config.avail_rod_ids, rod.id)
		config.potential_energy -= surface_energy
		config.n_rods--
		config.deletion_successes++
	}
	config.deletion_attempts++
}
