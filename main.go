package main

import (
	"fmt"
	"math/rand"
	"time"
)

var (
	MAX_RODS int = 2000
)

type GridSpace struct {
	grid_neighbors [9]int
	rod_neighbors  []int
}

type Rod struct {
	id                int
	loc               []float64
	orientation       float64
	length            float64
	width             float64
	length_by_2       float64
	width_by_2        float64
	long_axis         []float64
	short_axis        []float64
	rot_mat           []float64
	vertical_vertices []float64
	rotated_vertices  []float64
	grid_id           int
	exists            bool
}

func main() {
	fmt.Println("starting program...")
	// seed random number generator
	rand.Seed(time.Now().UnixNano())

	// reading parameter file
	fmt.Println("reading params...")
	config, err := ReadConfig("config.dat")
	Check(err)

	// initialize grid
	fmt.Println("initializing grid...")
	grid := make([]*GridSpace, config.n_grids)
	for i := 0; i < config.n_grids; i++ {
		x := i % config.n_bins
		y := i / config.n_bins
		grid[i] = &GridSpace{}
		GridInit(x, y, config.n_bins, grid[i])
	}

	// randomly place rods on surface
	fmt.Println("placing rods...")
	rods := make([]*Rod, config.n_rods, MAX_RODS)
	for i := 0; i < config.n_rods; i++ {
		// initialize rod
		rods[i] = &Rod{}
		rods[i].id = i

		// loop until rod is placed with no overlaps
		no_overlaps := false
		for !no_overlaps {
			// get random loc and orientation, set grid_id and vertices
			RodInit(&config, rods[i])

			// check nearest neighbor list to see if there are any overlaps
			no_overlaps = CheckNeighborOverlaps(rods[i], grid, rods, &config)
		}
		rods[i].exists = true

		// add rod to neighbor list of each neighboring grid
		gridspace := grid[rods[i].grid_id]
		for j := 0; j < len(gridspace.grid_neighbors); j++ {
			neighbor_id := gridspace.grid_neighbors[j]
			grid[neighbor_id].rod_neighbors = append(grid[neighbor_id].rod_neighbors, rods[i].id)
		}
	}

	MonteCarlo(&rods, grid, &config)

	fmt.Println(config.rotation_successes / config.rotation_attempts * 100)
	fmt.Println(config.translation_successes / config.translation_attempts * 100)
	fmt.Println(config.insertion_successes / config.insertion_attempts * 100)
	fmt.Println(config.deletion_successes / config.deletion_attempts * 100)

}
