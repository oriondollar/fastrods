package main

import (
	"fmt"
	"math/rand"
	"time"
)

var (
	MAX_RODS int     = 2000
	PI       float64 = 3.1415926535897932384626433832795028841971
)

func main() {
	// defer profile.Start(profile.TraceProfile, profile.ProfilePath(".")).Stop()

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

	// WriteTraj(rods, "rod_init.dat")

	// rods[0].loc = []float64{0.05, 0.05}
	// rods[0].orientation = (-PI / 6) * (180 / PI)
	// // rods[0].orientation = 0.
	// GetAxes(rods[0])
	// GetVertices(config.n_dim, config.n_vertices, rods[0])
	// fmt.Println("pos", rods[0].loc[0], rods[0].loc[1])
	// fmt.Println("long", rods[0].long_axis)
	// fmt.Println("short", rods[0].short_axis)
	// fmt.Println("rod mat", rods[0].rot_mat)
	// fmt.Println("og vertices", rods[0].vertical_vertices)
	// fmt.Println("vertices", rods[0].rotated_vertices)
	// width1 := math.Sqrt(math.Pow(rods[0].rotated_vertices[6]-rods[0].rotated_vertices[0], 2) + math.Pow(rods[0].rotated_vertices[7]-rods[0].rotated_vertices[1], 2))
	// width2 := math.Sqrt(math.Pow(rods[0].rotated_vertices[4]-rods[0].rotated_vertices[2], 2) + math.Pow(rods[0].rotated_vertices[5]-rods[0].rotated_vertices[3], 2))
	// fmt.Println("width1", float32(width1))
	// fmt.Println("width2", float32(width2))

	// fmt.Println("0 Select Rate -", config.n_rot_0/config.n_attempt_0*100)
	// fmt.Println("60 Select Rate -", config.n_rot_60/config.n_attempt_60*100)
	// fmt.Println("120 Select Rate -", config.n_rot_120/config.n_attempt_120*100)
	// fmt.Println("0 attempts -", config.n_attempt_0)
	// fmt.Println("60 attempts -", config.n_attempt_60)
	// fmt.Println("120 attempts -", config.n_attempt_120)

	MonteCarlo(&rods, grid, &config)

	fmt.Println(config.swap_successes / config.swap_attempts * 100)
	fmt.Println(config.rotation_successes / config.rotation_attempts * 100)
	fmt.Println(config.translation_successes / config.translation_attempts * 100)
	fmt.Println(config.insertion_successes / config.insertion_attempts * 100)
	fmt.Println(config.deletion_successes / config.deletion_attempts * 100)

}
