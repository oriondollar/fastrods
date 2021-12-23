package main

import (
	"fmt"
	"math/rand"
	"time"
	// _ "net/http/pprof"
)

var (
	MAX_RODS int     = 2000
	PI       float64 = 3.1415926535897932384626433832795028841971
)

func main() {
	// defer profile.Start(profile.TraceProfile, profile.ProfilePath(".")).Stop()
	// go func() {
	// 	log.Println(http.ListenAndServe("localhost:8080", nil))
	// }()

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
	for x := 0; x < config.n_bins[0]; x++ {
		for y := 0; y < config.n_bins[1]; y++ {
			i := x + y*config.n_bins[0]
			grid[i] = &GridSpace{}
			GridInit(x, y, config.n_bins, grid[i])
		}
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
		config.potential_energy += CalcSurfaceEnergy(rods[i], &config)

		// add rod to neighbor list of each neighboring grid
		gridspace := grid[rods[i].grid_id]
		for j := 0; j < len(gridspace.grid_neighbors); j++ {
			neighbor_id := gridspace.grid_neighbors[j]
			grid[neighbor_id].rod_neighbors = append(grid[neighbor_id].rod_neighbors, rods[i].id)
		}
	}

	// Prod
	// MonteCarlo(&rods, grid, &config)

	// WriteTraj(rods, "rod_init.dat")

	// Dev
	// var writer *bufio.Writer
	// out_file, err := os.Create("test_swap_k1.dat")
	// Check(err)
	// defer out_file.Close()
	// writer = bufio.NewWriter(out_file)
	// _, err = writer.WriteString("x,y,orientation\n")
	// Check(err)
	// for i := 0; i < 1000000; i++ {
	// 	rod := GetRandRod(rods)
	// 	// fmt.Printf("beta: %f\n", config.beta)
	// 	Swap(rod, grid, &config, rods)
	// 	_, err = writer.WriteString(fmt.Sprintf("%.3f,%.3f,%.3f\n", rod.loc[0], rod.loc[1], rod.orientation))
	// }
	// writer.Flush()
}
