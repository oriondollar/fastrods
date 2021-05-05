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

	// TROUBLESHOOTING VERTICES & ROD WIDTH AFTER ROTATIONS
	// rods[0].loc = []float64{5., 5.}
	// // rods[0].orientation = (-math.Pi / 6) * (180 / math.Pi)
	// rods[0].orientation = 120.
	// GetAxes(rods[0])
	// GetVertices(config.n_dim, config.n_vertices, rods[0])
	// short_axis_length := math.Sqrt(math.Pow(rods[0].short_axis[0], 2) + math.Pow(rods[0].short_axis[1], 2))
	// fmt.Println("pos", rods[0].loc[0], rods[0].loc[1])
	// fmt.Println("long", rods[0].long_axis)
	// fmt.Println("short", rods[0].short_axis)
	// fmt.Println("short length", short_axis_length)
	// fmt.Println("rod mat", rods[0].rot_mat)
	// fmt.Println("og vertices", rods[0].vertical_vertices)
	// fmt.Println("vertices", rods[0].rotated_vertices)
	// width1 := math.Sqrt(math.Pow(rods[0].rotated_vertices[6]-rods[0].rotated_vertices[0], 2) + math.Pow(rods[0].rotated_vertices[7]-rods[0].rotated_vertices[1], 2))
	// width2 := math.Sqrt(math.Pow(rods[0].rotated_vertices[4]-rods[0].rotated_vertices[2], 2) + math.Pow(rods[0].rotated_vertices[5]-rods[0].rotated_vertices[3], 2))
	// length1 := math.Sqrt(math.Pow(rods[0].rotated_vertices[2]-rods[0].rotated_vertices[0], 2) + math.Pow(rods[0].rotated_vertices[3]-rods[0].rotated_vertices[1], 2))
	// length2 := math.Sqrt(math.Pow(rods[0].rotated_vertices[4]-rods[0].rotated_vertices[6], 2) + math.Pow(rods[0].rotated_vertices[5]-rods[0].rotated_vertices[7], 2))
	// fmt.Println("width1", width1)
	// fmt.Println("width2", width2)
	// fmt.Println("length1", length1)
	// fmt.Println("length2", length2)

	// fmt.Println()
	// // correct width 1
	// offset := width1 - 1.0
	// offsetY := rods[0].short_axis[1] * offset
	// y := rods[0].rotated_vertices[7]
	// y_adjusted := y - offsetY
	// fmt.Printf("offsetY - %.120f\n", offsetY)
	// fmt.Printf("original y - %.120f\n", y)
	// fmt.Printf("adjusted y - %.120f\n", y_adjusted)
	// fmt.Println()

	// offsetX := rods[0].short_axis[0] * offset
	// x := rods[0].rotated_vertices[6]
	// x_adjusted := x - offsetX
	// fmt.Printf("offsetX - %.120f\n", offsetX)
	// fmt.Printf("original x - %.120f\n", x)
	// fmt.Printf("adjusted x - %.120f\n", x_adjusted)
	// fmt.Println()

	// fmt.Printf("y2 - %.120f\n", rods[0].rotated_vertices[7])
	// Y2 := rods[0].rotated_vertices[7] - offsetY
	// fmt.Printf("Y2 - %.120f\n", Y2)
	// fmt.Println(rods[0].rotated_vertices[0], rods[0].rotated_vertices[1])
	// fmt.Println()
	// // new_width1 := math.Sqrt(math.Pow(rods[0].rotated_vertices[6]-rods[0].rotated_vertices[0], 2) + math.Pow(rods[0].rotated_vertices[7]-rods[0].rotated_vertices[1], 2))
	// fmt.Println("width1", width1)
	// // fmt.Println("new width1", new_width1)
	// fmt.Println("width2", width2)
	// fmt.Println("length1", length1)
	// fmt.Println("length2", length2)

	// fmt.Println("0 Select Rate -", config.n_rot_0/config.n_attempt_0*100)
	// fmt.Println("60 Select Rate -", config.n_rot_60/config.n_attempt_60*100)
	// fmt.Println("120 Select Rate -", config.n_rot_120/config.n_attempt_120*100)
	// fmt.Println("0 attempts -", config.n_attempt_0)
	// fmt.Println("60 attempts -", config.n_attempt_60)
	// fmt.Println("120 attempts -", config.n_attempt_120)

	// TROUBLESHOOTING CHECK OVERLAP
	// r_short := 0.
	// r_long := 8.99999999999999
	// short_theta := math.Pi / 2
	// long_theta := math.Pi / 3
	// new_orientation := 120.

	// rods[0].orientation = new_orientation
	// rods[0].loc = []float64{50., 50.}
	// GetAxes(rods[0])
	// GetVertices(config.n_dim, config.n_vertices, rods[0])

	// v_short := [2]float64{r_short * math.Sin(short_theta), r_short * math.Cos(short_theta)}
	// v_long := [2]float64{r_long * math.Sin(long_theta), r_long * math.Cos(long_theta)}
	// new_loc := make([]float64, config.n_dim)
	// new_loc[0] = rods[0].loc[0] + v_short[0] + v_long[0]
	// new_loc[1] = rods[0].loc[0] + v_short[1] + v_long[1]
	// rods[1].orientation = new_orientation
	// rods[1].loc = new_loc
	// GetAxes(rods[1])
	// GetVertices(config.n_dim, config.n_vertices, rods[1])

	// fmt.Println(rods[0].loc)
	// fmt.Println(rods[1].loc)

	// overlap := CheckOverlapTS(rods[0], rods[1], &config)
	// fmt.Println(overlap)

	MonteCarlo(&rods, grid, &config)

	fmt.Println(config.swap_successes / config.swap_attempts * 100)
	fmt.Println(config.rotation_successes / config.rotation_attempts * 100)
	fmt.Println(config.translation_successes / config.translation_attempts * 100)
	fmt.Println(config.insertion_successes / config.insertion_attempts * 100)
	fmt.Println(config.deletion_successes / config.deletion_attempts * 100)

}
