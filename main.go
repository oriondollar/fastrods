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
}

func main() {
	fmt.Println("starting program...")
	// seed random number generator
	rand.Seed(time.Now().UnixNano())

	// reading parameter file
	fmt.Println("reading params...")
	config, err := ReadConfig("params.dat")
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
		RodInit(&config, rods[i], rods, grid)
	}

	WriteRodData(rods, "rod_locs.dat")

}
