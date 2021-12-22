package main

type GridSpace struct {
	grid_neighbors [9]int
	rod_neighbors  []int
}

func GetGridNeighbors(x int, y int, n_bins []int, grid_space *GridSpace) {
	x_neighbors := [3]int{x - 1, x, x + 1}
	y_neighbors := [3]int{y - 1, y, y + 1}
	x_bins := n_bins[0]
	y_bins := n_bins[1]
	if x_neighbors[2] > (x_bins - 1) {
		x_neighbors[2] -= x_bins
	}
	if x_neighbors[0] < 0 {
		x_neighbors[0] += x_bins
	}
	if y_neighbors[2] > (y_bins - 1) {
		y_neighbors[2] -= y_bins
	}
	if y_neighbors[0] < 0 {
		y_neighbors[0] += y_bins
	}
	i := 0
	for n := 0; n < 3; n++ {
		for l := 0; l < 3; l++ {
			grid_space.grid_neighbors[i] = x_neighbors[n] + y_neighbors[l]*x_bins
			i++
		}
	}
}

func GridInit(x int, y int, n_bins []int, grid_space *GridSpace) {
	GetGridNeighbors(x, y, n_bins, grid_space)
	grid_space.rod_neighbors = make([]int, 0, 80)
}
