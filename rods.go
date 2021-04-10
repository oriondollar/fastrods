package main

import (
	"math"
	"math/rand"

	"gonum.org/v1/gonum/mat"
)

func GetRandLoc(n_dim int, box_length float64, rod *Rod) {
	x := rand.Float64() * box_length
	y := rand.Float64() * box_length
	if n_dim == 2 {
		rod.loc = append((*rod).loc, x, y)
	} else if n_dim == 3 {
		z := rand.Float64() * box_length
		rod.loc = append((*rod).loc, x, y, z)
	}
}

func GetRandOrientation(restricted bool, rod *Rod) {
	if restricted {
		rand_w := rand.Float64()
		if rand_w < (1. / 3.) {
			rod.orientation = 0.
		} else if rand_w < (2. / 3.) {
			rod.orientation = 60.
		} else {
			rod.orientation = 120.
		}
	} else {
		rod.orientation = rand.Float64() * 180
	}
}

func GetGridID(box_length float64, n_bins int, grid_bins *[]float64, rod *Rod) {
	x := rod.loc[0]
	y := rod.loc[1]
	if x > box_length {
		x -= box_length
	} else if x < 0 {
		x += box_length
	}
	if y > box_length {
		y -= box_length
	} else if y < 0 {
		y += box_length
	}

	found_x := false
	found_y := false
	x_bin := n_bins - 1
	y_bin := n_bins - 1
	for i := 0; i < (n_bins - 1); i++ {
		bin := (*grid_bins)[i]
		if (x < bin) && (!found_x) {
			x_bin = i
			found_x = true
		}
		if (y < bin) && (!found_y) {
			y_bin = i
			found_y = true
		}
		if found_x && found_y {
			break
		}
	}
	rod.grid_id = x_bin + y_bin*n_bins
}

func GetAxes(rod *Rod) {
	v := [2]float64{0, 1}
	rod.long_axis = RotateVector(v, rod.orientation)
	rod.short_axis = RotateVector(v, rod.orientation-90)
	rod.rot_mat = rod.long_axis
	rod.rot_mat = append(rod.rot_mat, rod.short_axis...)
}

func GetVertices(n_dim int, n_vertices int, rod *Rod) {
	rot_mat := mat.NewDense(n_dim, n_dim, rod.rot_mat)
	vert_vx_mat := mat.NewDense(n_vertices, n_dim, rod.vertical_vertices)
	rot_vx_mat := mat.NewDense(n_vertices, n_dim, rod.rotated_vertices)
	rot_vx_mat.Mul(vert_vx_mat, rot_mat)
	for i := 0; i < len(rod.rotated_vertices); i++ {
		rod.rotated_vertices[i] += rod.loc[i%2]
	}
}

func RodInit(config *Config, rod *Rod) {
	rod.length = config.rod_length
	rod.width = config.rod_width
	rod.length_by_2 = rod.length / 2
	rod.width_by_2 = rod.width / 2
	rod.long_axis = make([]float64, config.n_dim)
	rod.short_axis = make([]float64, config.n_dim)
	rod.rot_mat = make([]float64, config.n_dim*config.n_dim)
	n_vertices := int(math.Pow(2, float64(config.n_dim)))
	rod.vertical_vertices = make([]float64, n_vertices*config.n_dim)
	all_vertices := []float64{-rod.length_by_2, rod.width_by_2, rod.length_by_2, rod.width_by_2, rod.length_by_2, -rod.width_by_2, -rod.length_by_2, -rod.width_by_2}
	rod.vertical_vertices = all_vertices
	rod.rotated_vertices = make([]float64, n_vertices*config.n_dim)
	GetRandLoc(config.n_dim, config.box_length, rod)
	GetRandOrientation(config.restrict_orientations, rod)
	GetGridID(config.box_length, config.n_bins, &config.grid_bins, rod)
	GetAxes(rod)
	GetVertices(config.n_dim, n_vertices, rod)
}
