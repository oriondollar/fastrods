package main

import (
	"fmt"
	"math"
	"math/rand"

	"gonum.org/v1/gonum/mat"
)

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
	static            bool
}

func GetRandRod(rods []*Rod) (rod *Rod) {
	exists := false
	for !exists {
		rod_id := rand.Intn(len(rods))
		rod = rods[rod_id]
		exists = rod.exists
	}
	return rod
}

func GetRandLoc(n_dim int, box_length float64, rod *Rod) {
	x := rand.Float64() * box_length
	y := rand.Float64() * box_length
	if n_dim == 2 {
		rod.loc = []float64{x, y}
	} else if n_dim == 3 {
		z := rand.Float64() * box_length
		rod.loc = []float64{x, y, z}
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
	v := [2]float64{1, 0}
	rod.long_axis = RotateVector(v, rod.orientation)
	rod.short_axis = []float64{-rod.long_axis[1], rod.long_axis[0]}
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

func RodRefresh(config *Config, rod *Rod) {
	GetRandLoc(config.n_dim, config.box_length, rod)
	GetRandOrientation(config.restrict_orientations, rod)
	GetGridID(config.box_length, config.n_bins, &config.grid_bins, rod)
	GetAxes(rod)
	GetVertices(config.n_dim, config.n_vertices, rod)
}

func RodInit(config *Config, rod *Rod) {
	// set rod struct variables
	rod.length = config.rod_length
	rod.width = config.rod_width
	rod.length_by_2 = rod.length / 2
	rod.width_by_2 = rod.width / 2
	rod.long_axis = make([]float64, config.n_dim)
	rod.short_axis = make([]float64, config.n_dim)
	rod.rot_mat = make([]float64, config.n_dim*config.n_dim)
	rod.vertical_vertices = make([]float64, config.n_vertices*config.n_dim)
	all_vertices := []float64{-rod.length_by_2, rod.width_by_2, rod.length_by_2, rod.width_by_2, rod.length_by_2, -rod.width_by_2, -rod.length_by_2, -rod.width_by_2}
	rod.vertical_vertices = all_vertices
	rod.rotated_vertices = make([]float64, config.n_vertices*config.n_dim)
	RodRefresh(config, rod)
}

func CheckOverlap(rod1 *Rod, rod2 *Rod, config *Config) bool {
	if rod1.id == rod2.id {
		return false
	}
	move_rod := false
	var r [2]float64
	r[0] = rod1.loc[0] - rod2.loc[0]
	r[1] = rod1.loc[1] - rod2.loc[1]
	for d := 0; d < len(rod1.loc); d++ {
		f := math.Round(r[d] / config.box_length)
		if math.Round(f) != 0 {
			r[d] = r[d] - config.box_length*f
			move_rod = true
		}
	}
	if move_rod {
		rod2 = RodDeepCopy(rod2)
		rod2.loc[0] = rod1.loc[0] - r[0]
		rod2.loc[1] = rod1.loc[1] - r[1]
		GetAxes(rod2)
		GetVertices(config.n_dim, config.n_vertices, rod2)
	}
	projections := make([][]float64, config.n_dim)
	for i := 0; i < config.n_dim; i++ {
		projections[i] = make([]float64, config.n_dim*2)
	}
	axes := make([][]float64, config.n_dim*2)
	for i := 0; i < config.n_dim*2; i++ {
		axes[i] = make([]float64, config.n_dim)
	}
	axes[:][0] = rod1.long_axis
	axes[:][1] = rod1.short_axis
	axes[:][2] = rod2.long_axis
	axes[:][3] = rod2.short_axis
	overlap := true
	for i := 0; i < config.n_dim*2; i++ {
		for j := 0; j < config.n_vertices; j++ {
			projections[0][j] = axes[i][0]*rod1.rotated_vertices[j*2] + axes[i][1]*rod1.rotated_vertices[(j*2)+1]
			projections[1][j] = axes[i][0]*rod2.rotated_vertices[j*2] + axes[i][1]*rod2.rotated_vertices[(j*2)+1]
		}

		if config.print_proj {
			fmt.Println(projections)
		}
		min_proj_1, max_proj_1 := MinMax(projections[0][:])
		min_proj_2, max_proj_2 := MinMax(projections[1][:])

		if (min_proj_1 > max_proj_2) || (max_proj_1 < min_proj_2) {
			overlap = false
			return overlap
		}
	}
	return overlap
}

// func CheckOverlapFrenkel(rod1 *Rod, rod2 *Rod, config *Config) bool {
// 	if rod1.id == rod2.id {
// 		return false
// 	}
// 	move_rod := false
// 	var r [2]float64
// 	r[0] = rod1.loc[0] - rod2.loc[0]
// 	r[1] = rod1.loc[1] - rod2.loc[1]
// 	for d := 0; d < len(rod1.loc); d++ {
// 		f := math.Round(r[d] / config.box_length)
// 		if math.Round(f) != 0 {
// 			r[d] = r[d] - config.box_length*f
// 			move_rod = true
// 		}
// 	}
// 	if move_rod {
// 		rod2 = RodDeepCopy(rod2)
// 		rod2.loc[0] = rod1.loc[0] - r[0]
// 		rod2.loc[1] = rod1.loc[1] - r[1]
// 		GetAxes(rod2)
// 		GetVertices(config.n_dim, config.n_vertices, rod2)
// 	}

// 	v_i := rod1.short_axis
// 	v_j := rod2.short_axis
// 	rhs := (math.Pow(config.rod_length, 2) / 4) * (1 - math.Pow(v_i[0]*v_j[0]+v_i[1]*v_j[1], 2))
// 	g_i := math.Pow(r[0]*v_i[0]+r[1]*v_i[1], 2) - rhs
// 	g_j := math.Pow(r[0]*v_j[0]+r[1]*v_j[1], 2) - rhs

// 	if (g_i < 0) && (g_j < 0) {
// 		return true
// 	} else {
// 		return false
// 	}
// }

func CheckNeighborOverlaps(rod *Rod, grid []*GridSpace, rods []*Rod, config *Config) bool {
	rod_neighbors := grid[rod.grid_id].rod_neighbors
	// rod_neighbor_count := len(rod_neighbors)
	// result := make(chan bool)
	// for i := 0; i < rod_neighbor_count; i++ {
	// 	go func(rod1 *Rod, rod2 *Rod, config *Config) {
	// 		overlap := CheckOverlap(rod1, rod2, config)
	// 		result <- overlap
	// 	}(rod, rods[grid[rod.grid_id].rod_neighbors[i]], config)
	// }
	// overlaps := false
	// for i := 0; i < rod_neighbor_count; i++ {
	// 	is_overlap := <-result
	// 	if is_overlap {
	// 		overlaps = true
	// 	}
	// }
	// return !overlaps
	no_overlaps := true
	for i := 0; i < len(rod_neighbors); i++ {
		overlap := CheckOverlap(rod, rods[grid[rod.grid_id].rod_neighbors[i]], config)
		if overlap {
			no_overlaps = false
			return no_overlaps
		}
	}
	return no_overlaps
}

func RemFromNeighborLists(rod *Rod, grid []*GridSpace) {
	for i := 0; i < len(grid[rod.grid_id].grid_neighbors); i++ {
		grid_id := grid[rod.grid_id].grid_neighbors[i]
		rod_neighbors := grid[grid_id].rod_neighbors
		cur_rod_idx := 0
		for j := 0; j < len(rod_neighbors); j++ {
			if rod_neighbors[j] == rod.id {
				cur_rod_idx = j
				break
			}
		}
		new_rod_neighbors := make([]int, 0)
		new_rod_neighbors = append(new_rod_neighbors, grid[grid_id].rod_neighbors[:cur_rod_idx]...)
		new_rod_neighbors = append(new_rod_neighbors, grid[grid_id].rod_neighbors[cur_rod_idx+1:]...)
		grid[grid_id].rod_neighbors = new_rod_neighbors
	}
}

func AddToNeighborLists(rod *Rod, grid []*GridSpace) {
	for i := 0; i < len(grid[rod.grid_id].grid_neighbors); i++ {
		grid_id := grid[rod.grid_id].grid_neighbors[i]
		grid[grid_id].rod_neighbors = append(grid[grid_id].rod_neighbors, rod.id)
	}
}
