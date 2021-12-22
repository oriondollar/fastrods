package main

import (
	"bufio"
	"fmt"
	"io"
	"math"
	"math/rand"
	"os"
	"strconv"
	"strings"
)

type Config struct {
	n_dim                 int
	n_rods                int
	rod_length            float64
	rod_width             float64
	aspect_ratio          float64
	n_vertices            int

	box_size           	  float64
	box_dims			  []float64
	V                     float64
	nn_cutoff             float64
	n_bins                []int
	grid_spacings         []float64
	n_grids               int
	grid_bins             [][]float64

	mc_alg                string
	cutoff_ratio          float64
	restrict_orientations bool
	restrict_translations bool

	lattice_pattern		  string
	facet_length		  float64
	lattice_x			  int
	lattice_y			  int
	lattice_grid		  [][][]float64
	lattice_moves		  [][]int

	temp                  float64
	kb                    float64
	beta                  float64
	mu                    float64

	n_cycles              int
	n_insert_deletes      int
	k                     int
	avail_rod_ids         []int
	next_unused_rod_id    int

	M               float64
	bias 			float64
	r_prime         float64
	overlap_penalty float64

	potential_energy	  float64
	swap_successes        float64
	swap_attempts         float64
	insertion_successes   float64
	insertion_attempts    float64
	deletion_successes    float64
	deletion_attempts     float64
	rotation_successes    float64
	rotation_attempts     float64
	translation_successes float64
	translation_attempts  float64

	write_CVs       bool
	write_traj      bool
	write_CV_freq   int
	write_traj_freq int
	CV_out          string
	traj_out        string

	print_proj bool
}

func Check(e error) {
	if e != nil {
		panic(e)
	}
}

func MinMax(array []float64) (float64, float64) {
	var max float64 = array[0]
	var min float64 = array[0]
	for _, value := range array {
		if max < value {
			max = value
		}
		if min > value {
			min = value
		}
	}
	return min, max
}

func RotateVector(v [2]float64, rot float64) (v_out []float64) {
	rad := rot * math.Pi / 180
	cos_rad := math.Cos(rad)
	sin_rad := math.Sin(rad)
	v_out = append(v_out, v[0]*cos_rad-v[1]*sin_rad)
	v_out = append(v_out, v[0]*sin_rad-v[1]*cos_rad)
	return
}

func SelectWeightedConfig(config_list []float64, value_list []float64, weights []float64, weights_sum float64, k int) (float64, float64) {
	rand_weight := rand.Float64() * weights_sum
	var running_weight float64 = 0
	idx_sel := k - 1
	for i := 0; i < k; i++ {
		running_weight += weights[i]
		if running_weight > rand_weight {
			idx_sel = i
			break
		}
	}
	selected_config := config_list[idx_sel]
	selected_value := value_list[idx_sel]
	return selected_config, selected_value
}

func RodDeepCopy(rod *Rod) *Rod {
	rodCopy := *rod
	rodCopy.loc = make([]float64, len(rod.loc))
	rodCopy.long_axis = make([]float64, len(rod.long_axis))
	rodCopy.short_axis = make([]float64, len(rod.short_axis))
	rodCopy.rotated_vertices = make([]float64, len(rod.rotated_vertices))
	copy(rodCopy.loc, rod.loc)
	copy(rodCopy.long_axis, rod.long_axis)
	copy(rodCopy.short_axis, rod.short_axis)
	copy(rodCopy.rotated_vertices, rod.rotated_vertices)
	return &rodCopy
}

func ReadConfig(fn string) (config Config, err error) {
	f, e := os.Open(fn)
	Check(e)
	defer f.Close()
	config.CV_out = "COLVAR.dat"
	config.traj_out = "traj.dat"

	reader := bufio.NewReader(f)

	for {
		line, e := reader.ReadString('\n')
		var ke error
		if equal := strings.Index(line, ":"); equal >= 0 {
			if key := strings.TrimSpace(line[:equal]); len(key) > 0 {
				value := ""
				if len(line) > equal {
					value = strings.TrimSpace(line[equal+1:])
				}
				if key == "n_dim" {
					config.n_dim, ke = strconv.Atoi(value)
					Check(ke)
				} else if key == "n_rods" {
					config.n_rods, ke = strconv.Atoi(value)
					Check(ke)
				} else if key == "rod_length" {
					config.rod_length, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "aspect_ratio" {
					config.aspect_ratio, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "box_size" {
					config.box_size, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "mc_alg" {
					config.mc_alg = value
				} else if key == "cutoff_ratio" {
					config.cutoff_ratio, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "restrict_orientations" {
					config.restrict_orientations, ke = strconv.ParseBool(value)
					Check(ke)
				} else if key == "restrict_translations" {
					config.restrict_translations, ke = strconv.ParseBool(value)
					Check(ke)
				} else if key == "lattice_pattern" {
					config.lattice_pattern = value
				} else if key == "facet_length" {
					config.facet_length, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "temp" {
					config.temp, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "mu" {
					config.mu, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "M" {
					config.M, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "bias"{
					config.bias, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "r_prime" {
					config.r_prime, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "n_cycles" {
					config.n_cycles, ke = strconv.Atoi(value)
					Check(ke)
				} else if key == "n_insert_deletes" {
					config.n_insert_deletes, ke = strconv.Atoi(value)
					Check(ke)
				} else if key == "k" {
					config.k, ke = strconv.Atoi(value)
					Check(ke)
				} else if key == "write_CVs" {
					config.write_CVs, ke = strconv.ParseBool(value)
					Check(ke)
				} else if key == "write_traj" {
					config.write_traj, ke = strconv.ParseBool(value)
					Check(ke)
				} else if key == "write_CV_freq" {
					config.write_CV_freq, ke = strconv.Atoi(value)
					Check(ke)
				} else if key == "write_traj_freq" {
					config.write_traj_freq, ke = strconv.Atoi(value)
					Check(ke)
				} else if key == "CV_out" {
					config.CV_out = value
				} else if key == "traj_out" {
					config.traj_out = value
				}
			}
		}
		if e == io.EOF {
			break
		}
		if e != nil {
			return config, e
		}
	}
	// Rod Geometry
	config.rod_width = config.rod_length / config.aspect_ratio
	config.n_vertices = int(math.Pow(2, float64(config.n_dim)))
	config.nn_cutoff = config.rod_length * config.cutoff_ratio

	// Lattice Geometry
	config.box_dims = make([]float64, config.n_dim)
	if config.restrict_translations == false {
		for i := 0; i < config.n_dim; i++ {
			config.box_dims[i] = config.box_size
		}
	} else if config.restrict_translations == true {
		if config.lattice_pattern == "triangular" {
			facet_height := math.Sqrt(3) / 2 * config.facet_length
			dx := config.facet_length
			dy := 2 * facet_height
			lattice_x := math.Floor(config.box_size / dx)
			lattice_y := 2 * math.Floor(config.box_size / dy)
			config.box_dims[0] = lattice_x * dx
			config.box_dims[1] = lattice_y * facet_height
			config.lattice_x = int(lattice_x)
			config.lattice_y = int(lattice_y)
			config.lattice_grid = make([][][]float64, config.lattice_x)
			for i := 0; i < config.lattice_x; i++ {
				config.lattice_grid[i] = make([][]float64, config.lattice_y)
				for j := 0; j < config.lattice_y; j++ {
					config.lattice_grid[i][j] = make([]float64, config.n_dim)
					var x_coord float64
					var y_coord float64
					if j % 2 == 0 {
						x_coord = float64(i) * dx + (dx / 4)
					} else {
						x_coord = (float64(i) * dx) + (dx / 2 + dx / 4)
					}
					y_coord = float64(j) * (dy / 2) + (dy / 4)
					config.lattice_grid[i][j][0] = x_coord
					config.lattice_grid[i][j][1] = y_coord
				}
			}
			lattice_moves := [11][2]int{{0, 0},
										{0, 2}, {0, -2},
										{1, 0}, {-1, 0},
										{0, 1}, {0, -1},
										{-1, 1}, {-1, -1},
										{1, 1}, {1, -1}}
			config.lattice_moves = make([][]int, 11)
			for i := 0; i < 11; i++ {
				config.lattice_moves[i] = make([]int, config.n_dim)
				config.lattice_moves[i][0] = lattice_moves[i][0]
				config.lattice_moves[i][1] = lattice_moves[i][1]
			}
		} else {
			fmt.Println("Grid pattern must be triangular")
		}
	}
	config.V = 1
	config.n_bins = make([]int, config.n_dim)
	config.grid_spacings = make([]float64, config.n_dim)
	config.grid_bins = make([][]float64, config.n_dim)
	config.n_grids = 1
	for i := 0; i < config.n_dim; i++ {
		config.V *= config.box_dims[i]
		n_bins := math.Floor(config.box_dims[i] / config.rod_length)
		config.grid_spacings[i] = config.box_dims[i] / n_bins
		config.n_bins[i] = int(n_bins)
		config.grid_bins[i] = make([]float64, config.n_bins[i])
		for j := 0; j < config.n_bins[i]; j++ {
			config.grid_bins[i][j] = config.grid_spacings[i]*float64(j+1)
		}
		config.n_grids *= config.n_bins[i]
	}

	// n_bins := math.Floor(config.box_size / config.rod_length)
	// config.grid_spacing = config.box_size / n_bins
	// config.n_bins = int(n_bins)
	// for i := 0; i < config.n_bins; i++ {
	// 	config.grid_bins = append(config.grid_bins, config.grid_spacing*float64(i+1))
	// }
	// config.n_grids = int(math.Pow(float64(config.n_bins), 2))
	config.kb = 1
	config.beta = 1 / (config.temp * config.kb)
	config.overlap_penalty = 1000000000000
	config.next_unused_rod_id = config.n_rods

	config.potential_energy = 0
	config.swap_successes = 0
	config.swap_attempts = 0
	config.insertion_successes = 0
	config.insertion_attempts = 0
	config.deletion_successes = 0
	config.deletion_attempts = 0
	config.rotation_successes = 0
	config.rotation_attempts = 0
	config.translation_successes = 0
	config.translation_attempts = 0

	config.print_proj = false

	return
}

func WriteTraj(rods []*Rod, fo string) {
	f, err := os.Create(fo)
	Check(err)
	defer f.Close()

	w := bufio.NewWriter(f)
	_, err = w.WriteString("id,x,y,orientation\n")
	Check(err)

	for i := 0; i < len(rods); i++ {
		rod := rods[i]
		id := rod.id
		x := rod.loc[0]
		y := rod.loc[1]
		orientation := rod.orientation
		if rod.exists {
			_, err = w.WriteString(fmt.Sprintf("%v,%.3f,%.3f,%.3f\n", id, x, y, orientation))
			Check(err)
		}

	}

	w.Flush()
}
