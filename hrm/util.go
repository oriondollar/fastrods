package hrm

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
	box_length            float64
	V                     float64
	nn_cutoff             float64
	n_bins                int
	grid_spacing          float64
	n_grids               int
	nn_window             int
	grid_bins             []float64
	mc_alg                string
	cutoff_ratio          float64
	restrict_orientations bool
	temp                  float64
	kb                    float64
	beta                  float64
	mu                    float64
	n_cycles              int
	n_insert_deletes      int
	k                     int
	avail_rod_ids         []int
	next_unused_rod_id    int

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
	rot = -rot
	rad := rot * math.Pi / 180
	cos_rad := math.Cos(rad)
	sin_rad := math.Sin(rad)
	v_out = append(v_out, v[0]*cos_rad-v[1]*sin_rad)
	v_out = append(v_out, v[0]*sin_rad-v[1]*cos_rad)
	return
}

func SelectWeightedConfig(config_list []float64, weights []float64, weights_sum float64, k int) float64 {
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
	return selected_config
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
				} else if key == "box_length" {
					config.box_length, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "mc_alg" {
					config.mc_alg = value
				} else if key == "cutoff_ratio" {
					config.cutoff_ratio, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "restrict_orientations" {
					config.restrict_orientations, ke = strconv.ParseBool(value)
					Check(ke)
				} else if key == "temp" {
					config.temp, ke = strconv.ParseFloat(value, 64)
					Check(ke)
				} else if key == "mu" {
					config.mu, ke = strconv.ParseFloat(value, 64)
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
	config.V = math.Pow(config.box_length, float64(config.n_dim))
	config.rod_width = config.rod_length / config.aspect_ratio
	config.n_vertices = int(math.Pow(2, float64(config.n_dim)))
	config.nn_cutoff = config.rod_length * config.cutoff_ratio
	n_bins := math.Round(config.box_length / (config.nn_cutoff / 2))
	config.grid_spacing = config.box_length / n_bins
	config.n_bins = int(n_bins)
	for i := 0; i < config.n_bins; i++ {
		config.grid_bins = append(config.grid_bins, config.grid_spacing*float64(i+1))
	}
	config.n_grids = int(math.Pow(float64(config.n_bins), 2))
	config.kb = 1
	config.beta = 1 / (config.temp * config.kb)
	config.next_unused_rod_id = config.n_rods

	config.insertion_successes = 0
	config.insertion_attempts = 0
	config.deletion_successes = 0
	config.deletion_attempts = 0
	config.rotation_successes = 0
	config.rotation_attempts = 0
	config.translation_successes = 0
	config.translation_attempts = 0

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
