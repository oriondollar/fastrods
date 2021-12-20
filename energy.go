package main

import (
	"math"
)

func CalcSurfaceEnergy(rod *Rod, config *Config) float64 {
	rad := rod.orientation * PI / 180
	v1 := (math.Cos(6*(rad+(PI/3))) + 1) / 8
	var v2 float64
	if (rod.orientation < 15.) || (rod.orientation > 45.) {
		v2 = 0.
	} else {
		v2 = config.bias * (math.Cos(12*(rad+(PI/12))) - 1) / 16
	}
	energy := v1 + v2
	return energy
}

func CalcSurfaceEnergy_MinMax(rod *Rod, config *Config) float64 {
	rot1 := 0.
	rot2 := 90 * (1 + config.r_prime)
	rot3 := 90 * (2 + (1 - config.r_prime))
	v1 := RotateVector([2]float64{1, 0}, rot1)
	v2 := RotateVector([2]float64{1, 0}, rot2)
	v3 := RotateVector([2]float64{1, 0}, rot3)

	up_energy := math.Pow(v1[0]*rod.long_axis[0]+v1[1]*rod.long_axis[1], 2)
	dr_energy := math.Pow(v2[0]*rod.long_axis[0]+v2[1]*rod.long_axis[1], 2)
	dl_energy := math.Pow(v3[0]*rod.long_axis[0]+v3[1]*rod.long_axis[1], 2)
	energy_arr := []float64{up_energy, dr_energy, dl_energy}
	_, max_energy := MinMax(energy_arr)
	return max_energy - 0.75
}
