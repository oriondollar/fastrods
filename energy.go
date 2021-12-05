package main

import (
	"math"
)

func CalcSurfaceEnergy(rod *Rod, config *Config) float64 {
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
	return max_energy - 1.
}
