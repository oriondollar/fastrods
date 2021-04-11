package hrm

import (
	"gonum.org/v1/gonum/mat"
)

func CalcDensity(config *Config) float64 {
	density := float64(config.n_rods) / config.V
	return density
}

func CalcS(rods []*Rod, config *Config) float64 {
	Q := []float64{0, 0, 0, 0}
	for i := 0; i < len(rods); i++ {
		long_axis := rods[i].long_axis
		if rods[i].exists {
			Q[0] += (2*long_axis[0]*long_axis[0] - 1)
			Q[1] += (2 * long_axis[0] * long_axis[1])
			Q[2] += (2 * long_axis[1] * long_axis[0])
			Q[3] += (2*long_axis[1]*long_axis[1] - 1)
		}
	}
	Q[0] /= float64(config.n_rods)
	Q[1] /= float64(config.n_rods)
	Q[2] /= float64(config.n_rods)
	Q[3] /= float64(config.n_rods)
	var Q_mat mat.Matrix
	Q_mat = mat.NewDense(2, 2, Q)
	var eig mat.Eigen
	eig.Factorize(Q_mat, 1)
	var S float64 = 0
	for _, v := range eig.Values(nil) {
		if real(v) > S {
			S = real(v)
		}
	}
	return S
}
