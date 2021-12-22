package main

import (
	"bufio"
	"fmt"
	"math/rand"
	"os"
	"time"

	"github.com/cheggaaa/pb/v3"
)

func MonteCarlo(rods *[]*Rod, grid []*GridSpace, config *Config) {
	// set up writers
	var CV_writer *bufio.Writer
	var traj_writer *bufio.Writer
	var move_probs_writer *bufio.Writer
	var err error
	if config.write_CVs {
		CV_file, err := os.Create(config.CV_out)
		Check(err)
		defer CV_file.Close()

		CV_writer = bufio.NewWriter(CV_file)
		_, err = CV_writer.WriteString("cycle,S,density,energy\n")
		Check(err)
	}

	if config.write_traj {
		traj_file, err := os.Create(config.traj_out)
		Check(err)
		defer traj_file.Close()

		traj_writer = bufio.NewWriter(traj_file)
		_, err = traj_writer.WriteString("cycle,id,x,y,orientation\n")
		Check(err)
	}

	if config.write_move_probs {
		move_probs_file, err := os.Create(config.move_probs_out)
		Check(err)
		defer move_probs_file.Close()

		move_probs_writer = bufio.NewWriter(move_probs_file)
		_, err = move_probs_writer.WriteString("cycle,swap,rot,trans,insert,delete\n")
	}

	// MC loop
	bar := pb.StartNew(config.n_cycles)
	bar.SetRefreshRate(time.Second)
	for i := 0; i < config.n_cycles; i++ {
		for j := 0; j < config.n_rods; j++ {
			move_prob := rand.Float64()
			if move_prob < (1. / 3.) {
				rod := GetRandRod(*rods)
				Rotate(rod, grid, config, *rods)
			} else if move_prob < (2. / 3.) {
				rod := GetRandRod(*rods)
				Translate(rod, grid, config, *rods)
			} else {
				rod := GetRandRod(*rods)
				Swap(rod, grid, config, *rods)
			}
		}
		if config.mc_alg == "grand_canonical" {
			for j := 0; j < config.n_insert_deletes; j++ {
				remove_prob := rand.Float64()
				if remove_prob < (1. / 2.) {
					Insert(grid, config, rods)
				} else {
					if config.n_rods != 0 {
						rod := GetRandRod(*rods)
						Delete(rod, grid, config, *rods)
					}
				}
			}
		}
		// write results
		if (config.write_CVs) && ((i+1)%config.write_CV_freq == 0) {
			density := CalcDensity(config)
			S := CalcS(*rods, config)
			_, err = CV_writer.WriteString(fmt.Sprintf("%v,%.3f,%.3f,%.3f\n", i+1, S, density, config.potential_energy))
			Check(err)
		}
		if (config.write_traj) && ((i+1)%config.write_traj_freq == 0) {
			for j := 0; j < len(*rods); j++ {
				id := (*rods)[j].id
				x := (*rods)[j].loc[0]
				y := (*rods)[j].loc[1]
				orientation := (*rods)[j].orientation
				if (*rods)[j].exists {
					_, err = traj_writer.WriteString(fmt.Sprintf("%v,%v,%.16f,%.16f,%.3f\n", i+1, id, x, y, orientation))
					Check(err)
				}
			}
		}
		if (config.write_move_probs) && ((i+1)%config.write_move_probs_freq == 0) {
			swap_prob := config.swap_successes / config.swap_attempts
			rot_prob := config.rotation_successes / config.rotation_attempts
			trans_prob := config.translation_successes / config.translation_attempts
			insert_prob := config.insertion_successes / config.insertion_attempts
			delete_prob := config.deletion_successes / config.deletion_attempts
			_, err = move_probs_writer.WriteString(fmt.Sprintf("%v,%.5f,%.5f,%.5f,%.5f,%.5f\n", i+1, swap_prob, rot_prob, trans_prob, insert_prob, delete_prob))
			Check(err)
		}
		if config.write_CVs {
			CV_writer.Flush()
		}
		if config.write_traj {
			traj_writer.Flush()
		}
		if config.write_move_probs {
			move_probs_writer.Flush()
		}
		bar.Increment()
	}
	bar.Finish()
}
