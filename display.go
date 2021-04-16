package main

import (
	"github.com/gdamore/tcell/v2"
)

func GetDisplayGridPix(w, h int, x, y, box_length float64, xGridBins, yGridBins []float64) (int, int) {
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

	x_bin := w - 1
	y_bin := h - 1
	for i := 0; i < x_bin; i++ {
		bin := xGridBins[i]
		if x < bin {
			x_bin = i
			break
		}
	}
	for i := 0; i < y_bin; i++ {
		bin := yGridBins[i]
		if y < bin {
			y_bin = i
			break
		}
	}
	return x_bin, y_bin
}

func getRodPix(w, h int, rods []*Rod, config *Config) ([]int, []int) {
	xGridSpacing := float64(w) / config.box_length
	yGridSpacing := float64(h) / config.box_length
	xGridBins := make([]float64, w)
	yGridBins := make([]float64, h)
	for i := 0; i < w; i++ {
		xGridBins[i] = xGridSpacing * float64(i+1)
	}
	for i := 0; i < h; i++ {
		yGridBins[i] = yGridSpacing * float64(i+1)
	}
	minSpacing, _ := MinMax([]float64{xGridSpacing, yGridSpacing})
	nRodBins := int(config.rod_length / minSpacing)
	minSpacing = minSpacing / config.rod_length

	var xPix []int
	var yPix []int
	// will eventually loop through all rods
	for i := 0; i < len(rods); i++ {
		rod := rods[i]
		if rod.exists {
			x1 := rod.rotated_vertices[0] - ((rod.rotated_vertices[0] - rod.rotated_vertices[6]) / 2)
			y1 := rod.rotated_vertices[1] - ((rod.rotated_vertices[1] - rod.rotated_vertices[7]) / 2)
			x2 := rod.rotated_vertices[2] - ((rod.rotated_vertices[2] - rod.rotated_vertices[4]) / 2)
			y2 := rod.rotated_vertices[3] - ((rod.rotated_vertices[3] - rod.rotated_vertices[5]) / 2)
			v := [2]float64{x2 - x1, y2 - y1}
			dx := v[0] * minSpacing
			dy := v[1] * minSpacing
			wpix, hpix := GetDisplayGridPix(w, h, x1, y1, config.box_length, xGridBins, yGridBins)
			xPix = append(xPix, wpix)
			yPix = append(yPix, hpix)
			for j := 1; j < (nRodBins - 1); j++ {
				new_x := x1 + dx*float64(j)
				new_y := y1 + dy*float64(j)
				wpix, hpix := GetDisplayGridPix(w, h, new_x, new_y, config.box_length, xGridBins, yGridBins)
				xPix = append(xPix, wpix)
				yPix = append(yPix, hpix)
			}
			wpix, hpix = GetDisplayGridPix(w, h, x2, y2, config.box_length, xGridBins, yGridBins)
			xPix = append(xPix, wpix)
			yPix = append(yPix, hpix)
		}

	}
	return xPix, yPix
}

func MakeBox(s tcell.Screen, rodState []*Rod, config *Config) {
	w, h := s.Size()

	if w == 0 || h == 0 {
		return
	}

	// var xCOMS []int
	// var yCOMS []int

	box_end_w := w
	box_end_h := int(float64(h) * 0.95)
	pixw, pixh := getRodPix(box_end_w, box_end_h, rodState, config)

	st := tcell.StyleDefault
	st_rod := tcell.StyleDefault
	mainc := ' '
	if s.Colors() > 256 {
		rgb := tcell.NewRGBColor(111, 26, 8)
		st = st.Background(rgb)
		rgb_rod := tcell.NewRGBColor(156, 100, 20)
		st_rod = st_rod.Background(rgb_rod)
	} else if s.Colors() > 1 {
		rgb := tcell.NewRGBColor(111, 26, 8)
		st = st.Background(rgb)
		rgb_rod := tcell.NewRGBColor(156, 100, 20)
		st_rod = st_rod.Background(rgb_rod)
	}

	for row := 0; row < box_end_h; row++ {
		for col := 0; col < box_end_w; col++ {
			s.SetContent(col, row, mainc, nil, st)
		}
	}

	for i := 0; i < len(pixw); i++ {
		col := pixw[i]
		row := pixh[i]
		s.SetContent(col, row, mainc, nil, st_rod)
	}
	s.Show()
}
