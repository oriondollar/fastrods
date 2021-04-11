package display

import (
	"github.com/gdamore/tcell/v2"
)

func makebox(s tcell.Screen) {
	w, h := s.Size()

	if w == 0 || h == 0 {
		return
	}

	st := tcell.StyleDefault
	gl := ' '
	if s.Colors() > 256 {
		rgb := tcell.NewRGBColor(111, 26, 8)
		st := st.Background(rgb)
	}

	for row := 0; row < h; row++ {
		for col := 0; col < w; col++ {
			s.SetContent(col, row, gl, ' ', st)
		}
	}
	s.Show()
}
