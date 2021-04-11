package display

import (
	"github.com/gdamore/tcell/v2"
	"github.com/oriondollar/fastrods/hrm"
)

func MakeBox(s tcell.Screen, rodState []*hrm.Rod) {
	w, h := s.Size()

	if w == 0 || h == 0 {
		return
	}

	box_end_w := w
	box_end_h := int(float64(h) * 0.95)

	st := tcell.StyleDefault
	mainc := ' '
	if s.Colors() > 256 {
		rgb := tcell.NewRGBColor(111, 26, 8)
		st = st.Background(rgb)
	} else if s.Colors() > 1 {
		rgb := tcell.NewRGBColor(111, 26, 8)
		st = st.Background(rgb)
	}

	for row := 0; row < box_end_h; row++ {
		for col := 0; col < box_end_w; col++ {
			s.SetContent(col, row, mainc, nil, st)
		}
	}
	s.Show()
}
