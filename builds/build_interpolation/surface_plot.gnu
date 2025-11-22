# $ gnuplot surface_plot.gnu

set title "Triangular Surface Plot" font ",15"
set xlabel "x"
set ylabel "y"
set zlabel "z"

set grid

splot "triangular_surface.dat" using 1:2:3 with lines title 'Surface', \
      "points.dat" using 1:2:3 with points pointtype 7 pointsize 1.5 title 'Points', \
      "points.dat" using 1:2:3:4 with labels offset char 1,1 font ",15" notitle
pause -1
