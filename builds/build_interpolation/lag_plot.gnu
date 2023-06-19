set title "Lagrangian interpolation" font ",15"
set xlabel "x"
set ylabel "y"

set grid

plot  "lag_interpolation.dat" using 1:2 with lines title 'interplation', \
      "lag_exact.dat" using 1:2 with lines title 'exact', \
      "lag_data.dat" using 1:2 with points pointtype 7 pointsize 1.5 title 'Points', \
      "lag_data.dat" using 1:2:3 with labels offset char 1,1 font ",15" notitle
pause -1
