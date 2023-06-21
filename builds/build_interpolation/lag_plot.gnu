set title "Lagrangian interpolation" font ",15"
set xlabel "x"
set ylabel "y"

set key font ",15"

set style line 2 linecolor rgb '#0060ad' linetype 2 linewidth 2
set style line 2 linecolor rgb '#d160ad' linetype 2 linewidth 2

set grid

plot  "lag_exact_interpolation_derivative.dat" using 1:2 with lines linestyle 1 title 'exact sin', \
      "lag_exact_interpolation_derivative.dat" using 1:3 with lines linestyle 2 title 'exact cos', \
      "lag_exact_interpolation_derivative.dat" using 1:4 with lines title 'interplation sin', \
      "lag_exact_interpolation_derivative.dat" using 1:5 with lines title 'interplation D(sin)', \
      "lag_data.dat" using 1:2 with points pointtype 7 pointsize 1.5 title 'Points', \
      "lag_data.dat" using 1:2:3 with labels offset char 1,1 font ",15" notitle
pause -1
