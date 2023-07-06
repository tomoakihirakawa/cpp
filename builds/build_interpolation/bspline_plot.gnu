set title "Bspline interpolation" font ",15"
set xlabel "x"
set ylabel "y"

set key font ",15"

set style line 2 linecolor rgb '#0060ad' linetype 2 linewidth 2
set style line 2 linecolor rgb '#d160ad' linetype 2 linewidth 2

set grid

plot  "bspline_sample_data.dat" using 1:2 with lines linestyle 1 title 'sample data',\
      "bspline_interpolated_data.dat" using 1:2 with lines linestyle 1 title 'interpolated',\
      "bspline_interpolated_data.dat" using 1:2 with points pointtype 7 pointsize 1. title 'interpolated',\
      "bspline_interpolated_data.dat" using 1:3 with lines linestyle 1 title 'interpolated',\
      "bspline_interpolated_data.dat" using 1:3 with points pointtype 7 pointsize 1. title 'interpolated'
pause -1
