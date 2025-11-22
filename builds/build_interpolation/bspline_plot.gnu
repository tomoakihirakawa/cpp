# --------------------------------- plot 1 --------------------------------- */

input1 = "bspline_sample_data.dat"
input2 = "bspline_interpolated_data.dat"

set title "Bspline interpolation" font ",15"
set xlabel "x"
set ylabel "y"

set key font ",15"

set style line 2 linecolor rgb '#0060ad' linetype 2 linewidth 2
set style line 3 linecolor rgb '#d160ad' linetype 2 linewidth 2

set grid

plot  input1 using 1:2 with lines linestyle 1 title 'sample data',\
      input2 using 1:2 with lines linestyle 1 title 'interpolated',\
      input2 using 1:2 with points pointtype 7 pointsize 1. title 'interpolated',\
      input2 using 1:3 with lines linestyle 1 title 'interpolated',\
      input2 using 1:3 with points pointtype 7 pointsize 1. title 'interpolated'

pause -1

# --------------------------------- plot 2 --------------------------------- */

bodyA = "../build_pybind11/bodyA.dat"
bodyB = "../build_pybind11/bodyB.dat"
bodyC = "../build_pybind11/bodyC.dat"

intp_bodyA = "bspline_interpolated_bodyA.dat"
intp_bodyB = "bspline_interpolated_bodyB.dat"
intp_bodyC = "bspline_interpolated_bodyC.dat"

set xrange [0:1.5]
set yrange [-0.15:0.15]

plot  bodyA using 2:3 with points pointtype 7 pointsize 1. title 'original bodyA',\
      bodyB using 2:3 with points pointtype 7 pointsize 1. title 'original bodyB',\
      bodyC using 2:3 with points pointtype 7 pointsize 1. title 'original bodyC',\
      intp_bodyA using 2:3 with lines linestyle 1 title 'intp bodyA',\
      intp_bodyB using 2:3 with lines linestyle 1 title 'intp bodyB',\
      intp_bodyC using 2:3 with lines linestyle 1 title 'intp bodyC'

pause -1
