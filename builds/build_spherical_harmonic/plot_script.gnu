# set output format
set terminal pngcairo

# set labels
set xlabel "X"
set ylabel "Y"
set title "Error of G_approx from G"

# enable contouring
set pm3d
set view map
set contour base

# set automatic contour levels
set cntrparam levels auto 10

# set color palette
set palette rgb 33,13,10

# set colorbox range
set cbrange [-20:0]

# plot commands for different n values
set output 'output_n3_A_10_10_10.png'
splot 'output_n3_A_10_10_10.txt' using 1:2:4 with pm3d

set output 'output_n6_A_10_10_10.png'
splot 'output_n6_A_10_10_10.txt' using 1:2:4 with pm3d

set output 'output_n9_A_10_10_10.png'
splot 'output_n9_A_10_10_10.txt' using 1:2:4 with pm3d

# plot commands for different n values
set output 'output_n3_A_5_5_5.png'
splot 'output_n3_A_5_5_5.txt' using 1:2:4 with pm3d

set output 'output_n6_A_5_5_5.png'
splot 'output_n6_A_5_5_5.txt' using 1:2:4 with pm3d

set output 'output_n9_A_5_5_5.png'
splot 'output_n9_A_5_5_5.txt' using 1:2:4 with pm3d
