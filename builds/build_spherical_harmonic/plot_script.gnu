# set output format
set terminal pngcairo enhanced font "Times-New-Roman,10" size 400,400 

# set labels
set xlabel "x" font "Times-New-Roman:Italic,15"
set ylabel "y" font "Times-New-Roman:Italic,15"
set title "Error of \log_{10}(\|G_{approx}(x,a) - G\|/\|G\|)" font "Times-New-Roman,15"

# enable contouring
set pm3d
set view map
set contour base

# set automatic contour levels
set cntrparam levels auto 5

# set color palette
set palette rgb 33,13,10

# set colorbox range
set cbrange [-8:1]

# set colorbox and its label
set colorbox vertical user origin 0.9,0.2 size 0.02,0.6
set cblabel "Error" font "Times-New-Roman,12"

# set tics font
set xtics font "Times-New-Roman,10"
set ytics font "Times-New-Roman,10"

# plot commands for different n values
set output 'output_n4_A_10_10_10.png'
splot 'output_n4_A_10_10_10.txt' using 1:2:4 with pm3d

set output 'output_n5_A_10_10_10.png'
splot 'output_n5_A_10_10_10.txt' using 1:2:4 with pm3d

set output 'output_n6_A_10_10_10.png'
splot 'output_n6_A_10_10_10.txt' using 1:2:4 with pm3d

set output 'output_n7_A_10_10_10.png'
splot 'output_n7_A_10_10_10.txt' using 1:2:4 with pm3d

set output 'output_n8_A_10_10_10.png'
splot 'output_n8_A_10_10_10.txt' using 1:2:4 with pm3d

# plot commands for different n values
set output 'output_n4_A_5_5_5.png'
splot 'output_n4_A_5_5_5.txt' using 1:2:4 with pm3d

set output 'output_n5_A_5_5_5.png'
splot 'output_n5_A_5_5_5.txt' using 1:2:4 with pm3d

set output 'output_n6_A_5_5_5.png'
splot 'output_n6_A_5_5_5.txt' using 1:2:4 with pm3d

set output 'output_n7_A_5_5_5.png'
splot 'output_n7_A_5_5_5.txt' using 1:2:4 with pm3d

set output 'output_n8_A_5_5_5.png'
splot 'output_n8_A_5_5_5.txt' using 1:2:4 with pm3d


set title "Error of \log_{10}(\|grad G_{approx}(x,a) - grad G\|/\|grad G\|)" font "Times-New-Roman,15"

# plot commands for different n values
set output 'output_n4_A_10_10_10_grad.png'
splot 'output_n4_A_10_10_10_grad.txt' using 1:2:4 with pm3d

set output 'output_n5_A_10_10_10_grad.png'
splot 'output_n5_A_10_10_10_grad.txt' using 1:2:4 with pm3d

set output 'output_n6_A_10_10_10_grad.png'
splot 'output_n6_A_10_10_10_grad.txt' using 1:2:4 with pm3d

set output 'output_n7_A_10_10_10_grad.png'
splot 'output_n7_A_10_10_10_grad.txt' using 1:2:4 with pm3d

set output 'output_n8_A_10_10_10_grad.png'
splot 'output_n8_A_10_10_10_grad.txt' using 1:2:4 with pm3d

# plot commands for different n values
set output 'output_n4_A_5_5_5_grad.png'
splot 'output_n4_A_5_5_5_grad.txt' using 1:2:4 with pm3d

set output 'output_n5_A_5_5_5_grad.png'
splot 'output_n5_A_5_5_5_grad.txt' using 1:2:4 with pm3d

set output 'output_n6_A_5_5_5_grad.png'
splot 'output_n6_A_5_5_5_grad.txt' using 1:2:4 with pm3d

set output 'output_n7_A_5_5_5_grad.png'
splot 'output_n7_A_5_5_5_grad.txt' using 1:2:4 with pm3d

set output 'output_n8_A_5_5_5_grad.png'
splot 'output_n8_A_5_5_5_grad.txt' using 1:2:4 with pm3d
