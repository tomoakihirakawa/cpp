set xlabel 'Expansion Order'
set ylabel 'Error'
set logscale y
set title 'FMM Potential Approximation Error'
plot 'fmm_error_data.dat' using 1:2 with linespoints title 'Error'
