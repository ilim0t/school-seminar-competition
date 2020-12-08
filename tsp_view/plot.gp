#!/usr/bin/gnuplot

set terminal postscript
set out "image.eps"
set multiplot

plot "data/points.dat" with points lw 0.2 pt 6 ps 0.5, "data/lines.dat" with lines lt 1
# plot "data/points.dat" with points lw 0.2 pt 6 ps 0.5, "data/lines.dat" with lines lt 5