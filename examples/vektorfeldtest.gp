
set terminal postscript enhanced color eps
set output 'vektorfeldtest.eps';
plot 'vektorfeldtest.dat' using 2:3:($4/5):($5*30) with vector

