set datafile separator whitespace
set terminal qt size 800,600
set title "3D Points Plot"
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set grid
set key off
set view equal xyz
set xrange [0:5]
set yrange [0:5]
set zrange [0:5]

# Указываем наименования файлов для каждого набора данных
interpolated_file = "test.txt"

# Указываем, что у нас трехмерные данные (x, y, z)
splot interpolated_file using 1:2:3 with points pointtype 7 pointsize 1.5 lc rgb "red" title "Interpolated Points"

