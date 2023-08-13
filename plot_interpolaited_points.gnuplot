set datafile separator whitespace
set terminal qt size 800,600
set title "3D Points Plot"
set xlabel "X"
set ylabel "Y"
set zlabel "Z"
set grid
set key off
set view equal xyz
set xrange [-10:10]
set yrange [-10:10]
set zrange [-5:5]

# Указываем наименования файлов для каждого набора данных
interpolated_points = "real_points.txt"
aprox_points = "aprox_points.txt"
gavno_points = "GoraGovna_points.txt"

# Указываем, что у нас трехмерные данные (x, y, z)
splot gavno_points using 1:2:3 with points pointtype 7 pointsize 1.5 lc rgb "red" title "Interpolated Points"

