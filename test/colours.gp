set term pdf size 9,8 font "Bitstream Charter,12"
set output "pault-colour.pdf"

unset title
unset key
unset border
unset tics

set multiplot
# plot no. 1 ################################
set lmargin at screen 0
set rmargin at screen 1
set tmargin at screen 1
set bmargin at screen 1./2

set parametric
set samples 9
set trange [-pi*1./8:pi*15./8]

a = sqrt(2)/(2*cos(pi/8))
b = 0.07

set xrange [-1:8]
set yrange [-1:3]
set size ratio 4./9

set style fill solid 1.0 noborder
set style function filledcurves

set label 1 "51,34,136\n#88CCEE"   at 0, 0+b center front
set label 2 "68,170,153\n#44AA99"  at 1, 1+b center front
set label 3 "136,204,238\n#88CCEE" at 0, 2+b center front
set label 4 "17,119,51\n#117733"   at 2, 0+b center front
set label 5 "153,153,51\n#999933"  at 3, 1+b center front
set label 6 "221,204,119\n#DDCC77" at 4, 2+b center front
set label 7 "204,102,119\n#CC6677" at 5, 1+b center front
set label 8 "136,34,85\n#882255"   at 6, 0+b center front
set label 9 "170,68,153\n#AA4499"  at 7, 1+b center front

plot a*cos(t)     , a*sin(t)      lc rgb '#332288', \
     a*cos(t) +  1, a*sin(t) + 1  lc rgb '#44aa99', \
     a*cos(t)     , a*sin(t) + 2  lc rgb '#88ccee', \
     a*cos(t) +  2, a*sin(t)      lc rgb '#117733', \
     a*cos(t) +  3, a*sin(t) + 1  lc rgb '#999933', \
     a*cos(t) +  4, a*sin(t) + 2  lc rgb '#ddcc77', \
     a*cos(t) +  5, a*sin(t) + 1  lc rgb '#cc6677', \
     a*cos(t) +  6, a*sin(t)      lc rgb '#882255', \
     a*cos(t) +  7, a*sin(t) + 1  lc rgb '#aa4499'

# plot no. 2 ################################
unset label
set samples 45
set isosamples 20
set urange [0:pi]
set vrange [-pi:pi]
set xrange [-3.2:3.2]
set yrange [-3.2:3.2]
set size square
set view map
set pm3d interpolate 5,5
set style line 1 lc rgb "black" lw 1.5

set tmargin at screen 1./2
set bmargin at screen 0.01
set lmargin at screen 0.03
set rmargin at screen 0.48

rcol(x) = (0.472-0.567*x+4.05*x**2)/(1.+8.72*x-19.17*x**2+14.1*x**3)
gcol(x) = 0.108932-1.22635*x+27.284*x**2-98.577*x**3+163.3*x**4-131.395*x**5+40.634*x**6
bcol(x) = 1./(1.97+3.54*x-68.5*x**2+243*x**3-297*x**4+125*x**5)
set palette model RGB functions rcol(gray), gcol(gray), bcol(gray)
splot v*sin(u), v*cos(u), sqrt(6) * (v**2 * sin(2*u)) w l ls 1

# plot no. 3 ################################
set lmargin at screen 0.52
set rmargin at screen 0.97

rcol(x) = 0.237 - 2.13*x + 26.92*x**2 - 65.5*x**3 + 63.5*x**4 - 22.36*x**5
gcol(x) = ((0.572 + 1.524*x - 1.811*x**2)/(1 - 0.291*x + 0.1574*x**2))**2
bcol(x) = 1/(1.579 - 4.03*x + 12.92*x**2 - 31.4*x**3 + 48.6*x**4 - 23.36*x**5)
set palette model RGB functions rcol(gray), gcol(gray), bcol(gray)
splot v*sin(u), v*cos(u), sqrt(6) * (v**2 * sin(2*u)) w l ls 1

unset multiplot

set output

