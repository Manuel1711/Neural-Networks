reset

set terminal postscript eps enhanced color dashed font 'Arial, 23'
set output 'graph.eps'

#-1-
# 0.5**6*1/2/(3/4)**(1/8)/6

#set terminal pngcairo font "Arial, 12"
#set terminal png font "Arial, 12"
#set output 'graph.png'

#set terminal 'jpeg' font arial 12 
#set output 'graph.jpeg'

#set for [i=1:9] linetype i dashtype i    {/Arial=20 u} = 0.01   {/Arial=20 uL^2} = 6000
#set style data linespoints

#set multiplot layout 2,2 rowsfirst title "{/:Bold=15 title}"\
              #margins screen MP_LEFT=0., MP_RIGHT=0., MP_BOTTOM=0., MP_TOP=0. spacing screen MP_xGAP=1., MP_yGAP=1.       # set xlabel after a plot to have not repetition

set multiplot

#set key spacing 1
set key top right          # top, bottom, right, left, center
#set key at first 7.7, first 0.3
#set key width 0.1 height 0.1
#set key horizontal 
#set key maxcols 4 #maxrows
#set key samplen 3        # length of line's type
set key font 'Arial, 27'

#set lmargin at screen 0.15#{{at screen} <margin>} #set rmargin #set tmargin #set bmargin #set margins
set lmargin at screen 0.14
set bmargin at screen 0.129
set rmargin at screen 0.965
#set tmargin at screen 1

#set xtics axis                    # start, incr, end
set mxtics 4                                                # number of mtics in 1 interval
set xtics scale 2.3, 1.                                    # dimension major, minor tics   1=default
#set xtics border offset 0,0.5 -5,1,5
#set xtics offset 0,0.5

set mytics 4                                               
set ytics scale 2.3, 1.

#set  border lt 1 ls 2 lc 1 lw 2 dt 2

#set grid
#set xtics 0.05
#set ytics 20

#set ytics add ("[+0.135]" 0.135) 
set xtics add ("0" 0) 
set ytics add ("0" 0)

set ytics font "Arial,20"
set xtics  font "Arial,20"

#set format x "%2.0t{/Symbol \264}10^{%L}"
#set format y "%2.0t{/Symbol \264}10^{%L}"
#set format x "10^{%L}"

#set format y "%.1f"   # number format in axis
#set format x "%1.0t {/Bold=18 ~ {.4\.}}10^{%L}"
#set format y "%1.0t {/Bold=18 ~ {.4\.}}10^{%L}"

set ylabel font 'Arial, 35'
set xlabel font 'Arial, 35'
#set format y "%.1f"   # number format in axis

#set grid 

#set title "Mean value of particles number N with {/Symbol=25 k}=1, {/Symbol=25 g}=1 and {/Symbol=25 d}=0"
#set title "Quench : {/Symbol=25 k} _i = 0 {/Symbol -->}  {/Symbol=25 k} = 3,   {/Symbol=25 g} = 100"
#set title " {/Arial=20 u} = 0.01 "
#set title "{/Symbol=25 k} = 0,  {/Symbol=25 g} = 20"#           {/Symbol=25 m}=1.8 and 
#set title "{/Symbol=25 d} = 1 and {/Symbol=25 m} = 2"
set xrange [0:0.2]
#set yrange [-2:60]
#set xlabel "{/Symbol=45 q}" offset 0.0,0.5
set xlabel "{/Symbol w}" offset 0, 0.4
#set xlabel "{/Symbol=25 k}"
#set xlabel 'w' offset 0,0.8
set ylabel "{/Symbol r}" offset 0.5, 0		# L^{/Arial=20 y_{/Symbol f}}
#set ylabel "{/Times=50 J}" offset 0.8,0
#set ylabel "L [ n_{L/2}({/Symbol=22 q}) - n_{L/2}(0) ]"
#set ylabel "L < C^+_{L/3} C_{2L/3} + h.c. >"
#set ylabel "~G{.6\\~}_c" offset 0.3, 0
#set ylabel "L G_c" offset 0.,0
#set ylabel "N / N_o" offset 0.8,0

#set logscale x 10
#set logscale y 10

#set arrow from -2, graph 0 to -2, graph 1 nohead lt 6 dt '_ '
#set arrow from graph 0, first 0.25 to graph 1, first 0.25 nohead lt 4 lw 2 dt '_._. ' 
#set arrow from graph 0, first -0.05  to graph 1, first -0.05 nohead lt 7 lw 1 dt '_._.'
#set arrow from graph 0, first -0.5  to graph 1, first -0.5 nohead lt 7 lw 1 dt '_._.'
#set arrow from first 0, first -1.6  to first 0.4, first -1.4 lt 7 lw 1 dt 1
#set arrow from first 0.4, first 0.8  to first 0., first 0.6 lt 7 lw 1 dt 1


set label font "Arial,30"

#set label 2 at first 0.01, first 0.55 "y = 0.5" tc lt 7 # graph, screen or first
#set label 3 at graph 0.93, first -0.6 "- 0.5"
#set label 4 at graph 0.955, first -1.1 "- 1"
#set label 1 at graph 0.05, first 0.45 "D@^{} _{}" tc lt 4
#set label 1 at graph 0.02, graph 0.09 "~{/Symbol m}{.9/Bold=35\\_} = 0" tc lt 8
#set label 2 at graph 0.85, graph 0.05 "{/Arial=30 w = 1}" tc lt 8

#set label 3 at graph 0.78, graph 0.76 "{/Arial=30 {/Symbol Q}_{/ZapfDingbats=15 \110} = -1.5}" tc lt 8
#set label 1 at graph 0.78, graph 0.22 "{/Arial=30 h_i = {/Arial=28 -0.05}}" tc lt 8
#set label 3 at graph 0.78, graph 0.76 "{/Arial=30 {/Symbol S} = 0}" tc lt 8
#set label 3 at graph 0.78, graph 0.76 "{/Arial=30 h = 0.2}" tc lt 8
#set label 2 at graph 0.78, graph 0.86 "{/Arial=30 {/Symbol U} = 0.05}" tc lt 8

plot  "outputNNwidth16sec.dat" every 200 using ($1):($2*90**(0.)):($3*50) t 'w = 16    ' w yerrorlines pt 11 ps 1 lw 3. lt rgb 'blue' dt 8 pi 100,\
      "outputNNwidth512sec.dat" every 200 using ($1):($2*90**(0.)):($3) t 'w = 512  ' w yerrorlines pt 3 ps 1 lw 3. dt 3 pi 100 lt 1,\
      "test.ddat" using ($1):(($2)*30**(0.)) t 'smear' w l lw 3. lt rgb 'dark-green' dt 2,\
      "outputNNwidth1024sec1k_step.dat" every 200 using ($1):($2*90**(0.)):($3) t 'w = 1024' w yerrorlines pt 1 ps 1 lw 3. dt 5 pi 100 lt 6,\
      "outputNNwidth2048sec1k_step.dat" every 200 using ($1):($2*70**(0.)):($3) t 'w = 2048' w yerrorlines pt 7 ps 1 lw 3. lt rgb 'red' dt 4 pi 100,\
      "outputNNwidth4096sec.ddat" every 200 using ($1):($2*70**(0.)):($3) t 'w = 4096' w yerrorlines pt 9 ps 1 lw 3. lt rgb 'dark-green' dt 4 pi 100,\
      "data_sec0611ntk1000.dat" every 200 using ($1):($2*90**(0.)):($3) t 'NTKgp' w yerrorlines pt 5 ps 1 lw 3. lt 8 dt 1 pi 100,\
      "data_sec0611nngp1000.dat" every 200 using ($1):($2*70**(0.)):($3) t 'NNGP' w yerrorlines pt 7 ps 1 lw 1. lt rgb 'dark-green' dt 4 pi 100,\

# w linespoint pt 6 ps 2 pi 30 lt 1 dt '-' lw 2,\


#      "py.dat" using ($1):(($2)*30**(1/8.)) t 'L = 150' w l lw 3. lt 8 dt 1,\
#      "isC2DeS150l50.dat" using ($1):($2*50**(1./8.)) t 'L = 50  ' w l lw 5. lt 1 dt 4,\
       #      "isC2DeS150l250.dat" using ($1):($2*250.**(1./8.)) t 'L = 250' w l lw 2. lt rgb 'black' dt 6,\
#      "isC2DeS150l300.dat" using ($1):(($2)*300.**(1/8.)) t 'L = 300' w l lw 5. lt 8 dt 7,\
#      "data/tPyL8g1.dat" using ($1):($2) t 'L=4' pt 7 ps 0.8 lt 7,\
#      "data/CtL6g1.dat" using ($1):($2) t 'L=4' w l lw 5. lt rgb 'blue' dt 1 ,\
#      "data/kzUp1Sig-1L8.dat" using ($1):($2) t 'L=4' w l lw 5. lt rgb 'dark-green' dt 2 ,\
#      "data/kzUp1Sig-1L8.dat" using ($1):($2) t 'L=4' w l lw 5. lt 1 dt 6 ,\
#      "data/gsL6g1.dat" using ($1):($2*6.**(0.)) t 'L=6' w l lw 5. lt rgb 'orange' dt 5 ,\
#      "data/kzUp1Sig-1L7.dat" using ($1/7.):($2*7.**(0.)) t 'L=7' w l lw 5. lt rgb 'red' dt 6,\
#      "data/kzUp1Sig-1L8.dat" using ($1/8.):($2*8.**(0.)) t 'L=8' w l lw 5. lt 8 dt 1,\
#do for [N=1:5] {
#plot func(N, x)
#pause -1
#}

#set size 0.525,0.445
#set origin 0.454,0.12
unset margin

set size 0.525,0.4
set origin 0.45,0.58

unset xlabel
unset ylabel
unset title

set ytics font "Arial,23"
set xtics  font "Arial,23"
#set xtics 0.02#, _, _  # offset _, _
unset key
#set xrange
unset label
unset xrange
unset yrange

set ytics 0.2
set xtics 0.01

set xlabel font "Arial,26"
#set ylabel font "Arial,27"
set label font "Arial,25"

set ytics font "Arial,20"
set xtics  font "Arial,20"

#set format x "%.3f"
#set yrange [0.76:0.]
set xrange [0.00:0.018]
set xlabel "1/L" offset 6, 1.8
#set ylabel "D({/Arial=9 L/2, t}) - D({/Arial=9 L/2, 0})" offset 3.5,0
#set xlabel "t L^{- 3}" offset 0,0.8
#set ylabel "L < C^+_{3L/8} C_{5L/8} + h.c. >" offset 2.5,0

set label 1 at screen 0.82, screen 0.87 "{/Symbol=28 W} = 1"

f(x)=a*x + b 
#fit f(x) "estr-isC2DS15Y104.dat" u (1./$1):($2*$1**(1./8.)) via a,b


#plot  f(x) notitle w l dashtype '_' lt 1, "estr-isC2DS15Y104.dat" using (1./$1):($2*$1**(1./8.)) t '' pt 7 ps 2 # lt 8
#plot  "k-dissipationL30.dat" using 1:2 notitle

#set label 1 at screen 0.685, screen 0.58 "Asymptotic behavior" tc lt 4
#set arrow from graph 0, first -0.5 to graph 1, first -0.5 nohead lt 3 lw 2 dt ' _ _ '

#plot  "k-dissipationL25.dat" using (1./$1):2 notitle w l lw 1.0 lt 8 dt '- - ',\
#      "k-dissipationL20.dat" using (1./$1):2 notitle w l lw 1.0 lt 7 dt '- - ',\

#plot  "k-dissipationL20.dat" using 1:2 t 'L = 16,  {/Arial=20 u} = 0.01' w l lw 1.0 lt 8 dt 3,\
#      "k-dissipationL25.dat" using 1:2 t 'L = 32,  {/Arial=20 u} = 0.01' w l lw 1.0 lt 7 dt 3,\

unset multiplot
