;Johanna Blee model to simulate biofilms exposed to constant or perioidic antibiotic treatemnts
;Use the interface and sliders to change the values required
globals [
  itdiff
  t                 ;global time (s)
  Ov                ;used in the shoving algorithm (total overlapping cells)
  h                 ;maximal biofilm height
  mS mA mDs mDa mdv mdh mM mMD mBulk-yes mBulk-no  ;matrices used in the diffusion-reaction algorithm
  biovol            ;volume of the biofilm (sum of the volume of all cells' types, µm3/µm²)
  mean-thickness    ;mean thickness of the biofilm
  Ra                      ;roughtness coefficient of the biofilm
  antibiotic-activated    ;true when there is antibiotic
  up-persisters           ;number of persisters with ycor > mean-thickness/2
  down-persisters         ;number of persisters under with ycor < mean-thickness/2
  i;looping variable
  tottime;total time exposed to antibiotic
  ]
patches-own [
  S     ;substrate concentration
  A     ;antibiotic concentration
  Ds    ;substrate diffusion coefficient
  Da    ;antibiotic diffusion coefficient
  bulk  ;true if the patch belongs to the bulk
  ]
turtles-own [
  m               ;mass of the cell
  dv xv yv        ;used to store the shoving vector
  lineage         ;used to follow the lineage of the initial bacteria
  ]
breed [bacteria bacterium]     ;susceptible cells
breed [persisters persister]   ;persisters
breed [dead dead1]             ;dead cells
extensions [profiler matrix]

to setup
  clear-all
  ask patches [
    set S Cs-bulk set Ds D-substrate set Da D-antibiotic
    set pcolor white]
  if delta-t * D-substrate / (delta-l * delta-l) >= 0.25 [error "dt*D/(dl*dl) too hight !"]
  set start-treatment? false
  setup-bacteria
  bulk-position
  vert
  hori
  set-mS
  set-mDs
  set-mA
  set-mDa
  set-mM
  set-mMD
  reset-ticks
end

to go
   ;;; --- SIMULATION PROCESS --- ;;;

   treatment
   set itdiff 1
   repeat (delta-tb / delta-t) [diffusion-reaction-substrate if antibiotic-activated = true [diffusion-reaction-antibiotic]] ;diffusion-reaction operates at a different time step than bacteria
   ask patches [set S matrix:get mS pycor pxcor if S < 0 [error "Substrate concentration negative! Decrease delta-t or decrease µ."]] ;stop the simulation if the mass balance is not respected
   if antibiotic-activated = true [
   ask patches [set A matrix:get mA pycor pxcor]
   if mean [A] of patches < (MIC / 10000) [set antibiotic-activated false]] ;deactivate the procedure diffusion-reaction-antibiotic if the antibiotic concentration is small enough
   coloring
   bact-update
   bulk-position
   set-mDs
   set-mDa
   phenotypic-switches

   ;;; --- INDICATORS --- ;;;

   biovolume
   thickness-distribution
   persisters-position


   set t t + delta-tb ; update time (s)



  ;;; --- STOP CONDITIONS & EXPORTED DATA --- ;;;

  ;uncomment any required code to export images of the simulation

   ;if Video? [if ((ceiling (ticks / 10)) * 10) = ticks [export-view (word ticks ".png")]] ;uncommment to export an image of the view every 10 ticks, can be used to make a video (no GUI needed)
   ;export-view (word behaviorspace-run-number "images" ticks ".png") ]

   ;export data at the start of the antibiotic treatment
   ;if t = (auto-start-treatment * 3600) [
     ;export-world (word "World " auto-start-treatment "h " behaviorspace-run-number ".csv")
     ;export-view (word "View " auto-start-treatment "h " behaviorspace-run-number ".png")
     ;color-lineage export-view (word "Lineage " auto-start-treatment "h " behaviorspace-run-number ".png")]
   ;file-open "/Users/tt21663/OneDrive - University of Bristol/run-outputs.txt"
   ;file-print (word "---------- Tick Number: " ticks "-----------")
   ;file-print (ticks)
   ;file-print (a-max)]

   ;export data at the end of the antibiotic treatment
   ;if t = ((auto-start-treatment + treatment-duration) * 3600) [
     ;export-world (word "World " (auto-start-treatment + treatment-duration) "h " behaviorspace-run-number ".csv")
     ;export-view (word "View " (auto-start-treatment + treatment-duration) "h " behaviorspace-run-number ".png")
     ;color-lineage export-view (word "Lineage " (auto-start-treatment + treatment-duration) "h " behaviorspace-run-number ".png")]

   ;export data at the end of the recovery and stop the simulation
   ;if t > ((auto-start-treatment + treatment-duration + recovery-time) * 3600) [
     ;export-world (word "World end " (auto-start-treatment + treatment-duration + recovery-time) "h " behaviorspace-run-number ".csv")
     ;export-view (word "View end " (auto-start-treatment + treatment-duration + recovery-time) "h " behaviorspace-run-number ".png")
     ;color-lineage export-view (word "Lineage end " (auto-start-treatment + treatment-duration + recovery-time) "h " behaviorspace-run-number ".png")
     ;stop]

   ;stop the simulation if the bulk leaves the computational domain (increase the number of patches to avoid this)
   if (boundary-layer-lenght / delta-l + h) > max-pycor [
     ;export-world (word "World top " behaviorspace-run-number ".csv")
     export-view (word "View top " behaviorspace-run-number ".png")
     ;color-lineage export-view (word "Lineage top " behaviorspace-run-number ".png")
     stop]

   ;stop the simulation if there is not any live cell (susceptible or persister)
   if (count bacteria + count persisters) = 0 [
     ;export-world (word "World dead " behaviorspace-run-number ".csv")
     export-view (word "View dead " behaviorspace-run-number ".png")

    stop]

     ;color-lineage export-view (word "Lineage dead" behaviorspace-run-number ".png") stop]
;a-max b-max treatment-duration recovery-time auto-start-treatment
   tick ;tick every bacteria update
end


    ;;; --- BACTERIA --- ;;;

to setup-bacteria ;set up an initial number of bacteria on the surface with random positions
  create-bacteria initial-bacteria
  [setxy random-xcor (-0.5 + 1 / delta-l / 2)
    set shape "circle"
    set color green
    set m (m-div / 2 + random m-div / 2)
    set size ((4 * m / (pi * cell-density * delta-l)) ^ 0.5 / delta-l)
    set lineage who]
  repeat 10 [shoving]
end

to bact-update ;division, size, death and lysis of cells
  ask bacteria [
    set m m * exp(µmax * (S / (S + Ks)) / 3600 * delta-tb)
    if m > m-div [let div random-float 0.1 rt random 360 hatch 1 [fd (size / 2) set m m * (0.5 + div)] bk (size / 2) set m (m * 0.5 - div)]
    set size 2 * (m / (pi * cell-density * delta-l)) ^ 0.5 / delta-l
    if matrix:get mA pycor pxcor > MIC and random-float 1 < (k-max * matrix:get mA pycor pxcor / (matrix:get mA pycor pxcor + Kk * MIC) / 3600 * delta-tb) [hatch-dead 1 [set shape "circle" set color 4] die]
  ]
  ask persisters [
    if matrix:get mA pycor pxcor > MIC and random-float 1 < (k-max-pers * matrix:get mA pycor pxcor / (matrix:get mA pycor pxcor + Kk * MIC) / 3600 * delta-tb) [hatch-dead 1 [set shape "circle" set color 4] die]
  ]
  ask dead [
    set m m * exp(- lysis / 3600 * delta-tb)
    set size 2 * (m / (pi * cell-density * delta-l)) ^ 0.5 / delta-l * 0.8 ;dead cells have a reduced size to acount for their shrinking and to relieve the shoving algorithm
    if m < 0.05 * m-div [die]
  ]
  shoving
  set-mM ;update the mass of active bacteria per patch
  set-mMD ;update the mass of dead bacteria per patch
end

     ;;; --- SHOVING --- ;;;

to calculate-Ov ;Number of overlapping cells
  set Ov 0
  ask turtles [
    if (count turtles in-radius (2.5 * size) with [distance myself > 0 and distance myself / ([size] of myself / 2 + [size] of self / 2) < 0.99 ]) > 0 ;count cells that overlap more than 1% of their size and that are not themselves
      [set Ov Ov + 1
       ]
  ]
end

to shoving
  let it 1
  calculate-Ov
  while [(Ov / count turtles * 100) > max-moving-cells] [ ;if Ov represents more than a particular percentage the algorithm continues
    ;;; --- Shoving
    repeat 10 [
      ask turtles [
        ask turtles in-radius (2.5 * size) ;select cells in my surroudings
        with [distance myself - [size] of myself / 2 - [size] of self / 2 < 0 and distance myself > 0] ;select overlapping cells that are not me
        [set xv xv + [xcor] of self - [xcor] of myself ;keep in memory the relative coordinates where the cell must head (direction of the shoving vector)
         set yv yv + [ycor] of self - [ycor] of myself
         set dv dv + (size / 2 + [size] of myself / 2 - distance myself)]] ;keep in memory the distance to move (lenght of the vector)
      ask turtles [facexy (xcor + xv) (ycor + yv) fd (dv / 2) ;dv/2 because each overlapping cell will move in opposite directions
        set dv 0 set xv 0 set yv 0]]
    ;;; --- Update Ov
    calculate-Ov
    set it it + 1
    if it > 100 [stop] ;stop the algorithm above 1000 iterations (the simulation keeps going)
  ]
end

      ;;; --- BULK --- ;;;

to bulk-position
  set h max [ycor] of turtles ;maximum biofilm thickness (*patches)
  ask patches [
    ifelse count turtles-here > 0 [set Ds D-substrate * D* set Da D-antibiotic * D*] [set Ds D-substrate set Da D-antibiotic] ;update diffusion rates
    ifelse pycor > (boundary-layer-lenght / delta-l + h) [set bulk 1] [set bulk 0] ;update bulk position
    ]
  set-mBulk
end

to set-mBulk ; create a matrix mBulk with bulk of patches
  set mBulk-yes matrix:make-constant (max-pycor + 1) (max-pxcor + 1) 0 ;bulk set to 1
  set mBulk-no matrix:make-constant (max-pycor + 1) (max-pxcor + 1) 0 ;bulk set to 0
  let x 0
  let y 0
  while [x <= max-pxcor] [
    while [y <= max-pycor] [
      matrix:set mBulk-yes y x [bulk] of patch x y
      matrix:set mBulk-no y x [abs (bulk - 1)] of patch x y
      set y y + 1]
    set y 0
  set x x + 1]
end

      ;;; --- ANTIBIOTIC TREATMENT --- ;;;

to treatment
  if t / 3600 > auto-start-treatment [set start-treatment? true]
  set i 1

  repeat 10000[
   set tottime 60 * ( (t / 3600) - auto-start-treatment - ((i - 1) * recovery-time))
   if t / 3600 > auto-start-treatment + i * treatment-duration + (i - 1) * recovery-time [set start-treatment? false]
   if t / 3600 > auto-start-treatment + i * treatment-duration + i * recovery-time [set start-treatment? true]
   if start-treatment? = true [set antibiotic-activated true]
    if start-treatment? = true [ask patches [set pcolor orange]]
   set i (i + 1)

  ]

end

      ;;; --- DIFFUSION-REACTION --- ;;;

to diffusion-reaction-substrate
  let convergence sum map [ ?1 -> sum ?1 ] matrix:to-row-list mS

  let m1 map [ ?1 -> map [ ??1 -> ??1 / (??1 + Ks) ] ?1 ] matrix:to-row-list mS ; computes S/(S+Ks) for the growth rates {µ=µmax*m*S/(S+Ks)}
  set m1 matrix:from-row-list m1 ; converts m1 into a matrix

  set mS matrix:times-element-wise mBulk-no (mS ; updates substrate concentrations during dt
  matrix:+ (delta-t / (delta-l * delta-l)) matrix:* (mdv matrix:* (matrix:times-element-wise mDs mS) matrix:+ (matrix:times-element-wise mDs mS) matrix:* mdh) ; diffusion
  matrix:- (matrix:times-element-wise (delta-t / (delta-l ^ 3)) (µmax / 3600) (1 / Yxs) mM m1)) ; substrate consumption by active cell (1fg/µm3 = g/L)
  matrix:+ (matrix:times-element-wise (delta-t / (delta-l ^ 3)) (lysis / 3600) mMD) ; substrate release by dead cells
  matrix:+ (matrix:times-element-wise mBulk-yes Cs-bulk) ;set the bulk to Cs-bulk

  ;set mM mM matrix:+ (matrix:times-element-wise delta-t (µmax / 3600) mM m1) ; updates the susceptible biomass (it doesn't update the cells!)
  ;set mMD mMD matrix:- (matrix:times-element-wise delta-t (lysis / 3600) mMD) ; update the dead biomass (it doesn't update the cells!)

  if convergence / sum map [ ?1 -> sum ?1 ] matrix:to-row-list mS <= 1 + 10 ^ -9 [stop] ;stop the algorithm when the solution has converged (negligible variations)
  if convergence / sum map [ ?1 -> sum ?1 ] matrix:to-row-list mS <= 1 - 10 ^ -9 [stop]

  set itdiff itdiff + 1
end

to diffusion-reaction-antibiotic
  set mA matrix:times-element-wise mBulk-no (mA
      matrix:+ (delta-t / (delta-l * delta-l)) matrix:* (mdv matrix:* (matrix:times-element-wise mDa mA) matrix:+ (matrix:times-element-wise mDa mA) matrix:* mdh)) ; diffusion
  if start-treatment? = true [set mA mA matrix:+ (matrix:times-element-wise mBulk-yes (Ca-bulk * MIC))] ; keeps the concentration of antibiotic in the bulk constant during treatments
end

      ;;; --- PHENOTYPIC SWITCHES --- ;;;

to phenotypic-switches
if strategy = 1 [ ;stochastic switches
  ask bacteria [
    if random-float 1 < (a-max / 3600 * delta-tb) [hatch-persisters 1 [set shape "circle" set color red] die]
    ]
  ask persisters [
    if random-float 1 < (b-max / 3600 * delta-tb) [hatch-bacteria 1 [set shape "circle" set color green] die]
    ]
]

if strategy = 2 [ ;substrate-dependent switches
  ask bacteria [
    if random-float 1 < (a-max / 3600 * delta-tb * (1 - S / (S + K))) [hatch-persisters 1 [set shape "circle" set color red] die]
      ]
  ask persisters [
    if random-float 1 < (b-max / 3600 * delta-tb * (S / (S + K))) [hatch-bacteria 1 [set shape "circle"] die]
  ]
]

if strategy = 3 [ ;antibiotic-dependent switching rates - Strategy 3
  ask bacteria [
    if random-float 1 < (a-max / 3600 * delta-tb * (A / (A + K'))) [hatch-persisters 1 [set shape "circle" set color red] die]
    ]
  ask persisters [
    if random-float 1 < (b-max / 3600 * delta-tb * (1 - A / (A + K'))) [hatch-bacteria 1 [set shape "circle"] die]
  ]
]

if strategy = 4 [ ;antibiotic and nutrient dependent switching rates - Strategy 4 - one used in our study
  ask bacteria [
    if random-float 1 < (a-max / 3600 * delta-tb * (A / (A + K'))) [hatch-persisters 1 [set shape "circle" set color red] die]
    ]
  ask persisters [
    if random-float 1 < (b-max / 3600 * delta-tb * (1 - A / (A + K'))) [hatch-bacteria 1 [set shape "circle"] die]
    ]
  ask bacteria [
    if random-float 1 < (a-max / 3600 * delta-tb * (1 - S / (S + K))) [hatch-persisters 1 [set shape "circle" set color red] die]
     ]
  ask persisters [
    if random-float 1 < (b-max / 3600 * delta-tb * (S / (S + K))) [hatch-bacteria 1 [set shape "circle"] die]
    ]
]
end

      ;;; --- COLORING --- ;;;

to color-lineage
  ask bacteria [set color 25 + lineage * 10]
end

to coloring
  if show-s? [ask patches [set S matrix:get mS pycor pxcor set pcolor scale-color blue (S / (Cs-bulk + 0.0001)) 1.7 0]]
  ;f show-ATB? [ask patches [set A matrix:get mA pycor pxcor set pcolor orange]
  if color-µ? [ask bacteria [set color scale-color green (m * µmax / 3600 * matrix:get mS pycor pxcor / (matrix:get mS pycor pxcor + Ks)) ((m-div * µmax / 3600 * Cs-bulk / (Cs-bulk + Ks)) + 0.5) -0.1]]
end


      ;;; --- INDICATORS --- ;;;

to biovolume ; calculate the total volume of cells, indicator "biovol"
  set biovol count turtles * pi * (sum [(size / delta-l / 2) ^ 2] of turtles) / (delta-l ^ 2 * max-pxcor) ; total volume of cells / total bottom area (µm3/µm²)
end

to thickness-distribution ;report the mean thickness and the roughness coefficient of the biofilm, indicators "mean-thickness" and "Ra"
  let mat matrix:make-constant 1 (max-pxcor * 2 + 1) 0
  let x 0
  while [x <= max-pxcor * 2] [
    carefully [matrix:set mat 0 x (((max [ycor] of turtles with [xcor >= (x / 2 - 0.5) and xcor <= (x / 2)]) + 0.5) * delta-l)]
    []
    set x x + 1]
  set mean-thickness mean matrix:get-row mat 0
  let mylist map [ ?1 -> abs (?1 - mean-thickness) ] matrix:get-row mat 0
  set Ra mean mylist
end

to persisters-position
  set up-persisters count persisters with [(ycor * delta-l) > (mean-thickness / 2)]
  set down-persisters count persisters with [(ycor * delta-l) < (mean-thickness / 2)]
end

     ;;; --- RESTORE DEFAULT PARAMETERS'VALUES --- ;;;

to restore-default-values
  set initial-bacteria 27
  set Cs-bulk 0.4                ;g/L
  set boundary-layer-lenght 40   ;µm
  set delta-l 4                  ;µm
  set delta-t 0.0035             ;s
  set D* 0.8

  set µmax 1.25           ;h-1
  set Ks 0.0035           ;g/L
  set Yxs 0.21            ;ratio
  set cell-density 200    ;fg/µm3
  set m-div 500           ;fg
  set delta-tb 67         ;s
  set lysis 0.05          ;h-1

  set max-moving-cells 5  ;%

  set start-treatment? false
  set auto-start-treatment 8  ;h

  set treatment-duration   2  ;h
  set D-antibiotic 900        ;µm2/s
  set MIC 0.000005            ;g/L
  set k-max 10                ;h-1
  set k-max-pers 0.1          ;h-1
  set Kk 6.4                  ;*MIC

  set a-max 1.0        ;h-1
  set b-max 0.1        ;h-1
  set K Ks             ;g/L
  set K' MIC           ;g/L
end

     ;;; --- MATRICES USED FOR DIFFUSION-REACTION --- ;;;

to vert ;vertical diffusion matrix
  ; boundaries y=0 and y=ymax --> no flux (surface and top of the computational domain)
  set mdv matrix:make-constant (max-pycor + 1) (max-pycor + 1) 0
  matrix:set mdv 0 0 -1
  matrix:set mdv 0 1 1
  matrix:set mdv max-pycor max-pycor -1
  matrix:set mdv max-pycor (max-pycor - 1) 1
  let y 1
  while [y <= (max-pycor - 1)] [
    matrix:set mdv y y -2
    matrix:set mdv y (y + 1)  1
    matrix:set mdv y (y - 1)  1
    set y y + 1]
end

to hori ;horizontal diffusion matrix
  ; boundaries x=0 and x=xmax --> cyclic (sides of the computational domain)
  set mdh matrix:make-constant (max-pxcor + 1) (max-pxcor + 1) 0
  matrix:set mdh 0 0 -2
  matrix:set mdh 0 1 1
  matrix:set mdh 0 max-pxcor 1
  matrix:set mdh max-pxcor 0 1
  matrix:set mdh max-pxcor (max-pxcor - 1) 1
  matrix:set mdh max-pxcor max-pxcor -2
  let x 1
  while [x <= (max-pxcor - 1)] [
    matrix:set mdh x x -2
    matrix:set mdh x (x + 1)  1
    matrix:set mdh x (x - 1)  1
    set x x + 1]
end

to set-mS ; create a matrix mS with S of patches
  set mS matrix:make-constant (max-pycor + 1) (max-pxcor + 1) 0
  let x 0
  let y 0
  while [x <= max-pxcor] [
    while [y <= max-pycor] [
      matrix:set mS y x [S] of patch x y
      set y y + 1]
    set y 0
  set x x + 1]
end

to set-mA ; create a matrix mA with A of patches
  set mA matrix:make-constant (max-pycor + 1) (max-pxcor + 1) 0
  let x 0
  let y 0
  while [x <= max-pxcor] [
    while [y <= max-pycor] [
      matrix:set mA y x [A] of patch x y
      set y y + 1]
    set y 0
  set x x + 1]
end

to set-mDs ; create a matrix mDs with Ds of patches
  set mDs matrix:make-constant (max-pycor + 1) (max-pxcor + 1) 0
  let x 0
  let y 0
  while [x <= max-pxcor] [
    while [y <= max-pycor] [
      matrix:set mDs y x [Ds] of patch x y
      set y y + 1]
    set y 0
  set x x + 1]
end

to set-mDa ; create a matrix mDa with Da of patches
  set mDa matrix:make-constant (max-pycor + 1) (max-pxcor + 1) 0
  let x 0
  let y 0
  while [x <= max-pxcor] [
    while [y <= max-pycor] [
      matrix:set mDa y x [Da] of patch x y
      set y y + 1]
    set y 0
  set x x + 1]
end

to set-mM ;create a matrix mM with the mass of susceptible cells of each grid cell/patch
  set mM matrix:make-constant (max-pycor + 1) (max-pxcor + 1) 0
  let x 0
  let y 0
  while [x <= max-pxcor] [
    while [y <= max-pycor] [
      matrix:set mM y x [sum [m] of bacteria-here] of patch x y
      set y y + 1]
    set y 0
  set x x + 1]
end

to set-mMD ;create a matrix mMD with the mass of dead cells of each grid cell/patch
  set mMD matrix:make-constant (max-pycor + 1) (max-pxcor + 1) 0
  let x 0
  let y 0
  while [x <= max-pxcor] [
    while [y <= max-pycor] [
      matrix:set mMD y x [sum [m] of dead-here] of patch x y
      set y y + 1]
    set y 0
  set x x + 1]
end
@#$#@#$#@
GRAPHICS-WINDOW
957
40
1682
766
-1
-1
27.6
1
10
1
1
1
0
1
0
1
0
25
0
25
1
1
1
min
30.0

BUTTON
6
10
69
43
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
73
10
136
43
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

SLIDER
166
220
321
253
delta-l
delta-l
0
10
4.0
1
1
µm
HORIZONTAL

SLIDER
166
255
321
288
delta-t
delta-t
0.0001
0.01
0.0035
0.0001
1
s
HORIZONTAL

SLIDER
164
88
321
121
D-substrate
D-substrate
1
1000
900.0
1
1
µm²/s
HORIZONTAL

MONITOR
7
198
73
243
dt.D/(dl²)
delta-t * D-substrate / (delta-l * delta-l)
3
1
11

SLIDER
164
53
321
86
Cs-bulk
Cs-bulk
0
0.5
0.4
0.01
1
g/L
HORIZONTAL

TEXTBOX
79
207
132
233
Have to be < 1/4 !
11
0.0
1

TEXTBOX
165
35
315
53
Substrate
12
0.0
1

SLIDER
523
88
665
121
D-antibiotic
D-antibiotic
0
1000
900.0
1
1
µm²/s
HORIZONTAL

SLIDER
340
51
489
84
initial-bacteria
initial-bacteria
0
200
27.0
1
1
NIL
HORIZONTAL

SLIDER
340
86
489
119
µmax
µmax
0
1.5
1.25
0.01
1
h-1
HORIZONTAL

SLIDER
339
122
489
155
Ks
Ks
0
0.01
0.0035
0.0001
1
g/L
HORIZONTAL

TEXTBOX
341
35
491
53
Cells
12
0.0
1

SLIDER
341
293
479
326
delta-tb
delta-tb
0
120
67.0
1
1
s
HORIZONTAL

SLIDER
340
230
489
263
m-div
m-div
0
1000
493.0
1
1
fg
HORIZONTAL

MONITOR
6
55
105
100
World width (µm)
(max-pxcor + 1) * delta-l
0
1
11

PLOT
13
658
468
812
Cells
Ticks
Cells
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Susceptible cells" 1.0 0 -13840069 true "" "plot count bacteria"
"Dead cells" 1.0 0 -11053225 true "" "plot count dead"

BUTTON
690
50
792
83
Color lineages
color-lineage
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

TEXTBOX
170
299
320
317
Shoving
12
0.0
1

BUTTON
690
86
793
119
Color in green
ask bacteria [set color green]
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

MONITOR
836
324
923
369
Time (h)
t / 3600
2
1
11

SWITCH
690
124
793
157
color-µ?
color-µ?
0
1
-1000

TEXTBOX
799
126
908
164
Green - cell growth rate
11
0.0
1

MONITOR
6
102
105
147
World heigh (µm)
delta-l * (max-pycor + 1)
0
1
11

MONITOR
6
150
113
195
Number of grid cells
count patches
0
1
11

SLIDER
165
185
321
218
boundary-layer-lenght
boundary-layer-lenght
0
60
40.0
1
1
µm
HORIZONTAL

TEXTBOX
691
33
780
51
Coloring
12
0.0
1

TEXTBOX
342
276
514
294
Frequency of cells' update
12
0.0
1

TEXTBOX
405
380
553
398
Phenotypic switches
12
0.0
1

SLIDER
404
433
546
466
a-max
a-max
0
1
1.0
0.01
1
h-1
HORIZONTAL

SLIDER
404
506
545
539
K
K
0
0.01
0.0035
0.0001
1
g/L
HORIZONTAL

SLIDER
404
542
545
575
K'
K'
0
0.00001
5.0E-6
0.000001
1
g/L
HORIZONTAL

SLIDER
404
469
546
502
b-max
b-max
0
1
0.11
0.01
1
h-1
HORIZONTAL

MONITOR
837
420
923
465
Persisters
count persisters
0
1
11

MONITOR
837
371
946
416
Susceptible cells
count bacteria
17
1
11

SLIDER
340
158
490
191
Yxs
Yxs
0
1
0.21
0.01
1
NIL
HORIZONTAL

SLIDER
522
124
665
157
MIC
MIC
0
0.01
5.0E-6
0.0001
1
g/L
HORIZONTAL

TEXTBOX
525
34
602
52
Antibiotic
12
0.0
1

SLIDER
522
52
665
85
Ca-bulk
Ca-bulk
0
1500
1000.0
1
1
*MIC
HORIZONTAL

SLIDER
523
183
665
216
k-max
k-max
0
15
10.0
0.01
1
h-1
HORIZONTAL

SLIDER
523
254
665
287
Kk
Kk
0
10
10.0
0.1
1
*MIC
HORIZONTAL

SWITCH
167
501
300
534
Start-treatment?
Start-treatment?
1
1
-1000

PLOT
479
848
865
1000
Mean solute concentrations out bulk
Ticks
g/L
0.0
10.0
0.0
0.01
true
true
"" ""
PENS
"Substrate" 1.0 0 -13345367 true "" "plot sum [S] of patches with [bulk = 0] / count patches with [bulk = 0]"
"Antibiotic" 1.0 0 -5825686 true "" "plot sum [A] of patches with [bulk = 0] / count patches with [bulk = 0]"

MONITOR
838
468
923
513
Dead cells
count dead
17
1
11

SWITCH
690
161
793
194
show-s?
show-s?
0
1
-1000

SWITCH
690
199
793
232
show-ATB?
show-ATB?
1
1
-1000

PLOT
22
1011
414
1146
% persisters
Ticks
%
0.0
1.0
0.0
0.1
true
true
"" ""
PENS
"%pers" 1.0 0 -16777216 true "" "plot count persisters / (count bacteria + count persisters) * 100"

TEXTBOX
799
164
920
195
Blue - substrate concentration
11
0.0
1

TEXTBOX
799
202
914
233
Orange - antibiotic concentration
11
0.0
1

SLIDER
340
195
490
228
cell-density
cell-density
20
400
200.0
1
1
fgS/µm3
HORIZONTAL

SLIDER
167
394
323
427
auto-start-treatment
auto-start-treatment
0
15
5.0
1
1
h
HORIZONTAL

SLIDER
523
219
665
252
k-max-pers
k-max-pers
0
1
0.1
0.01
1
h-1
HORIZONTAL

TEXTBOX
306
510
378
528
(Manual start)
11
0.0
1

TEXTBOX
168
377
318
395
Treatment
12
0.0
1

SLIDER
167
430
323
463
treatment-duration
treatment-duration
0
10
1.0
1
1
h
HORIZONTAL

SLIDER
164
149
321
182
D*
D*
0
1
0.8
0.1
1
NIL
HORIZONTAL

SLIDER
168
317
306
350
max-moving-cells
max-moving-cells
0
10
10.0
1
1
%
HORIZONTAL

PLOT
18
850
403
1004
Biofilm thickness and roughness
Ticks
µm
0.0
1.0
0.0
1.0
true
true
"" ""
PENS
"Mean thickness" 1.0 0 -8630108 true "" "plot mean-thickness"
"Roughness coefficient" 1.0 0 -955883 true "" "plot Ra"

SWITCH
8
262
111
295
Video?
Video?
1
1
-1000

BUTTON
690
252
889
285
Restore default parameters
restore-default-values
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

PLOT
479
660
892
810
Persisters
Ticks
Cells
0.0
10.0
0.0
10.0
true
true
"" ""
PENS
"Persisters" 1.0 0 -2674135 true "" "plot count persisters"

TEXTBOX
18
825
168
860
Biofilm structure
12
0.0
1

TEXTBOX
15
636
165
654
Cell population dynamics
12
0.0
1

SLIDER
167
466
322
499
recovery-time
recovery-time
0
10
1.0
1
1
h
HORIZONTAL

SLIDER
404
397
546
430
Strategy
Strategy
1
4
4.0
1
1
NIL
HORIZONTAL

SLIDER
523
313
665
346
Lysis
Lysis
0
0.1
0.085
0.001
1
h-1
HORIZONTAL

TEXTBOX
553
389
776
494
Strategy:\n1) Stochastic switches\n2) Substrate-dependent switches\n3) Antibiotic-dependent switches\n4) Antibiotic and substrate switches - used for our study
12
0.0
1

TEXTBOX
481
828
631
846
Solutes
12
0.0
1

TEXTBOX
9
302
159
358
If \"Video?\" is on, a view of the computational domain will be exported as a .png file every cell update
11
0.0
1

TEXTBOX
1009
10
1241
44
------- Computational domain -------
14
0.0
1

TEXTBOX
798
97
933
115
Color all cells in green
11
0.0
1

TEXTBOX
799
54
953
82
Attribute one color to the lineage of each initial cell
11
0.0
1

TEXTBOX
525
294
675
312
Lysis rate of dead cells
11
0.0
1

TEXTBOX
525
165
633
183
Death rates
11
0.0
1

TEXTBOX
390
604
567
624
--------- Graphs ---------
14
0.0
1

TEXTBOX
410
10
606
32
--------- Parameters ---------
14
0.0
1

PLOT
508
1024
708
1174
plot 1
ticks
tottime
0.0
10.0
0.0
10.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "plot tottime"

@#$#@#$#@
## WHAT IS IT?

This is an individual-based biofilm model. Cells can switch between two phenotypes: susceptible and persister. Persisters do not grow but are tolerant to the antibiotic. Different switching strategies between these two phenotypes are implemented. Antibiotic treatments can be simulated to study the survival and recovery of the biofilm.

For more information:
Carvalho et al. (2018), How do environment-dependent switching rates between susceptible and persister cells affect the dynamics of biofilms faced with antibiotics?, npj Biofilms and Microbiomes 4:6, doi:10.1038/s41522-018-0049-2

## HOW IT WORKS

Cells are initialized on a surface and grow, divide and shove each other to generate a biofilm. Dissolved substrate and antibiotic diffuse and react with the cells in a diffusion-reaction fashion. They diffuse from a bulk liquid above the biofilm, where concentrations are maintained constant, across a boundary layer toward the biofilm. In the model, bacteria can be susceptible, persister or dead.

## HOW TO USE IT

Setup and go!

## THINGS TO NOTICE

The dynamics of survival and recovery of the biofilm vary greatly depending on the switching strategy selected and the parameters a and b.

## THINGS TO TRY

(suggested things for the user to try to do (move sliders, switches, etc.) with the model)

## EXTENDING THE MODEL

(suggested things to add or change in the Code tab to make the model more complicated, detailed, accurate, etc.)

## NETLOGO FEATURES

(interesting or unusual features of NetLogo that the model uses, particularly in the Code tab; or where workarounds were needed for missing features)

## RELATED MODELS

(models in the NetLogo Models Library and elsewhere which are of related interest)

## CREDITS AND REFERENCES

Carvalho et al. (2018), How do environment-dependent switching rates between susceptible and persister cells affect the dynamics of biofilms faced with antibiotics?, npj Biofilms and Microbiomes 4:6, doi:10.1038/s41522-018-0049-2
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

sheep
false
15
Circle -1 true true 203 65 88
Circle -1 true true 70 65 162
Circle -1 true true 150 105 120
Polygon -7500403 true false 218 120 240 165 255 165 278 120
Circle -7500403 true false 214 72 67
Rectangle -1 true true 164 223 179 298
Polygon -1 true true 45 285 30 285 30 240 15 195 45 210
Circle -1 true true 3 83 150
Rectangle -1 true true 65 221 80 296
Polygon -1 true true 195 285 210 285 210 240 240 210 195 210
Polygon -7500403 true false 276 85 285 105 302 99 294 83
Polygon -7500403 true false 219 85 210 105 193 99 201 83

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

wolf
false
0
Polygon -16777216 true false 253 133 245 131 245 133
Polygon -7500403 true true 2 194 13 197 30 191 38 193 38 205 20 226 20 257 27 265 38 266 40 260 31 253 31 230 60 206 68 198 75 209 66 228 65 243 82 261 84 268 100 267 103 261 77 239 79 231 100 207 98 196 119 201 143 202 160 195 166 210 172 213 173 238 167 251 160 248 154 265 169 264 178 247 186 240 198 260 200 271 217 271 219 262 207 258 195 230 192 198 210 184 227 164 242 144 259 145 284 151 277 141 293 140 299 134 297 127 273 119 270 105
Polygon -7500403 true true -1 195 14 180 36 166 40 153 53 140 82 131 134 133 159 126 188 115 227 108 236 102 238 98 268 86 269 92 281 87 269 103 269 113

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.3.0
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
<experiments>
  <experiment name="biofilm search" repetitions="4" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>count turtles</metric>
    <metric>count bacteria</metric>
    <metric>count persisters</metric>
    <metric>count dead</metric>
    <metric>avg-thickness</metric>
    <metric>Ra</metric>
    <metric>up-persisters</metric>
    <metric>down-persisters</metric>
    <metric>Ov</metric>
    <enumeratedValueSet variable="Kk">
      <value value="6.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="3.5E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Profiler?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="3.5E-4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-s">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy2?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="2"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy3?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="S-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-c">
      <value value="1"/>
      <value value="0.1"/>
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-c">
      <value value="1"/>
      <value value="0.1"/>
      <value value="0.01"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-atb">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy4?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="ATB-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-atb">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="half-life">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-s">
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-atb">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-m?">
      <value value="false"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment6" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>a-max</metric>
    <metric>b-max</metric>
    <metric>recovery-time</metric>
    <metric>treatment-duration</metric>
    <metric>count persisters</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="0.167"/>
      <value value="0.33"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cs-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ca-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="493"/>
    </enumeratedValueSet>
    <steppedValueSet variable="b-max" first="0.1" step="0.3" last="1"/>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="67"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kk">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-antibiotic">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="27"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="0.167"/>
      <value value="0.33"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <steppedValueSet variable="a-max" first="0.1" step="0.3" last="1"/>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="smallervalues" repetitions="1" sequentialRunOrder="false" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>a-max</metric>
    <metric>b-max</metric>
    <metric>recovery-time</metric>
    <metric>treatment-duration</metric>
    <metric>count persisters</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="0.005"/>
      <value value="0.01"/>
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cs-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ca-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-max">
      <value value="0.1"/>
      <value value="0.4"/>
      <value value="0.7"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kk">
      <value value="6.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.05"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-antibiotic">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="27"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="0.005"/>
      <value value="0.01"/>
      <value value="0.03"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-max">
      <value value="0.1"/>
      <value value="0.4"/>
      <value value="0.7"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment5" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>a-max</metric>
    <metric>b-max</metric>
    <metric>recovery-time</metric>
    <metric>treatment-duration</metric>
    <metric>count persisters</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cs-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ca-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-max">
      <value value="0.1"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="67"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kk">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="1"/>
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-antibiotic">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="27"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-max">
      <value value="0.1"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experiment8" repetitions="1" runMetricsEveryStep="true">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>a-max</metric>
    <metric>b-max</metric>
    <metric>recovery-time</metric>
    <metric>treatment-duration</metric>
    <metric>count persisters</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment1">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="0.167"/>
      <value value="0.33"/>
      <value value="0.8"/>
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cs-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ca-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-max">
      <value value="0.2"/>
      <value value="0.55"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="67"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kk">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-antibiotic">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="27"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="0.167"/>
      <value value="0.33"/>
      <value value="0.8"/>
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-max">
      <value value="0.2"/>
      <value value="0.55"/>
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="experimentab" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>a-max</metric>
    <metric>b-max</metric>
    <metric>recovery-time</metric>
    <metric>treatment-duration</metric>
    <metric>count persisters</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="0.6"/>
      <value value="5"/>
      <value value="30"/>
      <value value="60"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cs-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ca-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-max">
      <value value="0.25"/>
      <value value="0.55"/>
      <value value="0.85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="67"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kk">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-antibiotic">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="27"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="0.3"/>
      <value value="10"/>
      <value value="30"/>
      <value value="48"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-max">
      <value value="0.25"/>
      <value value="0.55"/>
      <value value="0.85"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="largeexperiment1" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>a-max</metric>
    <metric>b-max</metric>
    <metric>recovery-time</metric>
    <metric>treatment-duration</metric>
    <metric>count persisters</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="0.33"/>
      <value value="0.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cs-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ca-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-max">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.4"/>
      <value value="0.55"/>
      <value value="0.7"/>
      <value value="0.85"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="67"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kk">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-antibiotic">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="27"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="0.166"/>
      <value value="0.33"/>
      <value value="0.5"/>
      <value value="0.66"/>
      <value value="0.833"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-max">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.4"/>
      <value value="0.55"/>
      <value value="0.7"/>
      <value value="0.85"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="largeexperiment2" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>a-max</metric>
    <metric>b-max</metric>
    <metric>recovery-time</metric>
    <metric>treatment-duration</metric>
    <metric>count persisters</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="0.66"/>
      <value value="0.833"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cs-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ca-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-max">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.4"/>
      <value value="0.55"/>
      <value value="0.7"/>
      <value value="0.85"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="67"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kk">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-antibiotic">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="27"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="0.166"/>
      <value value="0.33"/>
      <value value="0.5"/>
      <value value="0.66"/>
      <value value="0.833"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-max">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.4"/>
      <value value="0.55"/>
      <value value="0.7"/>
      <value value="0.85"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="largeexperiment3" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>a-max</metric>
    <metric>b-max</metric>
    <metric>recovery-time</metric>
    <metric>treatment-duration</metric>
    <metric>count persisters</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="1"/>
      <value value="1.16"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cs-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ca-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-max">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.4"/>
      <value value="0.55"/>
      <value value="0.7"/>
      <value value="0.85"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="67"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kk">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-antibiotic">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="27"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="0.166"/>
      <value value="0.33"/>
      <value value="0.5"/>
      <value value="0.66"/>
      <value value="0.833"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-max">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.4"/>
      <value value="0.55"/>
      <value value="0.7"/>
      <value value="0.85"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="largeexperiment4" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>a-max</metric>
    <metric>b-max</metric>
    <metric>recovery-time</metric>
    <metric>treatment-duration</metric>
    <metric>count persisters</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="1.33"/>
      <value value="1.5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cs-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ca-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="500"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-max">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.4"/>
      <value value="0.55"/>
      <value value="0.7"/>
      <value value="0.85"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="67"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kk">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-antibiotic">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="27"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="0.166"/>
      <value value="0.33"/>
      <value value="0.5"/>
      <value value="0.66"/>
      <value value="0.833"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-max">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.4"/>
      <value value="0.55"/>
      <value value="0.7"/>
      <value value="0.85"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
  <experiment name="largeexperiment" repetitions="1" runMetricsEveryStep="false">
    <setup>setup</setup>
    <go>go</go>
    <metric>ticks</metric>
    <metric>a-max</metric>
    <metric>b-max</metric>
    <metric>recovery-time</metric>
    <metric>treatment-duration</metric>
    <metric>count persisters</metric>
    <metric>count bacteria</metric>
    <enumeratedValueSet variable="Ks">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-t">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="cell-density">
      <value value="200"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="recovery-time">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
      <value value="40"/>
      <value value="50"/>
      <value value="60"/>
      <value value="70"/>
      <value value="80"/>
      <value value="90"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K">
      <value value="0.0035"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Cs-bulk">
      <value value="0.4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Ca-bulk">
      <value value="1000"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="m-div">
      <value value="493"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="b-max">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.4"/>
      <value value="0.55"/>
      <value value="0.7"/>
      <value value="0.85"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="max-moving-cells">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-tb">
      <value value="67"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Kk">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="color-µ?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D*">
      <value value="0.8"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="delta-l">
      <value value="4"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Lysis">
      <value value="0.085"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="boundary-layer-lenght">
      <value value="40"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="auto-start-treatment">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Yxs">
      <value value="0.21"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-ATB?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="MIC">
      <value value="5.0E-6"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max">
      <value value="10"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="µmax">
      <value value="1.25"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-antibiotic">
      <value value="900"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="show-s?">
      <value value="true"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="k-max-pers">
      <value value="0.1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="initial-bacteria">
      <value value="27"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="K'">
      <value value="1.0E-5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="treatment-duration">
      <value value="5"/>
      <value value="10"/>
      <value value="20"/>
      <value value="30"/>
      <value value="40"/>
      <value value="50"/>
      <value value="60"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="a-max">
      <value value="0.1"/>
      <value value="0.25"/>
      <value value="0.4"/>
      <value value="0.55"/>
      <value value="0.7"/>
      <value value="0.85"/>
      <value value="1"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Start-treatment?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Strategy">
      <value value="5"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="Video?">
      <value value="false"/>
    </enumeratedValueSet>
    <enumeratedValueSet variable="D-substrate">
      <value value="900"/>
    </enumeratedValueSet>
  </experiment>
</experiments>
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
