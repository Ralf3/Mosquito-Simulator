# Mosquito-Simulator
Simulates the spread of WestNilVirus a combination of ODE and spatial movement (Julia implementation)

## General Remarks
The Mosquito Simulator is completely written in Julia. The version presented here 
contains executable code, with the restriction that of the available regions only the region Magdeburg is available. For this region the weather data is also provided. If you want to test other regions, please send me a personal message.
### Installation
The simulator has been tested under Julia version 1.8.3, but should also run with subsequent versions. To install the software please follow the steps below:

- Install Julia on your computer (I use Linux as OS).
- simu.install contains all needed packages for the simulator. It can directly install the packages install.

## Simulation control
The file "control.jl" contains the central parameters for controlling a simulation: <br>
to start the simulation, please include inside of Julia:
**include("control.jl")**

### General simulation parameters

|Parameter|Explanation|
| :------------ |:---------------:|
| R="Dueben"  |  Selection of the region: Bitterfeld, Dueben, Radebeul, Magdeburg, Trebnitz|
| _from_=20190101| loads the wheater data (DWD) for the years: _from _to|
| _to_ = 20201231| use the format year month day|
| species="pipiens"| pipiens is the only one in this release |
| birds="sperling" | birds (German): sperling,gruenfink,amsel,elster,habicht,kraehe|
| |		   please use German for sparrow,greenfinch,blackbird,magpie,hawk,crow|
| model=false | if the habitat model is used: true|
| remark="None" | a short (12 letter) comment |

### Dynamic parameters of the simulation

|Parameter|Explanation|
| :------------ |:---------------:|
| u0b=[500.0,0.0,20.0,0.0,0.0]| initial values of bird(t=0) (SEIRD)|
| u0m=[200.0,200.0,10.0,10.0]| initial values of mosquitoes (LSEI)|
| steps=200| number of simulation steps: steps |
| infstepIB=220| step when infected (E) birds come in (start:start+steps)|
| addIB=10.0| number of infected (E) birds Float64|
| infstepIM=200 | step when infectious (I) mosquitoes come in (start:start+steps|
| addIM=20.0| number of infectious (I) mosquitoes Float64|
| walkn=1| number of walks/day should be 1|
|initial_targets="model"| where the simulation starts from ["siedlung","model","river"]|

### Description of the simulation

The step size (dt) of the simulation corresponds to one day. In the function "simulation3"
is running the MainLoop of the simulation. It controls the movement of the mosquitoes by
distributing the mosquitoes in the region from the solution of the ODE distributing the
mosquitoes in the region. The distribution of the birds follows a similar pattern. The
treatment of the mosquitoes is implemented in the module mosquito_struct.jl, that of
the birds in bird_struct.jl.

#### Selected commands

**add_m**(n::Int64,DoY::Int64,gx::Grid,state::Int64,limit=30.0f0) <br>
add (n>0) or remove mosquitoes (n<0) <br> 
n: number of mosquitoes; gx: Grid: initial_target: ["siedlung","river","model"] <br>
state: 2,3,4 (SEI) <br>
limit: is the maximum distance to a perfect habitat <br>

**add_b**(n::Int64,DoY::Int64,gx::Grid,state::Int64,limit=3.0f0) <br>
similar procedure for birds: <br>
if the birds come in the region a place with a distance<limit will be selected <br>
if the bird is already in: <br>
RB ==> converts an infectious bird into a recovered bird <br>
DB ==> converts an infectious bird in a dead bird <br>

The movement of the mosquitoes is carried out on a square grid of 25km*25km with a
cell size of 100*100m. The mosquito can move over a number of time steps (walkn)
move freely in the region. <br>

All mosquitoes are stored in the global array: array_of_mosqu=[] (see util.jl)
all birds are stored in the global array: array_of_birds=[] (see util.jl)
If a mosquito dies or leaves the model region its state is set zombie,
which is numeric as 1Z with Z=[SM=2,EM=3,IM=4] hence a 12 encodes a "zombie SM".
Since birds stay near their nesting place, zombies play only a minor role. <br>

## Evaluation of the dynamics of a simulation

An overview of the development of mosquitoes and birds is shown in the following functions:

**show_resm()**  shows the simulation development (Mosquito) LS (top) and EI (bottom) over time <br>
**show_resb()**  shows the simulation development (Bird) SI (top) and RD (bottom) over time. <br>

### Some control functions

The **solj2df()** function calculates a DataFrame that displays the development
(L,SM,EM,IM,SB,EB,IB,RB,DB) as a table. dfx=solj2df(); is called at the end of control.jl.

#### df1 displays lines from 1:steps as table.

| Row | L | SM | EM | IM | SB | EB | IB | RB | DB |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
|   1 | 100.0  | 500.0 | 0.0   | 50.0  | 500.0 | 30.0  | 30.0  | 0.0   | 0.0   |
|   2 | 129.84 | 488.0 | 1.0   | 49.0  | 516.0 | 11.0  | 42.0  | 2.0   | 10.0  |
|   3 | 125.38 | 475.0 | 1.0   | 45.0  | 535.0 | 7.0   | 29.0  | 4.0   | 27.0  |
|   4 | 90.979 | 453.0 | 1.0   | 41.0  | 554.0 | 3.0   | 25.0  | 5.0   | 38.0  |

#### selected visualization commands of dfx

**@df dfx plot([:L,:SM])** allows to plot the history of selected columns ([:L,:SM]) of dfx <br>
**@df dfx boxplot(label=["IM" "IB"],[:IM,:IB])** visualizes IM and IB as boxplot <br>
**@df dfx histogram([:IB],bins=10)**  displays a histogram of aselected column (IB) <br> 
**@df dfx density([:IB,:IM],label=["IB" "IM"])** visualizes the density of columns (IB,IM) <br>
**@df dfx marginalkde(:SB,:SM)** visualizes two-dimensional density splats (SB,SB) <br>
**@df dfx corrplot(cols(1:4), grid = true)** vsualizes a multidimensional correlation plot <br>
**@df dfx corrplot(cols([:SB,:IB,:RB,:DB]), grid = true)** further correlation plot <br>

More options can be found at https://github.com/JuliaPlots/StatsPlots.jl <br>

## Evaluation of the spatial distribution of mosquitoes and birds

### Visualization of the distribution of mosquitoes and birds

#### Visualization of the distribution of mosquitoes

The mosquitoes allow a complex pattern of movement: they can fly in an undirected manner,
they can fly directionally e.g. to search for food or blood feeding, they can leave or enter
the region. At contrast, birds can only enter the area, but then stay close to their nesting site.
We start with the most important features for spatial visualization of mosquitoes: <br>

**show_map_m(sel::Int64**) sel={SM=2,EM=3,IM=4} <br>
displays a map with the mosquitoes on the desired grid <br>
if you want to display corresponding zombies, you can use the function <br>
**gen_hist_m()** ==> Dict  about all possible stages of the mosquitoes <br>
**show_map_m(sel::Int64)** sel={2,12,3,13,4,14} corresponding to the stages including the zombies
<br>

This basic function can be adapted to the desired visualization by a number of modifications:<br>

**switch(sel::String)** sel=["settlement", "model", "river"] allows to change the background map <br>
**switch(t::Int64)**  sel {1:steps} allows the selection of a time step e.g. sel=125 <br>

#### Visualization of the distribution of birds

The birds are considered stationary, so a switch(sel::String) is ignored, i.e.
the birds always have the settlement as background map. <br>

**show_map_b(sel::Int64)** sel in {SB=1,EB=2,IB=3,RB=4;DB=5} <br>
for zombies use sel in {11,12,13,14,15} <br>
**gen_hist_b() ==> Dict** over all possible stages of the birds <br>
**switch(t::Int64)** sel {1:steps} is also effective for the birds

Remark: storing each time step requires a lot of memory and slows down the simulation. If you can
live without the switch at any point of the simulation, please comment out the lines
ABM[k]=deepcopy(array_of_mosqu) and BRD[k]=deepcopy(array_of_birds). Then a call to
switch(sel::Int64), but also (gen_hist_m and gen_hist_b) would lead to an error message.
The simulator could then run even on weak PCs.

## Store of the simulations (controls) in a database

**df,dfx=save_control()**; creates a new entry at the end of a simulation
into the database "simulations.csv", stores the control_{ID}.jl under controls.jl and
the result of the simulation under results/result_{ID}.csv. The df and the dfx are usable
afterwards in the simulation environment.

### DataFrame can be adjusted and saved:

**switch_valid(ID::Int64)** inverts the valid flag of the line with the ID to false or true <br>
**set_remark(ID::Int64,s::String)** overwrites the remark with s <br>
**sort_df()** sorts the table order=[:R, :species, :start]) <br>
**save_df(df::DataFrame)**  overwrite the old simulation.csv with the modified df<br>
**clean_controls()** removes data with valid=false from both df and <br>
from the controls. These cannot be restored! <br>
**activate(ID::Int64)**  overwrites the current control.jl with the control_{ID}.jl<br>

Remark: Working with the database is often limited to setting the valid flag to false and
the subsequent removal of the selected simulation from the database by clean_controls.

