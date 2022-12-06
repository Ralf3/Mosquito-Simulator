sensitivity.jl is used to analyze the parameters for the DGL
============================================================

sensitivity.jl loads all important modules used by the simulator. This enables a comprehensive, 
realistic analysis. However, the focus is on the analysis of the dynamic behavior of the simulation, 
which is realized by the module DGL.jl in the simulator. DGL.jl is emulated in sensitivity.jl. 
Data from other modules, which are needed for the simulation, e.g. the DWD data are loaded by 
DWD.jl exactly as in the simulator. 

# include all general parameters
# ==============================

The general parameters are used to define the simulation as a whole. This includes e.g. the 
initial conditions, the mosquito species, the bird species, the simulation period, etc.
Please modify this according your usecase.

u0m=[100.0,500.0,0.0,50.0]    # start values for mosquitoes (LSEI)
u0b=[500.0,30.0,30.0,0.0,0.0] # start values for birds (SEIRD)
u0bx=copy(u0b)                # copy start values for birds
u0mx=copy(u0m)                # copy start values for mosquitoes
R = "Dueben"                  # Bitterfeld, Dueben, Radebeul, Magdeburg, Trebnitz
_from_=20190101               # year month day for DWD
_to_ = 20201231               # year month day for DWD
start=20200430                # start of simulation
species="pipiens"             # define the species (pipiens, vexans)
birds="elster"                # bird: sperling,gruenfink,amsel,elster,habicht,kraehe
                              # sparrow,greenfinch,blackbird,magpie,hawk,crow
			      # please use german names 
master="Jump"                 # the branch of the git repository Jump is the only valid one.
remark="testX"                # 12 letters limit
steps=200                     # number of simulation steps
infstepIB=0.                  # set infstepIB day or 0. ==> no step
addIB= 20                     # default only 1 step
infstepIM=0.                  # set infstepIM day or 0. ==> no step
addIM = 10                    # default only 1 step
walkn=1                       # number of walks/day

# Parameter dictionary for analysis
# =================================

The dictionary contains all parameters which can be chnages in the simulation. 
This includes parameters and initial conditions. Don't change this!

param=OrderedDict()
param["αB"]=αB          # recovery rate of bird
param["λ"]=λ            # additional factor for small birds  exposed_rateB
param["γB"]=γB          # incubition rate of bird
param["νB"]=νB          # is killed by the WNV
param["mB"]=mB          # mortality of a bird
param["KM"]=KM          # carrying capacity mosquito
param["KB"]=KB          # carrying capacity bird
param["pM"]=pM          # prob. transmission WNV from mos to the bird
param["pB"]=pB          # prop. transmission WNV from bird to mos ∈ [0.125,0.2]
param["ϕB"]=ϕB          # Mosquito-to-bird-ratio
param["u1"]=u0[1]       # initial values for mos and birds
param["u2"]=u0[2]	
param["u3"]=u0[3]
param["u4"]=u0[4]
param["u5"]=u0[5]
param["u6"]=u0[6]
param["u7"]=u0[7]
param["u8"]=u0[8]
param["u9"]=u0[9]

# Modify and start the simulation
# ===============================

Before you can start a valid simulation with sensitivity of parameters please call
reset_simu() # to make the simulation clean.

The simulation itself is started using 
simulate(n,dx)         # The number of simulation runs n=50 and dx=0.1 should be adapted.
                       # n=10 is a good starting
to check the siulation and the ranges of the result. n>100 is used to
get reliable results but it takes longer.

In the following some parameter are used for sensitivity analysis

1.) linear modifications over simulation time
pM=param["pM"]*(1+0.01*i/n)     # the prob. transmission WNV from mosquito to the
                                # bird is linear chaned with time
KM=param["KM"]*(1-1/n+i/n)      # the carrying capacity mosquito is linear chnaged with time 

2.) Stochastic modification of a parameter
KM=param["KM"]*(1.0-dx*randn()) # stochastic KM

Please adapt one or more parameters accrording your needs for analysis

# Analysis of two simulation results
# ==================================

For simple visualization the are three functions:
show_resm()      # shows the simulation results without any change (blue) versus the 
                 # last simulation result with changes (red) for the mosquites (LSEI)
show_resb()      # shows the simulation results without any change (blue) versus the
                 # last simulation result with changes (red) for the birds(BIRD)
compare_sols(nr) # shows a selected nr of the LM=1 SM=2, ... see param
                 # (blue=reset; red=simulation)

# Statistical analysis of the simulation 
# ======================================

For each simulation run (n) the means of parameters are stored in the DataFrame df1

 Row │ Lm       SMm      EMm      IMm      SBm      EBm      IBm      RBm      DBm     
     │ Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64  Float64 
─────┼─────────────────────────────────────────────────────────────────────────────────
   1 │ 1475.08  2365.65  3.99005  15.7065  648.109  1.42786  2.98507  99.0945  111.532
   2 │ 1391.93  2210.73  3.77612  15.1095  637.08   1.36816  2.75622  95.7512  109.214
   3 │ 1380.83  2157.83  1.63184  12.1791  627.607  1.13433  2.88557  80.6617  105.801
   4 │ 1406.09  2202.74  2.80597  15.0547  617.632  1.46766  2.80597  94.204   112.303
   5 │ 1391.45  2228.92  2.34826  15.7114  598.751  1.52239  2.8408   93.4876  111.965
   6 │ 1252.92  1952.81  1.86567  14.3781  599.114  1.44279  2.81592  89.0697  110.403
...

The following visualizations are example for analysis of the df1

@df df1 boxplot(title="RBm versus DBm",label=["RBm" "DBm"],[:RBm,:DBm]) 
# it shows the boxplot comparing the recovered birds (RBm) with the dead birds (DBm)

@df df1 boxplot(label=["IMm" "IBm"],[:IMm,:IBm])
# It shows the boxplot comparing the infectious mosquitoes (IMm) with the infectious birds (IBm)

@df df1 histogram([:SBm],bins=10)
# it shows a histogram of the birds

@df df1 density([:IBm,:IMm])
# it shows a density plot of IBm (blue) versus IMm (red)

@df df1 marginalkde(:SBm,:SMm)
# it shows a marginalkde between SBm and SMm

@df df1 marginalscatter(:SBm,:SMm)
# it shows a marginal scatter between SBm and SMm

@df df1 marginalhist(:SBm,:SMm)
# ist shows the marginalhist between SBm and SMn

Other visualizations based on a DataFrame can be found at:
https://github.com/JuliaPlots/StatsPlots.jl
