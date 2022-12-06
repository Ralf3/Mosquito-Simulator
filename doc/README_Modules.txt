Individual modules used for simulation
======================================

bird_parameter.jl
-----------------
contains all paramters for birds like:
kreahe_bB() ==> bB  for the Crow, Hawk, Blackbird, Magpie, Greenfinch, Sparrow

image.jl
--------
handling of OSM images: help function not used after preprocessing

database.jl
-----------
store and handle df in simulations.csv and control_ID.jl under controls

data.jl
-------
loads the spatial data 800d and scales it to 250*250
loads the DWD-data


grid.jl
-------
one of the main modules based on Samt2: 
handels ASCII grids (read/write) many data transformation from numpy array as npz; npz_Grid
grid operations:    norm_Grid!(G::Grid); znorm_Grid!(G::Grid); invert_Grid!(G::Grid);
                    scale_Grid(G::Grid,nrows::Int64,ncols::Int64); enlarge_Grid(G::Grid,factor::Int64);
                    random_Points_Grid(N::Int64,G::Grid);random_Points_Grid(N::Int64,sel::Float32,G::Grid;
                    add_Points_Grid(points::Vector{Point},G::Grid,value=1.0);
                    set_Points_Grid(points::Vector{Point},G::Grid,value=1.0);
                    sample_Point_Grid(P::Vector{Point},G::Grid);
                    distance_Points_Grid(P::Vector{Point},G::Grid);distance_Grid_Sel(G::Grid,sel::Float32);
                    diffx_Grid(G::Grid); diffy_Grid(G::Grid);
                    
util.jl
-------
loads the Julia packages using Pkg
implements plot functions:
    show_resm(); show_resb();
    solj2df() ==> DataFrame;
    switch(sel::String)
loads all modules for simulation: 
    include("DWD.jl")               # base data to handle the DWD
    include("grid.jl")              # base data to handle the OSM
    include("data.jl")              # loads the spatial OSM data and the LT from DWD
    include("DGL.jl")               # load the paramter of the DGL models SIRM SIRB
    include("bird_struct.jl")       # load the bird ABM
    include("mosquito_struct.jl")   # load the mosquito ABM
    include("database.jl")          # store the results in simulation.csv

mosquito_struct.jl
------------------
implements a mutable struct Mosquito including some functions about it:
    Mosquito(p0::Point,DoY::Int64,state::Int64)     # constructor
    Base.copy(m::Mosquito)                          # copy
    describe(M::Mosquito,flag=true)                 # alpha text
    show_map_m(state)                               # visualization
    gen_hist_m()                                    # alpha hist selects ABM[t]
    
implements the main model functions like:
    random_walk_mosquito(id::Int64,to::Float64) # random
    walk_mosquito(id::Int64,to::Float64)        # directed
    walks_mosquito(n::Int64, sel::String, to::Float64, random=true) # for loop
    add_m(n::Int64,DoY::Int64,gx::Grid,state::Int64,limit=10.0f0)   # add a new mosquito
    

bird_struct.jl
--------------
implements a mutable struct Bird including some functions about it:
    Bird(p0,DoY,state)    # constructor
    Base.copy(b::Bird)    # copy
    describe(b::Bird)     # alpha text
    show_map_b(b::Bird)   # visualization
    gen_hist_b(t=0)       # alpha hist t selects BRD[t]
 
implements the main model functions like:   
    add_b(n::Int64,DoY::Int64,gx::Grid,state::Int64,limit=3.0f0) # find a place for new bird
    
DWD.jl
------
DWD is a modul to load DWD station as a table
    load_DWD(selector,from=20200101,to=20201231)        # temp
    load_Wind(selector,from=2020010100,to=2020123100)   # power and diretion 
some important functions that are implemented here but used by other modules:
    daylength(dayOfYear,Station)            
    δM(dayOfYear,Station)=1.0-1.0/(1.0+1775.7*exp(1.559*(daylength(dayOfYear,Station)-18.177)))
    
DGL.jl
------
most important function of the simulator:
    simulation()  # solves the jump problem 
Due to this importance and the complexity of DGL.jl a module sensitivity.jl was implemented which emulate the DGL.jl
For the sensitivity.jl exists a README_sensitivity.txt with a detailed description of the modul and of DGL.jl.


