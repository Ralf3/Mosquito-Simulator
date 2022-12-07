using Pkg
Pkg.activate("/home/ralf/ABM")
using DifferentialEquations
using Distributions
using Distributed
using DataStructures
using DataFrames
using CSV
using Random
using Interpolations
using SpecialFunctions
using Dates
using Base
using Images
using ImageView
using ImageFiltering 
using Colors
using FileIO
using ImageIO
using Statistics
using StatsPlots
using NPZ
using JLD2
using Dierckx
# using Plots

# save the infstepIM and infstepIB
infstepIM0=Int64(infstepIM)
infstepIB0=Int64(infstepIB)

# plot the resm and resb
# ======================

function show_resm()
    global solj   # result of simulation mosquitoes
    p1=plot(solj[1,:],xlabel = "t [d]")
    p2=plot(solj[2,:],xlabel = "t [d]")
    p3=plot(solj[3,:],xlabel = "t [d]")
    p4=plot(solj[4,:],xlabel = "t [d]")
    plot(p1, p2, p3, p4, layout = (2, 2), label=["L" "S" "E" "I"])
end

function show_resb()
    global solj   # result of simulation birds
    p1=plot(solj[5,:],xlabel = "t [d]")
    p2=plot(solj[7,:],xlabel = "t [d]")
    p3=plot(solj[8,:],xlabel = "t [d]")
    p4=plot(solj[9,:],xlabel = "t [d]")
    plot(p1, p2, p3, p4, layout = (2, 2), label=["S" "I" "R" "D"])
end

# Statistical analysis of the DGL results
# =======================================
# switch to databases.jl

# Dates for simulaton
# ===================
#
# start=20190430       # an integer 
Dates.dayofyear(s::Int64)=dayofyear(Date(string(start),dateformat"yyyymmdd"))
get_Year(s::Int64)=year(Date(string(start),dateformat"yyyymmdd"))
_year_=get_Year(start) # is used by data

# ----- calculated values for control do not touch it -------------
DoY=dayofyear(start)          # start DoY
u0bx=copy(u0b)                # copy start values for birds
u0mx=copy(u0m)                # copy start values for mosquitoes

# include the modules
# ===================
include("DWD.jl")    # base data to handle the DWD
include("grid.jl")   # base data to handle the OSM
include("data.jl")   # loads the spatial OSM data and the LT from DWD
include("DGL.jl")    # load the paramter of the DGL models SIRM SIRB
include("bird_struct.jl")     # load the bird ABM
include("mosquito_struct.jl") # load the mosquito ABM
include("database.jl")        # store the results in simulation.csv

# set the initial_target and walk_target for DataFrame
# ====================================================
initial_target=model250d
if(initial_targets=="river")
    initial_target=river250d 
end
if(initial_targets=="siedlung")
    initial_target=siedlung250d
end

# DOKUMENTATION ERGEBNISSE
# ========================

array_of_mosqu=[]    # from control.jl
array_of_birds=[]    # from control.jl
ABM=Dict()        # dict to store the mosquito movement
BRD=Dict()        # dict to store the birds
function switch(sel::Int64) # 1,2,...,10 
    global ABM,BRD,array_of_birds,array_of_mosqu
    array_of_birds=BRD[sel];
    array_of_mosqu=ABM[sel];
    reset_stay()
    println("nmos:", length(array_of_mosqu), " birds: ", length(array_of_birds))
end

function reset_stay()
    # reset the state=stay to state=suspected for all mosquitoes
    for k in 1:length(array_of_mosqu)
        if(array_of_mosqu[k].state==stay)
            array_of_mosqu[k].state=suspected
        end # if
    end # for
end # reset_stay()
    
gx=siedlung250d   # define the grid 
gxx=siedlung250dx # define the dx
gxy=siedlung250dy # define the dy
function switch(sel::String) # river,model,siedlung
    global gx,gxx,gxy,array_of_mosqu
    if(sel=="river")
        gx=river250d
        gxx=river250dx
        gxy=river250dy
    end
    if(sel=="model")
        gx=model250d
        gxx=model250dx
        gxy=model250dy
    end
    if(sel=="siedlung")
        gx=siedlung250d
        gxx=siedlung250dx
        gxy=siedlung250dy
    end
end
    
