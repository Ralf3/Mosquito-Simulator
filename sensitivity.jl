# adapt the pM (prob. transmition of des. from mosquito to the bird)
# adpat the pB (prop. from bird to mosquito 0.125 to 0.2)

# include all packages
# ====================

using Pkg
Pkg.activate("/home/ralf/ABM")
using DifferentialEquations
# using DiffEqJump
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

# include all general parameters
# ==============================

u0m=[100.0,500.0,0.0,50.0]    # start values for mosquitoes (LSEI)
u0b=[300.0,30.0,30.0,0.0,0.0] # start values for birds (SEIRD)
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
master="Jump"                 # the branch of the git repository do not change it:w
remark="testX"                # 12 letters limit
steps=200                     # number of simulation steps
infstepIB=0.                  # set infstepIB day or 0. ==> no step
addIB= 20                     # default only 1 step
infstepIM=0.                  # set infstepIM day or 0. ==> no step
addIM = 10                    # default only 1 step
walkn=1                       # number of walks/day

# include parameter and weather files
# ===================================

include("bird_parameter.jl")
include("DWD.jl")

# load the DWD-data
# =================
DWD_data=Dict()
DWD_data["Bitterfeld"]=Bitterfeld
DWD_data["Trebnitz"]=Trebnitz
DWD_data["Dueben"]=Dueben
DWD_data["Radebeul"]=Radebeul
DWD_data["Trebnitz"]=Trebnitz
DWD_data["Magdeburg"]=Magdeburg
LT=load_DWD(DWD_data[R], _from_, _to_);

# Bird parameter for sensitivity/uncertainty
# ==========================================

αB = 0.4        # recovery rate of bird
γB = 1.0        # incubition rate of bird
νB = 0.7        # portion of dead bird
KB = 500.0      # carrying capacity bird (will not be used in simulate)
mB = 0.00034    # mortality rate of bird original
λ = 1.0         # additional factor for small birds exposed_rateB

# select the bB according to the birds
if(birds == "sperling")
    bB=sperling_bB()
    KB=1000.0
end
if(birds == "gruenfink")
    bB=gruenfink_bB()
    KB=900.0
end
if(birds == "amsel")
    bB=amsel_bB()
    KB=1500.0
end
if(birds == "elster")
    bB=elster_bB()
    KB=800.0
end
if(birds == "kraehe")
    bB=kraehe_bB()
    KB=800.0
end
if(birds == "habicht")
    bB=habicht_bB()
    KB=200.0
end

if(birds=="sperling")
    #rate_a=exp(log(0.28)/365) # katrin
    rate_a=exp(log(0.26)/365)  # ralf
    mB = 1-rate_a
    λ = 0.5 # additional factor for small birds  exposed_rateB
end
if(birds=="gruenfink")
    #rate_a=exp(log(0.28)/365) # katrin
    rate_a=exp(log(0.26)/365)  # ralf 
    mB = 1-rate_a
    λ = 0.5 # additional factor for small birds  exposed_rateB
end
if(birds=="habicht")
    rate_a=exp(log(0.66)/365)
    mB = 1-rate_a
end
if(birds=="elster")
    rate_a=exp(log(0.84)/365)
    mB = 1-rate_a
    νB = 0.43    # according to de Oya et al. (2018)
end
if(birds=="kraehe")
    rate_a=exp(log(0.83)/365)
    mB = 1-rate_a
end
if(birds=="amsel")
    rate_a=exp(log(0.27)/365)
    mB = 1-rate_a
    νB = 0.9     # is killed by the WNV 
end

# from data.jl 
# vexans dict
# ============

VD=Dict("Dueben" => 0.65f0, "Bitterfeld" => 0.8f0,
        "Magdeburg" => 0.95f0, "Radebeul" => 0.95f0, "Trebnitz" => 0.65f0)
# pipiens dict
PD=Dict("Dueben" => 0.65f0, "Bitterfeld" => 0.9f0,
        "Magdeburg" => 0.9f0, "Radebeul" => 0.95f0, "Trebnitz" => 0.75f0)
# select dict
SD=Dict("vexans" => VD, "pipiens" => PD)

KM = 10000*SD[species][R] # KM is calculated from model

Station=Dict("Dueben" => Dueben, "Bitterfeld" => Bitterfeld,
             "Magdeburg" => Magdeburg, "Radebeul" => Radebeul,
             "Trebnitz" => Trebnitz) 

# Paramter from DGL.jl
# ====================

# transfer from bird to mosquito and vise versa USUTU
pM = 1.0     # prob. transmition of des. from mosquito to the bird
pB = 0.20    # prop. from bird to mosquito 0.125 to 0.2
ϕB = 30      # Mosquito-to-bird ratio

param=OrderedDict()
param["αB"]=αB   # recovery rate of bird
param["λ"]=λ     # additional factor for small birds  exposed_rateB
param["γB"]=γB   # incubition rate of bird
param["νB"]=νB   # is killed by the WNV 
param["mB"]=mB   # mortality of a bird 
param["KM"]=KM   # carrying capacity mosquito
param["pM"]=pM   # prob. transmission WNV from mosquito to the bird
param["pB"]=pB   # prop. transmission WNV from bird to mosquito 0.125 to 0.2
param["ϕB"]=ϕB   # Mosquito-to-bird ratio

biting_rate(T)=0.344/(1+1.231*exp(-0.184*(T-20)))
bL(T)=2.325*biting_rate(T)
bM(T)=bL(T)/10
mL(T)=(0.0025*T^2-0.094*T+1.0257)
mM(T)=0.1*mL(T)
βM(T)=biting_rate(T)*pM
βB(T)=biting_rate(T)*pB

function γM(T)
    if T<=15.0
        return 0
    end
    return 0.0093*T-0.1352
end

# Jump Model definition
# =====================

function SEIRJ(du,u,p,t)
    LM,SM,EM,IM,SB,EB,IB,RB,DB = u
    Nm=SM+EM+IM
    NB=SB+RB
    KM,KB = p
    T=LT[t]
    DoY=Int64(round(t))
    du[1] = (bL(T)*δM(DoY,Station[R])*Nm-mL(T)*LM)*(1.0-LM/KM) # dLarvae
end

# jump model mosquito
# ===================

birth_rateM(u,p,t) = bM(LT(t))*u[1]
function birth_affectM!(integrator)
    integrator.u[1] -= 1
    integrator.u[2] += 1
end
birth_jumpM = VariableRateJump(birth_rateM, birth_affectM!)

exposed_rateM(u,p,t) = δM(t,Station[R])*biting_rate(LT(t))*pB*u[7]/p[2]*u[2]
function exposed_affectM!(integrator)
    integrator.u[2] -= 1
    integrator.u[3] += 1
end
exposed_jumpM=VariableRateJump(exposed_rateM,exposed_affectM!)

infected_rateM(u,p,t) = γM(LT(t))*u[3]
function infected_affectM!(integrator)
    integrator.u[3] -= 1
    integrator.u[4] += 1
end
infected_jumpM = VariableRateJump(infected_rateM,infected_affectM!)

deadS_rateM(u,p,t) = mM(LT(t))*u[2]
function deadS_affectM!(integrator)
    integrator.u[2] -= 1
end
deadS_jumpM = VariableRateJump(deadS_rateM,deadS_affectM!)

deadE_rateM(u,p,t) = mM(LT(t))*u[3]
function deadE_affectM!(integrator)
    integrator.u[3] -= 1
end
deadE_jumpM = VariableRateJump(deadE_rateM,deadE_affectM!)

deadI_rateM(u,p,t) = mM(LT(t))*u[4]
function deadI_affectM!(integrator)
    integrator.u[4] -= 1
end
deadI_jumpM = VariableRateJump(deadI_rateM,deadI_affectM!)

# jump model bird
# ===============

birth_rateB(u,p,t) = bB(t)*u[5] 
function birth_affectB!(integrator)
    integrator.u[5] += 1
end
birth_jumpB = VariableRateJump(birth_rateB, birth_affectB!)

exposed_rateB(u,p,t) = λ*(δM(t,Station[R])*biting_rate(LT(t))*pM*u[4]/p[1])*ϕB*u[5]
function exposed_affectB!(integrator)
    integrator.u[5] -= 1
    integrator.u[6] += 1
end
exposed_jumpB=VariableRateJump(exposed_rateB,exposed_affectB!)

infected_rateB(u,p,t) = γB*u[6]
function infected_affectB!(integrator)
    if(integrator.u[6]>=1)
        integrator.u[6] -= 1
        integrator.u[7] += 1
    end
end
infected_jumpB = ConstantRateJump(infected_rateB,infected_affectB!)

recover_rateB(u,p,t) = (1-νB)*αB*u[7]
function recover_affectB!(integrator)
    integrator.u[7] -= 1
    integrator.u[8] += 1
end
recover_jumpB = ConstantRateJump(recover_rateB,recover_affectB!)

dead_rateB(u,p,t) = νB*αB*u[7] 
function dead_affectB!(integrator)
    if(integrator.u[7]>=1)
        integrator.u[7] -= 1
        integrator.u[9] += 1
    end
end
dead_jumpB = ConstantRateJump(dead_rateB,dead_affectB!)

deadS_rateB(u,p,t) = mB*u[5]
function deadS_affectB!(integrator)
    #if(integrator.u[5]>=1)
        integrator.u[5] -= 1
        integrator.u[9] += 1
    #end
end
deadS_jumpB = ConstantRateJump(deadS_rateB,deadS_affectB!)

deadE_rateB(u,p,t) = mB*u[6]
function deadE_affectB!(integrator)
    integrator.u[6] -= 1
    integrator.u[9] += 1
end
deadE_jumpB = ConstantRateJump(deadE_rateB,deadE_affectB!)

deadI_rateB(u,p,t) = mB*u[7]
function deadI_affectB!(integrator)
    #if(integrator.u[7]>=1)
        integrator.u[7] -= 1
        integrator.u[9] += 1
    #end
end
deadI_jumpB = ConstantRateJump(deadI_rateB,deadI_affectB!)

deadR_rateB(u,p,t) = mB*u[8]
function deadR_affectB!(integrator)
    integrator.u[8] -= 1
    integrator.u[9] += 1
end
deadR_jumpB = ConstantRateJump(deadR_rateB,deadR_affectB!)

# jump model for set  
# ==================

set_rateM(u,p,t) = 1.0
function set_affectM!(integrator)
    global infstepIM
    if(infstepIM==0.0)
        return
    end
    if(abs(integrator.t-(DoY+infstepIM))<=1.0)  
        integrator.u[4] += addIM
        println("M: ",ceil(integrator.t), " :", DoY+infstepIM, "  ",
                addIM, " ==> ", integrator.u[4])
        infstepIM=0.0
    end
end
set_jumpM = ConstantRateJump(set_rateM,set_affectM!)

set_rateB(u,p,t) = 1.0
function set_affectB!(integrator)
    global infstepIB
    if(infstepIB==0.0)
        return
    end
    if(abs(integrator.t-(DoY+infstepIB))<=1.0)  
        integrator.u[6] += addIB
        println("B: ",ceil(integrator.t), " :", DoY+infstepIB, "  ",
                addIB, " ==> ", integrator.u[6])
        infstepIB=0
    end
end
set_jumpB = ConstantRateJump(set_rateB,set_affectB!)


# Dates for simulaton
# ===================

Dates.dayofyear(s::Int64)=dayofyear(Date(string(start),dateformat"yyyymmdd"))
get_Year(s::Int64)=year(Date(string(start),dateformat"yyyymmdd"))
_year_=get_Year(start) # is used by data
DoY=dayofyear(start)          # start DoY
tspan=(Float64(DoY),Float64(DoY+steps))
# ----- calculated values for control do not touch it -------------

# Do the analysis
# ===============

function reset_simu()
    global DoY,tspan,RES,sols,df1
    sols=[]
    p=[KM,KB]
    u0 = vcat(u0mx,u0bx)
    param["u1"]=u0[1]
    param["u2"]=u0[2]
    param["u3"]=u0[3]
    param["u4"]=u0[4]
    param["u5"]=u0[5]
    param["u6"]=u0[6]
    param["u7"]=u0[7]
    param["u8"]=u0[8]
    param["u9"]=u0[9]
    prob = ODEProblem(SEIRJ,u0,tspan,p)
    jump_prob = JumpProblem(prob, Direct(),
                            birth_jumpM,
                            exposed_jumpM,
                            infected_jumpM,
                            deadS_jumpM,
                            deadE_jumpM,
                            deadI_jumpM,
                            birth_jumpB,
                            exposed_jumpB,
                            infected_jumpB,
                            recover_jumpB,
                            dead_jumpB,
                            deadS_jumpB,
                            deadE_jumpB,
                            deadI_jumpB,
                            deadR_jumpB,
                            set_jumpM,
                            set_jumpB
                            )                        


    # sensitivity analysis based on many simulations
    # ==============================================
    solj = solve(jump_prob, Tsit5())
    solj = solj(DoY:(DoY+steps))
    RES=OrderedDict()
    RES["Lm"]=mean(solj[1,:])
    RES["SMm"]=mean(solj[2,:])
    RES["EMm"]=mean(solj[3,:])
    RES["IMm"]=mean(solj[4,:])
    RES["SBm"]=mean(solj[5,:])
    RES["EBm"]=mean(solj[6,:])
    RES["IBm"]=mean(solj[7,:])
    RES["RBm"]=mean(solj[8,:])
    RES["DBm"]=mean(solj[9,:])
    df1=DataFrame(RES)
    sols=push!(sols,solj)  # save some solutions for the first and the last run
    true
end

function simulate(n,dx)
    """ adapt the parameters """
    reset_simu() 
    global param, df1,sols, ϕB,pB,λ,γB,νB,mB,KM,KB,pM,pB,ϕB
    for i in 1:n
        αB=param["αB"]   # recovery rate of bird
        λ=param["λ"]     # additional factor for small birds  exposed_rateB
        γB=param["γB"]   # incubition rate of bird
        νB=param["νB"]   # is killed by the WNV 
        mB=param["mB"]   # mortality of a bird 
        KM=param["KM"]   # carrying capacity mosquito
        pM=param["pM"]   # prob. transmission WNV from mosquito to the bird
        pB=param["pB"]   # prop. transmission WNV from bird to mosquito
        ϕB=param["ϕB"]   # Mosquito-to-bird ratio
	# αB=param["αB"]*(1+dx*randn())
        # νB=param["νB"]*(1+dx*randn())
        # αB=param["αB"]*(0.5+dx*i/n)
        # νB=param["νB"]*(0.5+dx*i/n)
	# println(i,"\t",αB,"\t",νB)
        RES=OrderedDict()
        DoY=dayofyear(start)          # start DoY
        tspan=(Float64(DoY),Float64(DoY+steps))
        p=[KM,KB]
        u0 = vcat(u0mx,u0bx)
        u0[2]=param["u2"]*(0.5+2*i/n)
        println(i,"\t",u0)
        prob = ODEProblem(SEIRJ,u0,tspan,p)
        jump_prob = JumpProblem(prob, Direct(),
                                birth_jumpM,
                                exposed_jumpM,
                                infected_jumpM,
                                deadS_jumpM,
                                deadE_jumpM,
                                deadI_jumpM,
                                birth_jumpB,
                                exposed_jumpB,
                                infected_jumpB,
                                recover_jumpB,
                                dead_jumpB,
                                deadS_jumpB,
                                deadE_jumpB,
                                deadI_jumpB,
                                deadR_jumpB)
        solj = solve(jump_prob, Tsit5())
        solj = solj(DoY:(DoY+steps))
        if(i==n)
            push!(sols,solj)
        end
        # println(i)
        RES["Lm"]=mean(solj[1,:])
        RES["SMm"]=mean(solj[2,:])
        RES["EMm"]=mean(solj[3,:])
        RES["IMm"]=mean(solj[4,:])
        RES["SBm"]=mean(solj[5,:])
        RES["EBm"]=mean(solj[6,:])
        RES["IBm"]=mean(solj[7,:])
        RES["RBm"]=mean(solj[8,:])
        RES["DBm"]=mean(solj[9,:])
 	push!(df1,RES)
    end
    df1
end

# Visualization of solj
# =====================

function show_resm()
    # result of simulation mosquitoes
    plot(sols[1][1,:],xlabel = "t [d]")
    p1=plot!(sols[2][1,:],xlabel = "t [d]")
    plot(sols[1][2,:],xlabel = "t [d]")
    p2=plot!(sols[2][2,:],xlabel = "t [d]")
    plot(sols[1][3,:],xlabel = "t [d]")
    p3=plot!(sols[2][3,:],xlabel = "t [d]")
    plot(sols[1][4,:],xlabel = "t [d]")
    p4=plot!(sols[2][4,:],xlabel = "t [d]")
    plot(p1, p2, p3, p4, layout = (2, 2))
end

function show_resb()
    # result of simulation birds
    plot(sols[1][5,:],xlabel = "t [d]")
    p1=plot!(sols[2][5,:],xlabel = "t [d]")
    plot(sols[1][7,:],xlabel = "t [d]")
    p2=plot!(sols[2][7,:],xlabel = "t [d]")
    plot(sols[1][8,:],xlabel = "t [d]")
    p3=plot!(sols[2][8,:],xlabel = "t [d]")
    plot(sols[1][9,:],xlabel = "t [d]")
    p4=plot!(sols[2][9,:],xlabel = "t [d]")
    plot(p1, p2, p3, p4, layout = (2, 2))
end

function compare_sols(nr)
    global sols
    plot(sols[1][nr,:],xlabel = "t [d]")
    plot!(sols[2][nr,:],xlabel = "t [d]")
end
    
"""
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
"""
True
