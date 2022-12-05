# define the parameter for the DGLs
# =================================

include("bird_parameter.jl")
lat=DWD_data[R]["lat"]
lon=DWD_data[R]["lon"]

αB = 0.4        # recovery rate of bird
γB = 1.0        # incubition rate of bird
νB = 0.7        # portion of dead bird
KB = 500.0      # carrying capacity bird

# bird parameter
mB = 0.00034    # mortality rate of bird original
λ = 1.0         # additional factor for small birds exposed_rateB
# assumption the number of birds: NoB(1)==NoB(365)
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

# mosquito part of parameter 
# this should be paramters for the modelling interface 
KM = 10000*SD[species][R]  # carrying capccity of mosquitoes depends on habitat

# transfer from bird to mosquito and vise versa USUTU
pM = 1.2     # prob. transmition of des. from mosquito to the bird
pB = 0.2     # prop. from bird to mosquito 0.125 to 0.2
# transfer from bird to mosquito and vise versa WNV
pM = 1.0     # prob. transmition of des. from mosquito to the bird
pB = 0.20    # prop. from bird to mosquito 0.125 to 0.2
ϕB = 30      # Mosquito-to-bird ratio

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

#function bB(d)
#    """ birth rate of the bird """
#    gx=Gamma{Float64}(42.710230749049124, 3.050791290864214)
#    pdf(gx,d)
#end

# original function
#function bB1(d)
#    α=86.4
#    β=1.4
#    ((d/β)^(α-1.0)*exp(-d/β))/(β*gamma(α))
#end

# select the dB according to the birds
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
    KB=1000.0
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

y(x)=1.2e7*pdf(Gamma(100.0),x-150)  # estimate the number of mosquitoes (nmos)
# y(x)=1.2e7*cdf(Gamma(100.0),x-150)  # estimate the number of mosquitoes (nmos)

# Mosquito model SIRM
# -------------------
function SIRM(du,u,p,t)
    LM,SM,EM,IM = u
    Nm=SM+EM+IM
    KM,KB,IB=p
    T=LT[t]
    DoY=Int64(round(t))
    lBM=δM(DoY,DWD_data[R])*biting_rate(T)*pB*IB/KB # lamda Bird Mosquito
    du[1] = (bL(T)*δM(DoY,DWD_data[R])*Nm-mL(T)*LM)*(1.0-LM/KM)-bM(T)*LM # dLarvae
    du[2] = -lBM*SM + bM(T)*LM - mM(T)*SM # dSM 
    du[3] = lBM*SM - γM(T)*EM - mM(T)*EM  # dEM
    du[4] = γM(T)*EM - mM(T)*IM           # dIM
end

# Bird model SIRB
# ---------------

function SIRB(du,u,p,t)
    """
    rB=bB-mB # (rubel 2008)
    dSB/dt = rb*(1-NB/KB)*NB=(bB-(bB-mB)*NB/KB)*NB-mB*SB
    """
    SB,EB,IB,RB=u
    NB=SB+RB
    DoY=Int64(round(t))
    T=LT[t]
    IM,KM=p  # connect to mosquito model
    # most important parameter from mosqito to birds
    lMB= δM(DoY,DWD_data[R])*biting_rate(T)*pM*IM/KM*ϕB # lambda Mosquito Bird
    # lMB= δM(DoY)*biting_rate(T)*ldata[t]/KM*0.04  #  anpassen
    du[1] = (bB(DoY)-(bB(DoY)-mB)*NB/KB)*NB-lMB*SB-mB*SB # SB
    du[2] = lMB*SB - γB*EB - mB*EB          # dEB
    du[3] = γB*EB - αB*IB - mB*IB           # dIB
    du[4] = (1-νB)*αB*IB - mB*RB            # dRB
    du[5] = νB*αB*IB + mB*EB +mB*IB+mB*RB   # dDB
end

# Combination Mosquito Bird DEQ
# =============================

function SEIR(du,u,p,t)
    LM,SM,EM,IM,SB,EB,IB,RB,DB = u
    Nm=SM+EM+IM
    NB=SB+RB
    KM,KB = p
    T=LT[t]
    DoY=Int64(round(t))
    lBM=δM(DoY,DWD_data[R])*biting_rate(T)*pB*IB/KB # lamda Bird Mosquito
    du[1] = (bL(T)*δM(DoY)*Nm-mL(T)*LM)*(1.0-LM/KM)-bM(T)*LM # dLarvae
    du[2] = -lBM*SM + bM(T)*LM - mM(T)*SM # dSM 
    du[3] = lBM*SM - γM(T)*EM - mM(T)*EM  # dEM
    du[4] = γM(T)*EM - mM(T)*IM           # dIM
    lMB= δM(DoY)*biting_rate(T)*pM*IM/KM*ϕB # lambda Mosquito Bird
    # lMB= δM(DoY)*biting_rate(T)*ldata[t]/KM*0.04  #  anpassen
    du[5] = (bB(DoY)-(bB(DoY)-mB)*NB/KB)*NB-lMB*SB-mB*SB # SB
    du[6] = lMB*SB - γB*EB - mB*EB          # dEB
    du[7] = γB*EB - αB*IB - mB*IB           # dIB
    du[8] = (1-νB)*αB*IB - mB*RB            # dRB
    du[9] = νB*αB*IB + mB*EB +mB*IB+mB*RB   # dDB
end

#tspan=(Float64(DoY),Float64(DoY+steps))
#p=[KM,KB]
#u0 = vcat(u0m,u0b)
#prob = ODEProblem(SEIR,u0,tspan,p)
#sol2 = solve(prob,Tsit5(),u0=u0,abstol=1e-8,reltol=1e-8,saveat=1.0)

# Jump Model for combination of Mosquito and Bird
# ===============================================

function SEIRJ(du,u,p,t)
    LM,SM,EM,IM,SB,EB,IB,RB,DB = u
    Nm=SM+EM+IM
    NB=SB+RB
    KM,KB = p
    T=LT[t]
    DoY=Int64(round(t))
    du[1] = (bL(T)*δM(DoY,DWD_data[R])*Nm-mL(T)*LM)*(1.0-LM/KM) # dLarvae
end

# jump model mosquito
birth_rateM(u,p,t) = bM(LT(t))*u[1]
function birth_affectM!(integrator)
    integrator.u[1] -= 1
    integrator.u[2] += 1
end
birth_jumpM = VariableRateJump(birth_rateM, birth_affectM!)

exposed_rateM(u,p,t) = δM(t,DWD_data[R])*biting_rate(LT(t))*pB*u[7]/p[2]*u[2]
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
birth_rateB(u,p,t) = bB(t)*u[5] 
function birth_affectB!(integrator)
    integrator.u[5] += 1
end
birth_jumpB = VariableRateJump(birth_rateB, birth_affectB!)

# exposed_rateB(u,p,t) = (δM(t)*biting_rate(LT(t))*pB*u[4]/KB)*u[5]
exposed_rateB(u,p,t) = λ*(δM(t,DWD_data[R])*biting_rate(LT(t))*
    pM*u[4]/p[1])*ϕB*u[5]
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
#
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

# unsing a very simple procedure due to callback did not work
# ===========================================================

tspan=(Float64(DoY),Float64(DoY+steps))
p=[KM,KB]
u0 = vcat(u0m,u0b)
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

function simulation()
    global infstepIB,infstepIM,DoY,addIB,addIM,jump_prob,u0
    if(infstepIM<=DoY && infstepIB<=DoY)
        solj = solve(jump_prob, Tsit5())
        return solj
    end
    if(infstepIM>DoY)
        tspan1=(Float64(DoY),Float64(infstepIM))
        prob1= ODEProblem(SEIRJ,u0,tspan1,p) # new definition of prob
        jump_prob1=JumpProblem(prob1, Direct(),
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
                        deadE_jumpB
                        )                        

        solj1 = solve(jump_prob1, Tsit5())
        tspan2=(Float64(infstepIM),Float64(DoY+steps))
        solj1.u[end][4]+=addIM
        u1=solj1.u[end]
        prob2= ODEProblem(SEIRJ,u1,tspan2,p)
        jump_prob2=remake(jump_prob1,prob=prob2)
        solj2 = solve(jump_prob2, Tsit5())
        return solj1,solj2
    else
        tspan1=(Float64(DoY),Float64(infstepIB))
        prob1= ODEProblem(SEIRJ,u0,tspan1,p) # new definition of prob
        jump_prob1=remake(jump_prob,prob=prob1)
        solj1 = solve(jump_prob, Tsit5())
        tspan2=(Float64(infstepIB),Float64(DoY+steps))
        solj.u[end][4]+=addIB
        u1=solj.u[end]
        prob2= ODEProblem(SEIRJ,u1,tspan2,p)
        jump_prob1=remake(jump_prob1,prob=prob2)
        solj2 = solve(jump_prob, Tsit5())
        return solj1,solj2
    end
end



solj = solve(jump_prob, Tsit5())
# solj = solve(jump_prob, Tsit5(), callbacks=cbS)
# jprob.prob.u0 .= [499,1,0]?
#sol1,sol2=simulation()
