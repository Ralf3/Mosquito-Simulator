# define the region and simulation time
# =====================================

R="Magdeburg"            # Bitterfeld, Dueben, Radebeul, Magdeburg, Trebnitz
_from_=20190101       # year month day for DWD
_to_ = 20201231       # year month day for DWD
start=20200430        # start of simulation
species="pipiens"      # define the species (pipiens, vexans)
birds="kraehe"        # bird: sperling,gruenfink,amsel,elster,habicht,kraehe
master="Jump"         # the branch of the git repository
remark="Abschluss"    # 12 letters limit

# define the parameters of the simulation
# ---------------------------------------
u0b=[300.0,0.0,0.0,0.0,0.0] # start values for birds (SEIRD)
u0m=[200.0,200.0,0.0,50.0] # start values for mosquitoes (LSEI)
steps=200                   # number of simulation steps
infstepIB=0.                # set infstepIB day or 0. ==> no step 
addIB=10                    # default only 1 step
infstepIM=0.                # set infstepIM day or 0. ==> no step
addIM = 10                   # default only 1 step
walkn=1                     # number of walks/day
walk_targets="model"     # the target of the walk s=string!
initial_targets="siedlung"     # where the simulation starts from

# use the utils to make the definition persistent
# ===============================================
global array_of_mosqu=[] 
global array_of_birds=[]  
include("util.jl")

# initial location of  birds
# --------------------------
add_b(Int64(round(u0b[1])),DoY,siedlung250d,1)   # set new birds at random pos
add_b(Int64(round(u0b[3])),DoY,siedlung250d,3)   # set new birds at random pos
# add_infected_b(125,125,60,60,Int64(round(u0b[3])),1,3) # birds i,j

# initial location of mosquitoes
# ------------------------------
add_m(Int64(round(u0mx[2])),DoY,initial_target,2)   # add mosquitoes
add_m(Int64(round(u0mx[3])),DoY,initial_target,3)   # add expected
add_m(Int64(round(u0mx[4])),DoY,initial_target,4)   # add infectious

# add_infected_m(125,125,60,60,Int64(round(u0mx[4])),DoY) # Mosquitoes

# heartbeat of the model calling SIRM and SIRB
# ============================================

u0 = vcat(u0m,u0b)
function simulation3(steps)
    global u0b, u0m, u0, ABM, BRD, array_of_mosqu, array_of_birds, gx, solj
    # switch("siedlung")   # start with siedlung as map
    k=0                  # selector für die ABM/BRD
    for i in 1:steps
        walks_mosquito(walkn,walk_targets,1.0,true)  # direct from the river
        # release_lost_m()                           # remove  lost mosquitoes
        replace_lost_m()                             # replace lost mosquitoes
        
        # --------------------- Do Not Change --------------------------------
        if(i%1==0)
            k+=1                         # set the selector k
            ABM[k]=deepcopy(array_of_mosqu)  # save the array
            h=gen_hist_m()               # address the mosquitoes
            # ------------------mosquitoes---------------------------------
            """
            println(k,
                    "\tsolj2:\t",Int64(round(solj(i+DoY)[2])),"\th[2]:\t",h[2], 
                    "\tsolj3:\t",Int64(round(solj(i+DoY)[3])),"\th[3]:\t",h[3],
                    "\tsolj4:\t",Int64(round(solj(i+DoY)[4])),"\th[4]:\t",h[4])
            """
            add_m(Int64(round(solj(i+DoY)[2]-h[2])),DoY+i,initial_target,2)
            add_m(Int64(round(solj(i+DoY)[3]-h[3])),DoY+i,initial_target,3) 
            add_m(Int64(round(solj(i+DoY)[4]-h[4])),DoY+i,initial_target,4) 
            # ------------------birds--------------------------------------
  
            BRD[k]=deepcopy(array_of_birds)   # copy the array_of_birds
            hb=gen_hist_b()
            
            println(k,
                    " tsolj5: ",Int64(round(solj(i+DoY)[5]))," h[5]: ",hb[1], 
                    " tsolj6: ",Int64(round(solj(i+DoY)[6]))," h[6]: ",hb[2],
                    " tsolj7: ",Int64(round(solj(i+DoY)[7]))," h[7]  ",hb[3],
                    " tsolj8: ",Int64(round(solj(i+DoY)[8]))," h[8]: ",hb[4],
                    " tsolj9: ",Int64(round(solj(i+DoY)[9]))," h[9]: ",hb[5])
            
            add_b(Int64(round(solj(i+DoY)[5]-hb[1])),DoY+i,siedlung250d,1) #S
            add_b(Int64(round(solj(i+DoY)[6]-hb[2])),DoY+i,siedlung250d,2) #E
            add_b(Int64(round(solj(i+DoY)[7]-hb[3])),DoY+i,siedlung250d,3) #I
            add_b(Int64(round(solj(i+DoY)[8]-hb[4])),DoY+i,siedlung250d,4) #R
            add_b(Int64(round(solj(i+DoY)[9]-hb[5])),DoY+i,siedlung250d,5) #D
        end 
    end
    switch("siedlung")
end

function main()
    global steps
    simulation3(steps)
end

main()
solj=solj(DoY:(DoY+steps));  # clean the solj from DGL
df,dfx=save_control();       # save the control.jl and result_?.csv
show_resb()
