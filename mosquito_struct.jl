# data type of an mosquito
# ========================
mutable struct Mosquito
    p0::Point     # point of birth; point of infection
    px::Point     # position in the the region
    DoY::Int64    # day of birth
    state::Int64  # 1,2,3,4, stay=5, lost=6, dead=-999
                  # 12,13,14 ZS, ZE, ZI
    ostate::Int64 # state before lost or dead
    walk::Vector{Any} # save the steps
end

LM=1
SM=2
EM=3
IM=4
stay=5
lost=6
dead=-999

# constructor of an mosquito
function Mosquito(p0::Point,DoY::Int64,state::Int64)
    px=copy(p0)
    return Mosquito(p0,px,DoY,state,state,[]) # state and ostate
end
# copy 
function Base.copy(m::Mosquito)
    p0=copy(m.p0)
    px=copy(m.px)
    walk=copy(m.walk)
    return Mosquito(p0,px,m.DoY,m.state,ostate,walk)
end
# description
function describe(M::Mosquito,flag=true)
    println("Place of birth i: ", M.p0.i, " j: ",M.p0.j)
    println("Place of end simulation i: ", M.px.i, " j: ",M.px.j)
    println("DoY:",M.DoY, " length of walk: ", length(M.walk),
            " end state:", M.state)
end

function show_map_m(state)
    """ search the complete array """
    global array_of_mosqu, gx
    g1=copy(gx)  # copy of the grid
    z=round(maximum(g1.grid))+1
    counter=0
    end1=length(array_of_mosqu)
    for i in 1:end1
        if(array_of_mosqu[i].state==state)
            g1.grid[array_of_mosqu[i].px.i,array_of_mosqu[i].px.j]=z
            counter+=1
        end # if
    end # for
    println(counter)
    show(g1)  # show the birds
end # show_map

# very important function to check the states of all mosquitoes
# -------------------------------------------------------------
function gen_hist_m()
    """ select from active array_of_mosqu """
    Hist=Dict()
    Hist[1]=0  # L
    Hist[2]=0  # S
    Hist[3]=0  # E
    Hist[4]=0  # I
    Hist[5]=0  # stay
    Hist[6]=0  # lost
    Hist[10]=0 # dummy
    Hist[12]=0 # ZS 
    Hist[13]=0 # ZE
    Hist[14]=0 # ZI
    Hist[dead]=0
    for i in 1:length(array_of_mosqu)
        Hist[array_of_mosqu[i].state]+=1
    end
    Hist
end

function gen_hist_m(sel::Int64)
    """ select from an array stored in ABM """
    Hist=Dict()
    Hist[1]=0  # L
    Hist[2]=0  # S
    Hist[3]=0  # E
    Hist[4]=0  # I
    Hist[5]=0  # stay
    Hist[6]=0  # lost
    Hist[10]=0 # dummy
    Hist[12]=0 # ZS 
    Hist[13]=0 # ZE
    Hist[14]=0 # ZI
    Hist[dead]=0
    for i in 1:length(ABM[sel])
        Hist[ABM[sel][i].state]+=1
    end
    Hist
end

function select_infectious_m(t=0)
    array_of_m=array_of_mosqu
    if(t!=0)
        array_of_m=ABM[t]
    end
    l_mosqu=size(array_of_m)[1]
    infected=[]
    for i in 1:l_mosqu
        if(array_of_m[i].state==IM)
            push!(infected,i)
        end
    end
    infected
end

function select_exposed_m(t=0)
    array_of_m=array_of_mosqu
    if(t!=0)
        array_of_m=ABM[t]
    end
    l_mosqu=size(array_of_m)[1]
    exposed=[]
    for i in 1:l_mosqu
        if(array_of_m[i].state==EM)
            push!(exposed,i)
        end
    end
    exposed
end

function select_mosquitoes(t=0)
    array_of_m=array_of_mosqu
    if(t!=0)
        array_of_m=ABM[t]
    end
    l_mosqu=size(array_of_m)[1]
    mosquitoes=[]
    for i in 1:l_mosqu
        if(array_of_m[i].state==SM)
            push!(mosquitoes,i)
        end
    end
    mosquitoes
end

# some simple simulation functions
# ================================
function random_step_possible(id::Int64,to::Float64)
    """ check if the mosquito with id is able to perform a step """
    global array_of_mosqu,walk_targets,gx  # gx by util.jl
    M=array_of_mosqu[id]
    if(M.state==lost || M.state==dead || M.ostate==lost || M.ostate==dead)
        return dead
    end
    nrows=gx.nrows 
    ncols=gx.ncols
    if(M.px.i<=1 || M.px.i>=nrows || M.px.j<=1 || M.px.j>=ncols)
        M.ostate=M.state
        M.state=lost            # if outside then stop the calc
        return lost       
    end
    if(to<0.0 && walk_targets=="model") # if to<0 then use model as target
        if(model250d.grid[M.px.i,M.px.j]<=2.0f0)
            return stay
        end
    end
    if(to<0.0 && walk_targets!="model") # if to<0 then use siedlung as target
        if(siedlung.grid[M.px.i,M.px.j]<=2.0f0)
            return stay
        end
    end
    if(to>0.0)                  # if to>0 then use gx as target
        if(gx.grid[M.px.i,M.px.j]<=2.0f0) 
            return stay
        end
    end
    if(to<0.0 && walk_targets!="model") # if to<0 then use siedlung as target
        if(siedlung.grid[M.px.i,M.px.j]<=2.0f0)
            return stay
        end
    end
    if(to>0.0)                  # if to>0 then use gx as target
        if(gx.grid[M.px.i,M.px.j]<=2.0f0) 
            return stay
        end
     return SM # is used as true value for can fly 
    end # if
end # function random_step_possible

function random_walk_mosquito(id::Int64,to::Float64,region="siedlung")
    """
    region ∈ ["siedlung","model",river"]
    random step with direction to the river or from the the river
    to==1.0 ==> from to==-1.0 ==> to   
    """
    global array_of_mosqu,walk_targets
    if(random_step_possible(id,to)!=SM)
        return
    end
    if(region=="siedlung")
        g1=siedlung250d
        g1x=siedlung250dx
        g1y=siedlung250dy
    end
     if(region=="model")
        g1=model250d
        g1x=model250dx
        g1y=model250dy
     end
    if(region=="river")
        g1=river250d
        g1x=river250dx
        g1y=river250dy
     end
    # random walk
    nrows=g1.nrows 
    ncols=g1.ncols
    M=array_of_mosqu[id]
    di=0 # set the gradient direction i
    dj=0 # set the gradient direction j
    r=0  # define the resulting direction
    prob1=ones(size(direction)[1])  # all direction are equal
    di=-Int64(round(to*g1y.grid[M.px.i,M.px.j]))
    dj=-Int64(round(to*g1x.grid[M.px.i,M.px.j]))
    # println(di," ",dj)
    for i in 1:9
        if(direction[i,1]==di && direction[i,2]==dj)
            prob1[i]+=1   # add a higer prob to the selected dir
            break
        end
    end
    prob1./=sum(prob1)              # normalize it
    prob1=cumsum(prob1)             # random walk with prefered dir
    if(M.state==SM || M.state==EM || M.state==IM) # larve can not flight
        r=1
        sel=rand()               
        for i in 1:9
            if(prob1[i]>=sel)       # select an random direction
                r=i
                break
            end
        end
        M.px.i+=direction[r,1]
        M.px.j+=direction[r,2]
        push!(M.walk,[M.px.i,M.px.j])
        if(M.px.i<1 || M.px.i>=nrows || M.px.j<1 || M.px.j>=ncols)
            M.ostate=M.state
            M.state=lost        # if outside then stop the calc
            return lost       
        end # if
    end # 100m
end # random_walk_mosquito

# directed walk to the target
# ---------------------------
function walk_mosquito(id::Int64,to::Float64,region="siedlung")
    global array_of_mosqu
    if(random_step_possible(id,to)!=SM)
        return
    end
    if(region=="siedlung")
        g1=siedlung250d
        g1x=siedlung250dx
        g1y=siedlung250dy
    end
     if(region=="model")
        g1=model250d
        g1x=model250dx
        g1y=model250dy
     end
    if(region=="river")
        g1=river250d
        g1x=river250dx
        g1y=river250dy
     end
    M=array_of_mosqu[id]
    nrows=g1.nrows
    ncols=g1.ncols
    # directed walk
    if(M.state==SM || M.state==EM || M.state==IM)  # larve can not flight
        pxi=Int64(round(to*g1y.grid[M.px.i,M.px.j]))
        pxj=Int64(round(to*g1x.grid[M.px.i,M.px.j]))
        M.px.i-=pxi
        M.px.j-=pxj
        if(M.px.i<1 || M.px.i>=nrows || M.px.j<1 || M.px.j>=ncols)
            M.ostate=M.state
            M.state=lost        # if outside then stop the calc
            return lost
        end # if
        push!(M.walk,[M.px.i,M.px.j])
    end # 
end # walk_mosquito

# control function to start walk
# ------------------------------
function walks_mosquito(n::Int64, sel::String, to::Float64, random=true)
    """ walk all mosquitoes n times
        
    """
    global array_of_mosqu
    # aendern um allgemein gueltig zu sein !!!
    switch(sel)   # river,model,siedlung
    if(random==false)
        for i in 1:n
            for k in 1:length(array_of_mosqu)
                walk_mosquito(k,to,walk_targets)
            end # for
        end #for
    else
        for i in 1:n
            for k in 1:length(array_of_mosqu)
                if(rand()>=0.3)
                    random_walk_mosquito(k,to,"siedlung") # random walk
                else
                    random_walk_mosquito(k,to,"model")    # random walk
                end
            end # for
        end #for
    end #if
end # function walks_mosquito

# handle lost mosquitoes
# ----------------------
function find_lost_dead()
    global array_of_mosqu
    bad_mosqu=[]
    for k in 1:length(array_of_mosqu)
        if(array_of_mosqu[k].state==dead || array_of_mosqu[k].state==lost)
            push!(bad_mosqu,(k,array_of_mosqu[k].ostate))
        end
    end
    bad_mosqu
end

function replace_lost_m()
    """ place the lost mosquitos with a new mosquito """ 
    global array_of_mosqu, model250
    nrows=model250.nrows
    ncols=model250.ncols
    for k in 1:length(array_of_mosqu)
        if(array_of_mosqu[k].state==lost) # set it to a random pos
            array_of_mosqu[k].px=Point(rand(1:nrows),rand(1:ncols))
            array_of_mosqu[k].walk=[] # reset walk
            array_of_mosqu[k].state=2 # reset the state suspected
        end
    end
end

function release_lost_m()
    """ set the mosquito to zombie """
    global array_of_mosqu, model250
    nrows=model250.nrows
    ncols=model250.ncols
    counter=0
    for k in 1:length(array_of_mosqu)
        if(array_of_mosqu[k].state==lost)
            array_of_mosqu[k].state=12 # zombie
            counter+=1
        end
    end
    counter
end # release_lost_m()
        
# visualization of the walk
# =========================

function show_way_mosquitoes(M::Vector{Any},G::Grid)
    gx=copy(G)  # make a copy
    nrows=G.nrows
    ncols=G.ncols
    znorm_Grid!(gx)   # normalize the grid
    for m in M
        for p in 1:length(m.walk)
            if(m.walk[p][1]<nrows && m.walk[p][2]<ncols &&
                m.walk[p][1]>=1 && m.walk[p][2]>=1)
                gx.grid[m.walk[p][1],m.walk[p][2]]=3.0
            end #if
        end # for
    end # for
    show(gx)
end # show_Way_Mosquitoes

function show_way_mosquitoes(M::Vector{Mosquito},G::Grid)
    gx=copy(G)  # make a copy
    nrows=G.nrows
    ncols=G.ncols
    znorm_Grid!(gx)   # normalize the grid
    for m in M
        for p in 1:length(m.walk)
            if(m.walk[p][1]<nrows && m.walk[p][2]<ncols &&
                m.walk[p][1]>=1 && m.walk[p][2]>=1)
                gx.grid[m.walk[p][1],m.walk[p][2]]=3.0
            end #if
        end # for
    end # for
    show(gx)
end # show_Way_Mosquitoes

# some function to add or remove mosquitoes
# -----------------------------------------
function gen_zombie_m(n,state)
    """ set n S in array_of_mosqu to zombie """
    global array_of_mosqu
    for k in 1:length(array_of_mosqu)
        if(array_of_mosqu[k].state==state) # S, E, I
            array_of_mosqu[k].state=10+state    # ZS,ZE,ZI
            n-=1
            if(n<=0)
                return
            end #if
        end #if
    end # for
end # gen_zombie

function add_m(n::Int64,DoY::Int64,gx::Grid,state::Int64,limit=10.0f0)
    """ adds n mosquitoes at random positions """
    global array_of_mosqu
    if(n==0) # || state==0 || state==10)
        return 
    end # if
    if(n<0)
        gen_zombie_m(-n,state)
        return
    end # if
    n1=copy(n)  # save the n
    nrows=gx.nrows
    ncols=gx.ncols
    count=1       # counts the number
    while(n1>=count)
        i=rand(1:nrows)
        j=rand(1:ncols)
        if(gx.grid[i,j]<limit)
            push!(array_of_mosqu,Mosquito(Point(i,j),DoY,state)) # new mosqito
            count+=1
        end # if
    end # while
 end #add

#function distance_m(i,j,k)
#    """ returns the manhatten distance between i,j array_of_birds[k] """
#    return (abs(array_of_mosqu[k].px.i-i)+abs(array_of_mosqu[k].px.j-j))
#end # distance

function distance_m(i,j,k)
""" returns distance between i,j array_of_birds[k] """
    return sqrt((array_of_mosqu[k].px.i-i)^2+abs(array_of_mosqu[k].px.j-j)^2)
end # distance

function min_distance_m(i,j)
    """ retruns the minimal manhatten distance """
    dmin=Inf
    k1=-999
    for k in 1:length(array_of_mosqu)
        if(array_of_mosqu[k].state==LM || array_of_mosqu[k].state==SM)
            d=distance_m(i,j,k)
            if(dmin>d) 
                dmin=d
                k1=k
            end #if
        end #if
    end#for
    k1
end # min_distance

mean_point_m=Dict()
mean_point_m[1]=(125,125)   # L
mean_point_m[2]=(125,125)   # S
mean_point_m[3]=(125,125)   # E
mean_point_m[4]=(125,125)   # I

function mean_point_of_infection_m(sel=4)
    global array_of_birds, array_of_mosqu, mean_point_m
    ax=array_of_mosqu
    end1=length(ax)
    i=[] # collects the i coord of ifected mosquitoes
    j=[] # collects the i coord of ifected mosquitoes
    for k in 1:end1
        if(ax[k].state==sel) # infected
            push!(i,ax[k].px.i)
            push!(j,ax[k].px.j)
        end # if
    end # for
    if(isempty(i))
        return mean_point_m[sel] # return to the middle of game
    end # if
    mean_point_m[sel]=(Int64(round(mean(i))),Int64(round(mean(j))))
end # mean_point_of_infection

function select_infected()
    global array_of_birds, array_of_mosqu, mean_point_m
    infected1=[]
    for i in 1:length(array_of_mosqu)
        if(array_of_mosqu[i].state==14)
            push!(infected1,i)
        end #if
    end# for
    infected1
end # function select_infected

# set a number if infected mosuitos near by i,j
function add_infected_m(i::Int64,j::Int64,dx::Int64,dy::Int64,
                        n1::Int64,DoY::Int64,state=4)
    """ set a defined number n as infected with a short distance """
    global array_of_birds, array_of_mosqu, siedlung250d, model250d
    for k in 1:n1
        di=rand(-dy:dy)
        dj=rand(-dx:dx)
        while(model250d.grid[i+di,j+dj]>10.0f0)
            di=rand(-dy:dy)
            dj=rand(-dx:dx)
        end # while
        push!(array_of_mosqu,Mosquito(Point(i+di,j+dj),DoY,state))
    end # for
    mean_point_of_infection_m()  # set the point of infection
end # set_infected_m

# changes a state from to state to
#function change_state_m(from,to,DoY,nr=1)
#    """ help function to change a state """
#    if(nr<=0)
#        return 0
#    end
#    global array_of_mosqu, mean_infection_m
#    ax=array_of_mosqu
#    end1=length(ax)
#    counter=nr
#    if(from==SM)
#        for k1 in 1:nr
#            i,j=mean_infection_m[3]
#            k2=min_distance_m(i,j)
#            if(k2==-999) 
#                return counter
#            end # if
#            ax[k2].state=to
#            ax[k2].DoY=DoY
#            counter-=1
#        end # for
#    else
#        for k1 in 1:nr
#            for k2 in 1:end1
#                if(ax[k2].state==from)
#                    ax[k2].state=to
#                    ax[k2].DoY=DoY
#                    counter-=1
#                    if(counter<=0)
#                        return 0
#                    end # if
#                end # if
#            end # for
#        end #for
#    end # else
#    counter
#end# change_state

