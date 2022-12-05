# data type of a bird
# ===================
mutable struct Bird
    px::Point    # position in the the region
    DoY::Int64   # time of arival
    state::Int64 # S=1,E=2,I=3,R=4,lost=6,dead=-999, zombie=[11,12,13,14]
end
SB=1
EB=2
IB=3
RB=4
DB=5

function Bird(p0,DoY,state)
    return Bird(p0,DoY,state)
end

function Base.copy(b::Bird)
    px=copy(b.px)
    return Bird(px,DoY,state)
end

function describe(b::Bird)
    println("point: ", b.px, " DoY: ", b.DoY, " state: ", b.state)
end

# some function to handle birds
# -----------------------------

array_of_birds=[] # global array to store all birds
array_of_infected_birds=[] # global array for infected birds over time

function show_map_b(state)
    """ search the complete array """
    global array_of_birds, siedlung250d
    gx=siedlung250d        # second name for the grid
    g1=copy(siedlung250d)  # copy of the grid
    z=round(maximum(g1.grid))+1
    counter=0
    end1=length(array_of_birds)
    for i in 1:end1
        if(array_of_birds[i].state==state)
            g1.grid[array_of_birds[i].px.i,array_of_birds[i].px.j]=z
            counter+=1
        end # if
    end # for
    println(counter)
    show(g1)  # show the birds
end # show_map

function show_inf_b()
    """ search the complete array """
    global array_of_birds, siedlung250d
    gx=siedlung250d        # second name for the grid
    g1=copy(siedlung250d)  # copy of the grid
    z=round(maximum(g1.grid))+1
    counter=0
    end1=length(array_of_infected_birds)
    for i in 1:end1
        g1.grid[array_of_infected_birds[i].px.i,
                array_of_infected_birds[i].px.j]=z
        counter+=1
    end # for
    println(counter)
    show(g1)  # show the birds
end # show_map

# Some statistics about birds
# ---------------------------

function count_of_infectious_b(sel::Int64)
    """ found all infected birds over all simu; sel=IB """
    global BRD
    nr=0
    k=length(BRD)    # number of stored bird_array
    res=zeros(k)     # store the selected state of the bird
    for i in 1:k
        stop=size(BRD[i])[1]
        for j in 1:stop
            if(BRD[i][j].state==sel)
                res[i]+=1
            end
        end
    end
    res
end

function gen_hist_b(t=0)
    """ t>0 selects a copy from BRD[t] """
    array_of_b=array_of_birds
    if(t!=0)
        array_of_b=BRD[t]
    end
    Hist=Dict()
    Hist[1]=0  # S solj=5
    Hist[2]=0  # E solj=6
    Hist[3]=0  # I solj=7
    Hist[4]=0  # R solj=8
    Hist[5]=0  # D solj=9
    Hist[11]=0 # ZS
    Hist[12]=0 # ZE
    Hist[13]=0 # ZI
    Hist[14]=0 # ZR
    Hist[-1]=0 # Z
    for i in 1:length(array_of_b)
        Hist[array_of_b[i].state]+=1
    end
    Hist
end

# calc function of birds
# ----------------------
   
function gen_zombie_b(n,state)
    """ set n S in array_of_birds to zombie """
    for k in 1:length(array_of_birds)
        if(array_of_birds[k].state==state) 
            array_of_birds[k].state=state+10  # zombie
            n-=1
            if(n<=0)
                return
            end
        end
    end # for
end # gen_zombie_b

function add_dead_b(n::Int64,DoY::Int64)
    """
    avoid that a mosquito.state==dead will die again
    """
    counter=0
    for k in 1:length(array_of_birds)
        if(array_of_birds[k].state!=DB)
            array_of_birds[k].state=DB
            counter+=1
            if(counter==n)
                return(counter)
            end # if
        end # if
    end# for
    counter
end # add_dead_b

function add_recovered_b(n::Int64,DoY::Int64)
    """ changes the infected birds to recoverd birds
        if n>infected_birds then add new mosquitoes not far from
        center of infection
    """
    if(n<=0)
        return(0)
    end
    counter=0
    for k in 1:length(array_of_birds)
        if(array_of_birds[k].state==IB)
            array_of_birds[k].state=RB
            counter+=1
            # println(n," ", counter)
            if(counter==n)
                return(counter)
            end # if
        end # if
    end# for
    counter
end # add_recovered_b

function add_infected_b(n::Int64,DoY::Int64)
    """ changes the exposed birds to infected birds
        if n>infected_birds then add new mosquitoes not far from
        center of infection
    """
    if(n<=0)
        return(0)
    end
    counter=0
    for k in 1:length(array_of_birds)
        if(array_of_birds[k].state==EB)
            array_of_birds[k].state=IB
            counter+=1
            # println(n," ", counter)
            if(counter==n)
                return(counter)
            end # if
        end # if
    end# for
    counter
end # add_infected_b

function add_exposed_b(n::Int64,DoY::Int64)
    """ changes the suceptible birds to exposed birds
        if n>infected_birds then add new mosquitoes not far from
        center of infection
    """
    if(n<=0)
        return(0)
    end
    counter=0
    for k in 1:length(array_of_birds)
        if(array_of_birds[k].state==SB)
            array_of_birds[k].state=EB
            counter+=1
            # println(n," ", counter)
            if(counter==n)
                return(counter)
            end # if
        end # if
    end# for
    counter
end # add_infected_b

function add_b(n::Int64,DoY::Int64,gx::Grid,state::Int64,limit=3.0f0)
    """ add n birds at random position destance < limit
        DoY and state=1,2,3,4  are possible
        returns a vector of new birds
    """
    if(state==IB && n>0)
        add_infected_b(n,DoY)
    end
    if(state==EB && n>0)
        add_exposed_b(n,DoY)
    end
    if(n==0)
        return
    end
    if(n<=0)
        gen_zombie_b(-n,state)
        return
    end # if
    if(n>0 && state==RB)
        add_recovered_b(n,DoY)
        return
    end # if recovered_b
    if(n>0 && state==DB)
        add_dead_b(n,DoY)
        return
    end # if dead_b
    n1=copy(n)  # save the n
    global array_of_birds,siedlung250d
    nrows=siedlung250d.nrows
    ncols=siedlung250d.ncols
    count=1
    while(n1>=count)
        i=rand(1:nrows)
        j=rand(1:ncols)
        if(gx.grid[i,j]<limit)
            push!(array_of_birds,Bird(Point(i,j),DoY,state)) # new bird
            if(state==IB)
                push!(array_of_infected_birds,Bird(Point(i,j),DoY,state))
            end
            count+=1
        end # if
    end # while
end #add

mean_point_b=Dict()
mean_point_b[1]=(125,125)   # save a new mean point of infection  S
mean_point_b[2]=(125,125)   # save a new mean point of infection  E
mean_point_b[3]=(125,125)   # save a new mean point of infection  I
mean_point_b[4]=(125,125)   # save a new mean point of infection  R
mean_point_b[5]=(125,125)   # save a new mean point of infection  D

function mean_point_of_infection_b(k,sel) # k=0 sel=13 
    global array_of_birds, BRD
    global mean_point_b
    array_b=array_of_birds
    if(k>0)
        array_b=BRD[k]
    end
    end1=length(array_b)
    i=[] # collects the i coord of ifected mosquitoes
    j=[] # collects the i coord of ifected mosquitoes
    for k in 1:end1
        if(array_b[k].state==sel) # infected
            push!(i,array_b[k].px.i)
            push!(j,array_b[k].px.j)
        end # if
    end # for
    if(isempty(i))
        return mean_point_b[sel] # return to the middle of game
    end # if
    mean_point_b[sel]=(Int64(round(mean(i))),Int64(round(mean(j))))
end # mean_point_of_infection

function distance_b(i,j,k)
    """ returns the  distance between i,j array_of_birds[k] """
    return sqrt((array_of_birds[k].px.i-i)^2+(array_array_of_birds[k].px.j-j)^2)
    end # distance 

function min_distance_b(i,j)
    """ retruns the minimal distance """
    dmin=Inf
    k1=-999
    for k in 1:length(array_of_birds)
        if(array_of_birds[k].state==1)
            d=distance_b(i,j,k)
            if(dmin>d) 
                dmin=d
                k1=k
            end #if
        end #if
    end#for
    k1
end # min_distance_b
    
