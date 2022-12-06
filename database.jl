# Definition of the database using DataFrames
# ===========================================
using DataStructures: OrderedDict

# Statistical analysis of the DGL results
# =======================================

function solj2df()
    """ generates a DataFrame from solj for statistical analysis """
    SOLJ=OrderedDict()
    SOLJ["L"]=solj[1,:]
    SOLJ["SM"]=solj[2,:]
    SOLJ["EM"]=solj[3,:]
    SOLJ["IM"]=solj[4,:]
    SOLJ["SB"]=solj[5,:]
    SOLJ["EB"]=solj[6,:]
    SOLJ["IB"]=solj[7,:]
    SOLJ["RB"]=solj[8,:]
    SOLJ["DB"]=solj[9,:]
    df=DataFrame(SOLJ)
end

# will be automatically called at the end of an simulation
function save_control()
    SAVE=OrderedDict()
    SAVE["ID"]=1
    SAVE["R"]=R
    SAVE["species"]=species
    SAVE["bird"]=birds
    SAVE["start"]=start
    SAVE["steps"]=steps
    SAVE["init"]=initial_targets     # where it start the simulation
    SAVE["LM"]=Int64(round(u0mx[1])) # number of larvae at t=0
    SAVE["M"]=Int64(round(u0mx[2]))  # number suspected mosquitoes at t=0
    SAVE["IM"]=Int64(round(u0mx[4])) # number of infectious mosquitoes at t=0
    SAVE["B"]=Int64(round(u0bx[1]))  # number of suspected birds at t=0
    SAVE["IB"]=Int64(round(u0bx[3])) # number of infectious birds at t=0
    SAVE["infstepIB"]=infstepIB      # simulation step of infection bird
    SAVE["addIB"]=addIB              # number of infectious birds
    SAVE["infstepIM"]=infstepIM      # simulation step of infection mosquites
    SAVE["addIM"]=addIM              # number of infectious mosquites
    SAVE["walkn"]=walkn              # number of walks/day
    SAVE["walk_target"]=walk_targets # target of the walk
    SAVE["remark"]=remark # 12 letter limit
    SAVE["valid"]=true    # default is true can be used to clean the data
    
    if(isfile("simulations.csv")==false) # if not so create it
        df=DataFrame(SAVE)
        CSV.write("simulations.csv",df,delim=';')
    else  # read the csv and append the new SAVE to it
        df = read_simulations()
        SAVE["ID"]=1+df[end,:].ID
        push!(df,SAVE)
        CSV.write("simulations.csv",df,delim=';')
    end
    # save the control
    if(isdir("controls")==false)
        mkdir("controls")
    end
    dst=string("controls/","control_",SAVE["ID"],".jl") # save the control.jl
    cp("control.jl",dst)
    _valid_=true   # reset the _valid_ after activation
    # save the dfx
    # ------------
    dfx=solj2df()
    if(isdir("results")==false)
        mkdir("results")
    end
    dst=string("results/","result_",SAVE["ID"],".csv")
    CSV.write(dst,dfx,delim=';')
    df,dfx
end # function save_control

# shows the stored simulations as DataFrame
function read_simulations(path="./")
    filename=string(path,"simulations.csv")
    if(isfile(filename))
        df = CSV.File(filename,delim=';',
                      types=[Int64,String,String,String,Int64,Int64,
                             String,Int64,Int64,Int64,Int64,Int64,
                             Float64,Int64,Float64,Int64,Int64,
                             String,String,Bool]) |> DataFrame
    else
        println("error: could not open: ", filename)
    end
    df
end

# set or unset the valid flag (true, false)
function switch_valid(ID)
    for i in 1:size(df)[1]
        if(df[i,:].ID==ID)
            df[i,:].valid= ! df[i,:].valid
            return df[i,:].valid
        end #if
    end #for
end # switch_valid

# set the remark use line number to address
function set_remark(ID,_s_)
    for i in 1:size(df)[1]
         if(df[i,:].ID==ID)
             df[row,:].remark=_s_
             return df[row,:].remark
         end # if
    end #for
end # set_remark
    
# order the DataFrame according to R, species and start
function sort_df(order=[:R, :species, :start])
    sort!(df,order)
end

# clean removes all valid=false from simulations.csv and
# remove the control_i.jl
# please handle with care

function clean_controls()
    for i in 1:size(df)[1]
        if(df[i,:].valid==false)
            filename=string("controls/","control_",df[i,:].ID,".jl")
            if(isfile(filename))
                rm(filename)
            end
        end
    end
    new_df=df[(df.valid).==true,:]
    CSV.write("simulations.csv",new_df,delim=';')
    read_simulations()
end

# store a adapted df (valid, remark)
function save_df(df)
    CSV.write("simulations.csv",df,delim=';')
end

# activate an control_?.jl
# it replaces the actual control.jl with control_?.jl
function activate(ID)
    global _valid_
    filename=string("controls/","control_",ID,".jl")
    cp(filename,"control.jl"; force=true)
    true
end
    
    
