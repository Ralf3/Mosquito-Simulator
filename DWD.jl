"""
DWD is a modul to load DWD station as a table
The load of the DWD data is organized as selection for each investigation
area. The modeler selects one region and the DWD data will be loaded 

Additionally the modul provides some functions

"""
# describe the metadata for the DWD station
#--------------------------------------------
Trebnitz=Dict()
Trebnitz["wind"]="Trebnitz_wind_stunde_02878"
Trebnitz["wetter"]="Trebnitz"
Trebnitz["lat"]=51.3909
Trebnitz["lon"]=11.8786
Trebnitz["alt"]=117.38
#--------------------------------------------
Dueben=Dict()
Dueben["wind"]="Bad_Dueben_wind_stunde_05629"
Dueben["wetter"]="Dueben"
Dueben["lat"]=51.4348
Dueben["lon"]=12.2396
Dueben["alt"]=105.00
#--------------------------------------------
Radebeul=Dict()
Radebeul["wind"]="Radebeul_Dresden_wind_stunde_01048"
Radebeul["wetter"]="Radebeul"
Radebeul["lat"]=51.1278
Radebeul["lon"]=13.7543
Radebeul["alt"]=228
#--------------------------------------------
Bitterfeld=Dict()
Bitterfeld["wind"]="Bitterfeld_wind_stunde_02932"
Bitterfeld["wetter"]="Bitterfeld"
Bitterfeld["lat"]=51.4348
Bitterfeld["lon"]=12.2396
Bitterfeld["alt"]=131
#--------------------------------------------
Magdeburg=Dict()
Magdeburg["wind"]="Magdeburg_wind_stunde_03126"
Magdeburg["wetter"]="Magdeburg"
Magdeburg["lat"]=52.1029
Magdeburg["lon"]=11.5827
Magdeburg["alt"]=79

# load the climate data
# ---------------------
function load_DWD(selector,from=20200101,to=20201231)
    file=string("../DWD/",selector["wetter"],".txt")
    table=CSV.File(file) |> DataFrame
    tmk=table[((table[!, "MESS_DATUM"] .>= from) .&
               (table[!, "MESS_DATUM"] .< to)), :" TMK"];
    LT=LinearInterpolation(1:length(tmk),tmk);
    return LT
end

function load_Wind(selector,from=2020010100,to=2020123100)
    file=string("../DWD/",selector["wind"],".txt")
    table=CSV.File(file) |> DataFrame
    table1=table[((table[!, "MESS_DATUM"] .>= from) .&
        (table[!, "MESS_DATUM"] .< to)), :]
    POW=[]
    DIR=[]
    step=Int64(size(table1)[1]/24-1) # daily step
    for i in 0:step
        push!(POW,table[18+24*i,4]) # 4: "  FF"
        push!(DIR,table[18+24*i,5]) # 5: "  DD"
    end
    return POW,DIR
end

#biting_rate(T)=0.344/(1+1.231*exp(-0.184*(T-20)))
#mL(T)=(0.0025*T^2-0.094*T+1.0257)       # mortality of larve
#mM(T)=0.1*mL(T)                         # mortality of mosquito 

function daylength(dayOfYear,Station)
    """
        the Station is Dict [Trebnitz, Dueben]
    """
    lat=Station["lat"]
    lon=Station["lon"]
    latInRad = deg2rad(lat)
    declinationOfEarth =
        23.45*sin(deg2rad(360.0*(283.0+dayOfYear)/365.0))
    if -tan(latInRad) * tan(deg2rad(declinationOfEarth)) <= -1.0
        return 24.0
    end
    if -tan(latInRad) * tan(deg2rad(declinationOfEarth)) >= 1.0
        return 0.0
    end
    hourAngle =
        rad2deg(acos(-tan(latInRad)*tan(deg2rad(declinationOfEarth))))
    return 2.0*hourAngle/15.0
end

# 
δM(dayOfYear,Station)=1.0-
    1.0/(1.0+1775.7*exp(1.559*(daylength(dayOfYear,Station)-18.177)))

# connection to the wind drift model of mosquitoes
# ================================================

direction=[-1 -1; -1 0; 0 -1; 0 0; 0 1; 1 0; -1 1; 1 -1; 1 1];
function  genWindPropability(dayOfYear,pow=POW,dir=DIR)
    """ 
    calculates the probs for wind using DWD wind data at 17:00 (POW,DIR)
    for one selected dayOfYear
    return the probs
    """
    probs=ones(size(direction)[1]) # define a new probs
    limit=10.0 # limit for wind driven flight
    if(POW[dayOfYear]<0 || POW[dayOfYear]>limit || DIR[dayOfYear]<0)
        return(probs) # wind power or direction is out of range
    end
    # definition of direction 0:360 degrees
    # north:1   270:90  
    # south:-1  90:270
    # east:1     0:180   
    # west:-1    180:0
    # calc windi,windj using the triangular
    α=DIR[dayOfYear]*π/180       # angle to 
    windj=sin(α)*POW[dayOfYear]  # j: x-direction
    windi=cos(α)*POW[dayOfYear]  # i: y-direction
    # println(POW[dayOfYear],"  ", DIR[dayOfYear]," i: ",windi," j: ",windj)
    wi=windi/limit #  normalize the wind
    wj=windj/limit #  normalize the wind
    # add the wi to the probs according to direction
    maxdir=0.0
    ik=1
    for i in 1:9
        val=direction[i,1]*wi+direction[i,2]*wj
        if(maxdir<val)
            maxdir=val
            ik=i
        end
    end
    # println(ik,"  ", maxdir)
    probs[ik]=sqrt(POW[dayOfYear])
    probs/=sum(probs)         # norm 
    return(cumsum(probs))     # cdf
end
