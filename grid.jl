# ideas from statistics
# =====================
using Formatting

function mediansorted(x::AbstractVector{T}, i::Integer, l::Integer)::T where T
    len = l - i + 1
    len > zero(len) || throw(ArgumentError("Array slice cannot be empty."))
    mid = i + len ÷ 2
    return isodd(len) ? x[mid] : (x[mid-1] + x[mid]) / 2
end

function fivenum(x::AbstractVector{T}) where T<:AbstractFloat
    r   = zeros(T,5)
    xs  = sort(x)
    mid = length(xs) ÷ 2
    lowerend = isodd(length(xs)) ? mid : mid - 1
    r[1] = xs[1]
    r[2] = mediansorted(xs, 1, lowerend)
    r[3] = mediansorted(xs, 1, length(xs))
    r[4] = mediansorted(xs, mid, length(xs))
    r[end] = xs[end]
    return r
end

# ideas from samt2
# ================

# simple points as mutable struct
# ===============================
mutable struct Point
       i::Int64
       j::Int64
end
create_Points(N)=Vector{Point}(undef,N)
function Base.copy(p::Point)
    px=Point(p.i,p.j)
    return px
end
    
# define a simple rect as mutable structure
# =========================================

mutable struct Rect
    lli::Int64 # lower left i 
    uri::Int64 # upper right i
    llj::Int64 # lower left j
    urj::Int64 # upper right j
end

# simple mutable structure of a grid
# ==================================
mutable struct Grid
    nrows::Int64
    ncols::Int64
    xllcorner::Float64
    yllcorner::Float64
    nodata::Float32
    csize::Float32
    grid::Matrix{Float32}    
end
# constructor of a Grid
function Grid(nrows::Int64,ncols::Int64)
    return Grid(nrows,ncols
                ,0.0,0.0,-9999.0f0,1.0f0,Float32.(zeros(nrows,ncols)))
end
# constructor of an empty Grid
function Grid()
    return Grid(1,1,0.0,0.0,-9999.0f0,1.0f0,Float32.(zeros(1,1)))
end

# Base function for visualization of Grids
Base.show(X::Grid) = Plots.heatmap(X.grid)
Plots.histogram(G::Grid)=histogram(reshape(G.grid,length(G.grid)))
function Plots.plot(G::Grid)
    rs=reshape(G.grid,G.nrows*G.ncols)
    rs1=reshape(Gray.(Float32.(rs)),(G.nrows,G.ncols))
    plot(rs1)
end

# some statistics
# ---------------
fivenum(X::Grid) = fivenum(reshape(X.grid,X.nrows*X.ncols))
# sample from grid
function sample_Grid(N::Int64,G::Grid,limit::Float32)
    """ 
    uses a 2D grid as virtually habitat model
    the limit cut all values below
    returns N samples produced fron the 2D grid
    """
    rows=zeros(G.nrows) # store the sum of all rows
    cols=zeros(G.ncols) # store the sum of all cols
    gx=copy(G)
    gx.grid[gx.grid.<limit].=0.0 # set all values = 0.0 below limit
    count=1             # counter for while loop
    res=create_Points(N)# allocate memory for points
    while count<=N
        i=rand(1:gx.nrows)
        j=rand(1:gx.ncols)
        if(gx.grid[i,j]>=limit)
            res[count]=Point(i,j)
            count+=1
        end
    end
    return(res,gx) 
end

# help function to visualize the sampled points  
function vis_Point_Grid(pp::Vector{Point},G::Grid)
    gx=copy(G)
    # gx.grid[:,:].=0.0
    for i in 1:length(pp)
        gx.grid[pp[i].i,pp[i].j]=1.0f0
    end
    show(gx)
end

# Base function for handling Grids
function Base.copy(G::Grid)
    g1=Grid()
    g1.nrows=G.nrows
    g1.ncols=G.ncols
    g1.xllcorner=G.xllcorner
    g1.yllcorner=G.yllcorner
    g1.nodata=G.nodata
    g1.csize=G.csize
    g1.grid=copy(G.grid)
    return g1
end

# hadle ASCII grids elementary
function load_ascii(filename)
    g1=Grid()               # define a new empty grid
    io=open(filename,"r")   # open the ascii grid
    line=readline(io)
    print(line)
    m = match(r"(\w+)\s+(\d+)",line)
    g1.ncols=parse(Int64,m[2])
    line=readline(io)
    print(line)
    m = match(r"(\w+)\s+(\d+)",line)
    g1.nrows=parse(Int64,m[2])
    line=readline(io)
    print(line)
    m = match(r"(\w+)\s+(\d+[\.\d*])",line)
    g1.xllcorner=parse(Float64,m[2])
    line=readline(io)
    print(line)
    m = match(r"(\w+)\s+(\d+[\.\d*])",line)
    g1.yllcorner=parse(Float64,m[2])
    line=readline(io)
    print(line)
    m = match(r"(\w+)\s+(\d*)",line)
    g1.csize=parse(Float32,m[2])
    line=readline(io)
    print(line)
    m = match(r"(\w+)\s+([+-]\d*)",line)
    g1.nodata=parse(Float32,m[2])
    g1.grid=zeros((g1.nrows,g1.ncols))
    for i in 1:g1.nrows
        line=readline(io)
        s=split(line[2:end]," ")
        for j in 1:g1.ncols
            g1.grid[i,j]=parse(Float32,s[j])
        end
    end
    close(io)  # do not forget it!
    g1
end
    
function save_ascii(filename::String,G::Grid) # ASCII format
    io=open(filename,"w")
    write(io,sprintf1("ncols     %d\n",G.ncols))
    write(io,sprintf1("nrows     %d\n",G.nrows))
    write(io,sprintf1("xllcorner     %f\n",G.xllcorner))
    write(io,sprintf1("yllcorner     %f\n",G.yllcorner))
    write(io,sprintf1("cellsize      %f\n",G.csize))
    write(io,sprintf1("NODATA_value     %d\n",G.nodata))
    for i in 1:G.nrows
        for j in 1:(G.ncols)
            write(io,sprintf1("%f ",G.grid[i,j]))
        end
        write(io,"\n") # end of line
    end
    close(io) # do not forget it!
end

function save_grid(filename::String,G::Grid) # jld2 format
    ncols=G.ncols
    nrows=G.nrows
    xllcorner=G.xllcorner
    yllcorner=G.yllcorner
    nodata=G.nodata
    csize=G.csize
    grid=G.grid
    jldsave(filename;ncols,nrows,xllcorner,yllcorner,nodata,csize,grid)
end

function load_grid(filename::String)  # jld2 format
    g1=Grid()
    xxx=load(filename)
    g1.nrows=xxx["nrows"]
    g1.ncols=xxx["ncols"]
    g1.xllcorner=xxx["xllcorner"]
    g1.yllcorner=xxx["yllcorner"]
    g1.nodata=xxx["nodata"]
    g1.csize=xxx["csize"]
    g1.grid=xxx["grid"]
    g1
end
    
function gen_random_Grid(nrows::Int64,ncols::Int64)
    grid=Grid(nrows,ncols)
    grid.grid=Float32.(rand(nrows,ncols))
    return grid
end

function load_xyz(filename,csize)
    # load a DEM from https://data.geobasis-bb.de/geobasis/daten/dgm/xyz/
    # do not forget to add "X" "Y" "Z" in top of the data
    data=CSV.read(filename,DataFrame)
    g1=Grid()
    g1.ncols=Int64(maximum(data[:,1])-minimum(data[:,1]))+1
    g1.nrows=Int64(maximum(data[:,2])-minimum(data[:,2]))+1
    g1.xllcorner=minimum(data[:,1])
    g1.yllcorner=minimum(data[:,2])
    g1.nodata=-9999
    g1.csize=csize
    g1.grid=zeros(Float32,g1.nrows,g1.ncols)
    # read in the grid (Matrix{Float32})
    for i in 1:g1.nrows
        for j in 1:g1.ncols
            g1.grid[i,j]=data[g1.nrows*(i-1)+j,3]
        end
    end
    return top_down_Grid(g1) # from top to down
end
                  
function load_Image_Grid(filename,rev=true)
    """
    load a grid from a numpy array stored as npy ore npz
    for all the formats from controls
    """
    data=nothing
    try
        data=npzread(filename)

    catch
        println("filename: ",filename, " not found")
        return noting
    end
    out=Grid(size(data)[1],size(data)[2])
    if(rev==true)
        for i in 1:size(data)[1]
            out.grid[size(data)[1]+1-i,:]=Float32.(data[i,:])
        end
    else
        for i in 1:size(data)[1]
            out.grid[i,:]=Float32.(data[i,:])
        end
    end
    out
end

function load_npz_Grid(filename,rev=false)
    """
    load a grid from a numpy array stored as npz
    """
    data=nothing
    try
        data=npzread(filename)["arr_0"]

    catch
        println("filename: ",filename, " not found")
        return noting
    end
    out=Grid(size(data)[1],size(data)[2])
    if(rev==true)
        for i in 1:size(data)[1]
            out.grid[size(data)[1]+1-i,:]=Float32.(data[i,:])
        end
    else
        for i in 1:size(data)[1]
            out.grid[i,:]=Float32.(data[i,:])
        end
    end
    out
end

function top_down_Grid(G::Grid)
    """ changes the top to down and returns a new grid """
    out=copy(G)
    for i in 1:size(G.grid)[1]
        out.grid[size(G.grid)[1]+1-i,:]=G.grid[i,:]
    end
    out
end
        
# some normalization of Grids
function norm_Grid!(G::Grid)
    max=maximum(G.grid)
    min=minimum(G.grid)
    G.grid=(G.grid.-min)./(max-min)
end

function znorm_Grid!(G::Grid)
    m=mean(G.grid)
    s=std(G.grid)
    G.grid=(G.grid.-m)./s
end

function invert_Grid!(G::Grid)
    max=maximum(G.grid)
    min=minimum(G.grid)
    for i in 1:G.nrows
        for j in 1:G.ncols
            G.grid[i,j]=max-G.grid[i,j]+min
        end
    end
end

# Enhanced functions of a Grid
function scale_Grid(G::Grid,nrows::Int64,ncols::Int64)
    """ 
        uses a grid G to scale it to the new nrows,ncols
        the new csize=G.csize*(G.nrows/nrows)
        result: out::Grid  with a a scaled grid
    """
    out=Grid()
    out.nrows=nrows
    out.ncols=ncols
    out.xllcorner=G.xllcorner
    out.yllcorner=G.yllcorner
    out.nodata=G.nodata
    out.csize=G.csize*(G.nrows/nrows)
    out.grid=imresize(G.grid,nrows,ncols)
    return(out)
end
    
function enlarge_Grid(G::Grid,factor::Int64)
    nrows=G.nrows
    ncols=G.ncols
    Y= Float32.(collect(0.0:nrows/(nrows-1):nrows)*factor)
    X= Float32.(collect(0.0:ncols/(ncols-1):ncols)*factor)
    spl=Spline2D(X,Y,G.grid)
    # new grid with factor size
    out=Grid(factor*G.nrows,factor*G.ncols)
    out.csize=G.csize/factor
    out.nodata=G.nodata
    for i in 1:out.nrows
        for j in 1:out.nrows
            out.grid[i,j]=spl(i,j)
        end
    end
    out
end

function varpart_Grid(G::Grid,nr=5000)
    """
    divides the grid into a set of rectangels according to the variance
    returns: a set of three dictionary: 
    var_dict : containing the variance of recangles
    mean_dict: containing the mean of the rectangles
    rect_dict: containing the rectangles itself
    """
    # defines the dicts for storage
    var_dict=Dict()  # stores the variance of each rectangle
    mean_dict=Dict() # stores the mean value of each rectangle
    rect_dict=Dict() # stores the rectangle itself
    # fill the whoole grid
    rect_dict[1]=Rect(1,G.nrows,1,G.ncols) # lli,uri,llj,urj
    r0=rect_dict[1]
    var_dict[1]=var(G.grid[r0.lli:r0.uri,r0.llj:r0.urj])
    mean_dict[1]=mean(G.grid[r0.lli:r0.uri,r0.llj:r0.urj])
    last_key=1 # last key in the dicts
    for l in 1:nr
        # find the rectangle with the maximum variance
        maxvar=0.0
        max_key=1 # key of var_dict[key] => max
        for k in 1:last_key
            if(var_dict[k]>maxvar)
                maxvar=var_dict[k]
                max_key=k
            end
        end
        # select the new rectangel
        r=rect_dict[max_key]
        dy=r.uri-r.lli # length in i direction
        dx=r.urj-r.llj # length in j direction
        if(dx<=3 && dy<=3)
            var_dict[max_key]=0.0
            continue
        end
        if(dy>3 && dy>dx)
            r1=Rect(r.lli,r.lli+Int64(round(dy/2)),r.llj,r.urj)
            r2=Rect(r.lli+Int64(round(dy/2)),r.uri,r.llj,r.urj)
        end
        if(dx>3 && dx>=dy)
            r1=Rect(r.lli,r.uri,r.llj,r.llj+Int64(round(dx/2)))
            r2=Rect(r.lli,r.uri,r.llj+Int64(round(dx/2)),r.urj)
        end
        var_dict[max_key]=var(G.grid[r1.lli:r1.uri,r1.llj:r1.urj])
        mean_dict[max_key]=mean(G.grid[r1.lli:r1.uri,r1.llj:r1.urj])
        rect_dict[max_key]=r1
        last_key+=1
        var_dict[last_key]=var(G.grid[r2.lli:r2.uri,r2.llj:r2.urj])
        mean_dict[last_key]=mean(G.grid[r2.lli:r2.uri,r2.llj:r2.urj])
        rect_dict[last_key]=r2
    end
    rect_dict,var_dict,mean_dict
end

function build_from_varpart(G::Grid,RD::Dict{Any, Any},MD::Dict{Any, Any})
    """  
    builds from the output of varpart: rect_dict and mean_dict a new Grid 
    returns
    the new grid (out), RMSE(G,out)
    """
    out=Grid(G.nrows,G.ncols)
    for k in keys(RD)
        r=RD[k]
        val=MD[k]
        out.grid[r.lli:r.uri,r.llj:r.urj].=val
    end
    RMSE=sqrt(sum((out.grid.-G.grid).^2)/length(out.grid))
    out,RMSE
end

# function with points
function random_Points_Grid(N::Int64,G::Grid)
    """
    generates N Points randomly located at Grid.grid
    returns:
    Vector{Points}
    """
    points=create_Points(N)
    for k in 1:N
        i=rand(1:G.nrows)
        j=rand(1:G.ncols)
        points[k]=Point(i,j)
    end
    points
end

function random_Points_Grid(N::Int64,sel::Float32,G::Grid)
    """
    generates N Points randomly located at Grid.grid
    returns:
    Vector{Points}
    """
    k=1
    points=create_Points(N)
    while true
        i=rand(1:G.nrows)
        j=rand(1:G.ncols)
        if(G.grid[i,j]==sel)
            points[k]=Point(i,j)
            k+=1
        end
        if(k>N)
            return points
        end
    end
end

function add_Points_Grid(points::Vector{Point},G::Grid,value=1.0)
    for p in points 
        G.grid[p.i,p.j]+=value
    end
end

function set_Points_Grid(points::Vector{Point},G::Grid,value=1.0)
    for p in points 
        G.grid[p.i,p.j]=value
    end
end

function sample_Point_Grid(P::Vector{Point},G::Grid)
    """ 
    uses Points to sample from the from a given Grid
    returns:
    Vector{Float32} of the values of the Grid
    """
    res=zeros(length(P))
    for i in 1:length(P)
        res[i]=G.grid[P[i].i,P[i].j]
    end
    res
end
    
function distance_Points_Grid(P::Vector{Point},G::Grid)
    g1=Grid(G.nrows,G.ncols)
    maxdist=Inf
    for i in 1:G.nrows
        for j in 1:G.ncols
            maxdist=Inf
            for k in 1:length(P)
                dist=sqrt((i-P[k].i)*(i-P[k].i)+(j-P[k].j)*(j-P[k].j))
                if(dist<maxdist)
                    maxdist=dist
                end
            end
            g1.grid[i,j]=maxdist
        end
    end
    g1
end

function distance_Grid_Sel(G::Grid,sel::Float32)
    """ for each grid the shortest distance to a selected value """
    g1=Grid(G.nrows,G.ncols)
    g1.grid[:,:].=Float32(sqrt(G.nrows*G.nrows+G.ncols+G.ncols))
    for i in 1:G.nrows
        # println(i)
        for j in 1:G.ncols
            if(G.grid[i,j]>=sel)
                for i1 in 1:G.nrows
                    for j1 in 1:G.ncols
                        d=sqrt((i-i1)*(i-i1)+(j-j1)*(j-j1))
                        if(g1.grid[i1,j1]>d)
                            g1.grid[i1,j1]=d
                        end
                    end
                end
            end
        end
    end
    return(g1)
end

# numerical differentiation
function diffx_Grid(G::Grid)
    gx=copy(G)
    dx=backdiffx(G.grid)
    gx.grid=dx
    return(gx)
end

function diffy_Grid(G::Grid)
    gx=copy(G)
    dy=backdiffy(G.grid)
    gx.grid=dy
    return(gx)
end

# test it
function test1()
    gx=Grid(400,400)
    points=random_Points_Grid(10,gx)
    add_Points_Grid(points,gx)
    distance_Points_Grid(points,gx)
    RD,VD,MD=varpart_Grid(gx)
    out,RMSE=build_from_varpart(gx,RD,MD)
    return gx,points,out
end

function save_River(dir,file,res=800)
    """
    calculates the distance to a river (coded with 1.0f0)
    and the diffx,diffy.
    all the files are saved.
    param:
    dir:directory name without /
    filename: with extension .npz
    res: resolution of the image
    result:
    
    """
    s=string(dir,"/",file)
    river=load_Image_Grid(s)
    river800=scale_Grid(river,res,res)
    river800d=distance_Grid_Sel(river800,1.0f0)
    river800dx= diffx_Grid(river800d)
    river800dy= diffy_Grid(river800d)
    filename=split(file,".")[1]
    s1=string(dir,"/",filename,res,"d",".jld2")
    s2=string(dir,"/",filename,res,"dx",".jld2")
    s3=string(dir,"/",filename,res,"dy",".jld2")
    save_grid(s1, river800d)
    save_grid(s2, river800dx)
    save_grid(s3, river800dy)
end
