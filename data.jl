# load the spatial data
# =====================

path=string("./",R,"/")
if(R=="Bitterfeld")
    siedlung800d=load_grid(string(path,"Bitterfeld_Siedlung800d.jld2"))
    siedlung800dx=load_grid(string(path,"Bitterfeld_Siedlung800dx.jld2"))
    siedlung800dy=load_grid(string(path,"Bitterfeld_Siedlung800dy.jld2"))
    river800d=load_grid(string(path,"Bitterfeld_River800d.jld2"))
    river800dx=load_grid(string(path,"Bitterfeld_River800dx.jld2"))
    river800dy=load_grid(string(path,"Bitterfeld_River800dy.jld2"))
    model_=load_Image_Grid(string(path,string(species,"_",_year_,".npy")))
    #blood_=load_Image_Grid(string(path,"Blood.npz"))
    #sugar_=load_Image_Grid(string(path,"Sugar.npz"))
end

if(R=="Magdeburg")
    siedlung800d=load_grid(string(path,"Magdeburg_Siedlung800d.jld2"))
    siedlung800dx=load_grid(string(path,"Magdeburg_Siedlung800dx.jld2"))
    siedlung800dy=load_grid(string(path,"Magdeburg_Siedlung800dy.jld2"))
    model_=load_Image_Grid(string(path,string(species,"_",_year_,".npy")))
    river800d=load_grid(string(path,"Magdeburg_River800d.jld2"))
    river800dx=load_grid(string(path,"Magdeburg_River800dx.jld2"))
    river800dy=load_grid(string(path,"Magdeburg_River800dy.jld2"))
    #blood_=load_Image_Grid(string(path,"Blood.npz"))
    #sugar_=load_Image_Grid(string(path,"Sugar.npz"))
end

if(R=="Trebnitz")
    siedlung800d=load_grid(string(path,"Trebnitz_Siedlung800d.jld2"))
    siedlung800dx=load_grid(string(path,"Trebnitz_Siedlung800dx.jld2"))
    siedlung800dy=load_grid(string(path,"Trebnitz_Siedlung800dy.jld2"))
    river800d=load_grid(string(path,"Trebnitz_River800d.jld2"))
    river800dx=load_grid(string(path,"Trebnitz_River800dx.jld2"))
    river800dy=load_grid(string(path,"Trebnitz_River800dy.jld2"))
    model_=load_Image_Grid(string(path,string(species,"_",_year_,".npy")))
    #blood_=load_Image_Grid(string(path,"Blood.npz"))
    #sugar_=load_Image_Grid(string(path,"Sugar.npz"))
end

if(R=="Dueben")
    path=string("../Geodaten/",R,"/")
    siedlung800d=load_grid(string(path,"Dueben_Siedlung800d.jld2"))
    siedlung800dx=load_grid(string(path,"Dueben_Siedlung800dx.jld2"))
    siedlung800dy=load_grid(string(path,"Dueben_Siedlung800dy.jld2"))
    river800d=load_grid(string(path,"Dueben_River800d.jld2"))
    river800dx=load_grid(string(path,"Dueben_River800dx.jld2"))
    river800dy=load_grid(string(path,"Dueben_River800dy.jld2"))
    model_=load_Image_Grid(string(path,string(species,"_",_year_,".npy")))
    #blood_=load_Image_Grid(string(path,"Blood.npz"))
    #sugar_=load_Image_Grid(string(path,"Sugar.npz"))
end

if(R=="Radebeul")
    siedlung800d=load_grid(string(path,"Radebeul_Siedlung800d.jld2"))
    siedlung800dx=load_grid(string(path,"Radebeul_Siedlung800dx.jld2"))
    siedlung800dy=load_grid(string(path,"Radebeul_Siedlung800dy.jld2"))
    river800d=load_grid(string(path,"Radebeul_River800d.jld2"))
    river800dx=load_grid(string(path,"Radebeul_River800dx.jld2"))
    river800dy=load_grid(string(path,"Radebeul_River800dy.jld2"))
    model_=load_Image_Grid(string(path,string(species,"_",_year_,".npy")))
    #blood_=load_Image_Grid(string(path,"Blood.npz"))
    #sugar_=load_Image_Grid(string(path,"Sugar.npz"))
end

# scale all maps to 250*250
siedlung250d=scale_Grid(siedlung800d,250,250)    # distance to the siedlung
siedlung250d.grid.=siedlung250d.grid./3.2        # scale it from 800*800
siedlung250dx=scale_Grid(siedlung800dx,250,250)  # dx(distanceto the siedlung)
siedlung250dy=scale_Grid(siedlung800dy,250,250)  # dy(distanceto the siedlung)
river250d=scale_Grid(river800d,250,250)          # distance to the river
river250d.grid.=river250d.grid/3.2               # scale it from 800*800 
river250dx=scale_Grid(river800dx,250,250)        # dx(distance to the river)
river250dy=scale_Grid(river800dy,250,250)        # dy(distance to the river)
model250=enlarge_Grid(model_,10)                 # model from OSM 25 ==> 250

# vexans dict
VD=Dict("Dueben" => 0.65f0, "Bitterfeld" => 0.8f0,
        "Magdeburg" => 0.95f0, "Radebeul" => 0.95f0, "Trebnitz" => 0.65f0)
# pipiens dict
PD=Dict("Dueben" => 0.65f0, "Bitterfeld" => 0.9f0,
        "Magdeburg" => 0.9f0, "Radebeul" => 0.95f0, "Trebnitz" => 0.75f0)
# select dict
SD=Dict("vexans" => VD, "pipiens" => PD)

norm_Grid!(model250)                             # norm before
model250d=distance_Grid_Sel(model250,SD[species][R])  # distance to best habitat
model250dy=diffy_Grid(model250d)                 # dx(model)
model250dx=diffx_Grid(model250d)                 # dy(model)

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
#.jl POW,DIR=load_Wind(DWD_data[R], _from_, _to_);
