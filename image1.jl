"""
Revised version according to large images: Katrin Kuhls 21.7.2021 
"""

# load an OSM image
#im1=load("images/data.png")
#imshow(im1)
# generate three maps for each color
function transfer_Image_Array(im)
    r=zeros(Float32,size(im))
    b=zeros(Float32,size(im))
    g=zeros(Float32,size(im))
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            r[i,j]=im[i,j].r
            g[i,j]=im[i,j].g
            b[i,j]=im[i,j].b
        end
    end
    r,g,b
end
# find four corners of an image which are not white
function find_Corner_Image(im)
    flag_ul=0
    flag_ur=0
    flag_ll=0
    flag_lr=0
    ul=(1,1)
    ur=(1,1)
    ll=(1,1)
    lr=(1,1)
    for i in 1:1:size(im)[1]
        for j in 1:1:size(im)[2]
            if(flag_ul==0 && flag_ur==0 && flag_ll==0 && flag_lr==0 &&
                im[i,j].r!=1.0 && im[i,j].b!=1.0 && im[i,j].g!=1.0)
                ul=(i,j) # find the first none white column
                flag_ul=1
                continue
            end
            if(flag_ul==1 && flag_ur==0 && flag_ll==0 && flag_lr==0 &&
                im[i,j].r==1.0 && im[i,j].b==1.0 && im[i,j].g==1.0)
                ur=(i,j-1) # find the first white column
                flag_ur=1
                break
            end
        end
    end
    
    for i in size(im)[1]:-1:1
        for j in 1:1:size(im)[2]
            if(flag_ul==1 && flag_ur==1 && flag_ll==0 && flag_lr==0 &&
                im[i,j].r!=1.0 && im[i,j].b!=1.0 && im[i,j].g!=1.0)
                ll=(i,j) # find the first none white row 
                flag_ll=1
                continue
            end
            if(flag_ul==1 && flag_ur==1 && flag_ll==1 && flag_lr==0 &&
                im[i,j].r==1.0 && im[i,j].b==1.0 && im[i,j].g==1.0)
                lr=(i,j)
                flag_lr=1
                return ul,ur,ll,lr
            end
        end
    end
end

function cut_Corner_Image(ul,ur,ll,lr,im,dl=0,dr=0,du=0,dd=0)
    """ cut from ul,ur,ll,lr add delta: dl, dr, du, dd """
     im[(ul[1]+dl):(ll[1]-dr),(ul[2]+du):(ur[2]-dd)]
end

#im2=cut_Corner_Image(ul,ur,ll,lr,im1,10,59,15,55)
#r,g,b=transfer_Image_Array(im2)

function select_fluss(im)
    river=zeros(Float32,size(im)) # define a grid with river=1 other=0
    R=152/255
    G=201/255
    B=216/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                river[i,j]=1.0
            end
         end
     end
     river
end

function select_see(im)
    lake=zeros(Float32,size(im))
    R=170/255
    G=211/255
    B=223/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                lake[i,j]=1.0
            end
        end
    end
    lake
end

function select_acker(im)
    acker=zeros(Float32,size(im))
    R=244/255
    G=236/255
    B=221/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                acker[i,j]=1.0
            end
        end
    end
    acker
end

function select_brache1(im)
    brache1=zeros(Float32,size(im))
    R=242/255
    G=239/255
    B=233/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                brache1[i,j]=1.0
            end
        end
    end
    brache1
end

function select_brache2(im)
    brache2=zeros(Float32,size(im))
    R=227/255
    G=230/255
    B=189/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                brache2[i,j]=1.0
            end
        end
    end
    brache2
end

function select_bauland(im)
    bau=zeros(Float32,size(im))
    R=199/255
    G=199/255
    B=181/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                bau[i,j]=1.0
            end
        end
    end
    bau
end

function select_buschland(im)
    busch=zeros(Float32,size(im))
    R=217/255
    G=227/255  # 216
    B=197/255  # 171
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                busch[i,j]=1.0
            end
        end
    end
    busch
end

function select_wald(im)
    wald=zeros(Float32,size(im))
    R=199/255
    G=224/255
    B=189/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                wald[i,j]=1.0
            end
        end
    end
    wald
end

function select_wiese(im)
    wiese=zeros(Float32,size(im))
    R=205/255
    G=235/255
    B=176/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                wiese[i,j]=1.0
            end
        end
    end
    wiese
end

function select_garten(im)
    garten=zeros(Float32,size(im))
    R=215/255
    G=233/255
    B=208/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(B-im[i,j].g) <= 0.01 &&
               abs(G-im[i,j].b) <= 0.01)
                garten[i,j]=1.0
            end
        end
    end
    garten
end

function select_friedhof(im)
    fried=zeros(Float32,size(im))
    R=199/255
    G=218/255
    B=200/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                fried[i,j]=1.0
            end
        end
    end
    fried
end

function select_park(im)
    park=zeros(Float32,size(im))
    R=211/255
    G=250/255
    B=214/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                park[i,j]=1.0
            end
        end
    end
    park
end

function select_sport(im)
    sport=zeros(Float32,size(im))
    R=223/255
    G=252/255
    B=226/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                sport[i,j]=1.0
            end
        end
    end
    sport
end

function select_stadion(im)
    stadion=zeros(Float32,size(im))
    R=170/255
    G=224/255
    B=203/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                stadion[i,j]=1.0
            end
        end
    end
    stadion
end

function select_siedlung(im)
    sied=zeros(Float32,size(im))
    R=221/255 # 233
    G=221/255 # 223
    B=221/255 # 223
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                sied[i,j]=1.0
            end
        end
    end
    sied
end

function select_industrie(im)
    indus=zeros(Float32,size(im))
    R=238/255
    G=225/255
    B=236/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                indus[i,j]=1.0
            end
        end
    end
    indus
end

function select_gewerbe1(im)
    gewerbe1=zeros(Float32,size(im))
    R=221/255
    G=221/255
    B=221/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                gewerbe1[i,j]=1.0
            end
        end
    end
    gewerbe1
end

function select_gewerbe2(im)
    gewerbe2=zeros(Float32,size(im))
    R=221/255
    G=221/255
    B=221/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                gewerbe2[i,j]=1.0
            end
        end
    end
    gewerbe2
end

function select_deponie(im)
    deponie=zeros(Float32,size(im))
    R=207/255
    G=207/255
    B=184/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                deponie[i,j]=1.0
            end
        end
    end
    deponie
end

function select_bauernhof(im)
    bauer=zeros(Float32,size(im))
    R=238/255
    G=226/255
    B=206/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                bauer[i,j]=1.0
            end
        end
    end
    bauer
end

function select_bergbau(im)
    berg=zeros(Float32,size(im))
    R=212/255
    G=212/255
    B=212/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                berg[i,j]=1.0
            end
        end
    end
    berg
end

function select_gebaeude(im)
    gebau=zeros(Float32,size(im))
    R=209/255
    G=198/255
    B=189/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                gebau[i,j]=1.0
            end
        end
    end
    gebau
end
               
function select_autobahn(im)
    auto=zeros(Float32,size(im))
    R=234/255
    G=128/255
    B=86/255
    R1=244/255
    G1=211/255
    B1=82/255
    for i in 1:size(im)[1]
        for j in 1:size(im)[2]
            if(abs(R-im[i,j].r) <= 0.01 &&
               abs(G-im[i,j].g) <= 0.01 &&
               abs(B-im[i,j].b) <= 0.01)
                auto[i,j]=1.0
            end
        end
    end
    auto
end

# read map
#im1=load("images/data.png")
#ul,ur,ll,lr=find_Corner_Image(im1)
#im2=cut_Corner_Image(ul,ur,ll,lr,im1,10,59,15,55)
# read all features

function gen_maps(im2)
    maps=Dict()
    # fill the maps
    maps["fluss"]=select_fluss(im2)
    maps["see"]=select_see(im2)
    maps["acker"]=select_acker(im2)
    maps["brache1"]=select_brache1(im2)
    maps["brache2"]=select_brache2(im2)
    maps["bauland"]=select_bauland(im2)
    maps["buschland"]=select_buschland(im2)
    maps["wald"]=select_wald(im2)
    maps["wiese"]=select_wiese(im2)
    maps["garten"]=select_garten(im2)
    maps["friedhof"]=select_friedhof(im2)
    maps["park"]=select_park(im2)
    maps["sport"]=select_sport(im2)
    maps["stadion"]=select_stadion(im2)
    maps["siedlung"]=select_siedlung(im2)
    maps["industrie"]=select_industrie(im2)
    maps["gewerbe1"]=select_gewerbe1(im2)
    maps["gewerbe2"]=select_gewerbe2(im2)
    maps["deponie"]=select_deponie(im2)
    maps["bauernhof"]=select_bauernhof(im2)
    maps["bergbau"]=select_bergbau(im2)
    maps["gebaeude"]=select_gebaeude(im2)
    maps["autobahn"]=select_autobahn(im2)
    return maps
end

# imgg = imfilter(maps["siedlung"], Kernel.gaussian(5))
# imshow(imgg.>0.1)
# npzwrite("/home/ralf/Julia/Pipiens/Geodaten/Magdeburg/Magdeburg_Siedlung.npz",
#          imgg.>0.1)


#selection=[]
#for key in keys(maps)
#    if(sum(maps[key])>100.0f0)
#        push!(selection,key)
#    end
#end
