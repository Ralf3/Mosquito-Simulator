#using Pkg
#Pkg.activate("/home/ralf/ABM")
#using Distributions
#using SpecialFunctions
#using DataStructures
#using Interpolations
#using Plots

DEBUG=false   # plot the bB of all birds

# ========================================================================
# data Amsel (Tomialojc 1994)
# ========================================================================

function amsel_bB()
    reprate_urban = 2.0 # urban habitat, reproduction rate per 1 bird per year
    gx=Gamma(48.004389717013744, 2.833490870352461) # from fit
    function bB_am(d)
        reprate_urban*pdf(gx,d)  # [105,...,185]
    end
    bB_am
end

# ========================================================================
# data Elster (Gorski & Kotlarz 1987)
# ========================================================================

function elster_bB()
    reprate_urban = 0.658 # urban habitat,
    gx=Gamma(75.43790459820515, 1.7039179532441344)
    function bB_el(d)
        reprate_urban*pdf(gx,d)  # [105,...,185]
    end
    bB_el
end

# ========================================================================
# data Nebelkraehe (Zduniak & Kuczynski 2003)
# ========================================================================

function kraehe_bB()
    reprate_rural = 0.919 # rural habitat
    gx=Gamma(290.52988577683277, 0.4201151274803995)
    function bB_kr(d)
        reprate_rural*pdf(gx,d)  # [105,...,185]
    end
    bB_kr
end

# ========================================================================
# data Habicht (Looft 2017)
# ========================================================================

function habicht_bB()
    reprate = 0.667 # reproduction rate per 1 bird per year
    gx=Gamma(477.31448641987004, 0.2946422201740774)
    function bB_ha(d)
        reprate*pdf(gx,d)  # [105,...,185]
    end
    bB_ha
end

# ========================================================================
# data Sperling (Wegrzynowicz 2017)
# ========================================================================

function sperling_bB()
    spx = [100, 105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165,
           170, 175, 180, 185]

    spy = [0.09, 0.17, 0.1, 0.065, 0.025, 0.015, 0.025, 0.04, 0.1, 0.1, 0.08,
           0.025, 0.025, 0.015, 0.015, 0.055, 0.04, 0.015]

    # add border elements
    dx=spx[2]-spx[1]     # 5
    dy=minimum(spy)/10   # 0.0015
    spx=vcat([spx[1]-dx],spx,[spx[end]+dx])
    spy=vcat([dy],spy,[dy])
    # Interpolation
    spi = CubicSplineInterpolation(95:5:190, spy)
    res=[]  # store the results of the interpolation
    for i in 95:1:190
        push!(res,spi(i))
    end
    scaler=sum(res) # scaler is used to sum(res)==1.0
    function bB_sp(d)
        reprate =  3.11
        if(d<95 || d>190)
            return 0.0
        end
        return reprate*spi(d)/scaler
    end
    return bB_sp
end

# ========================================================================
# data Gruenfink 
# ========================================================================

function gruenfink_bB()
    gfx = [105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165,
           170, 175, 180, 185, 190, 195, 200, 205, 210]

    gfy = [0.0037, 0.00741, 0.03333, 0.2037, 0.09926, 0.09259, 0.0963, 0.08148,
           0.04444, 0.04815, 0.02223, 0.05185, 0.05556, 0.07778, 0.06296,
           0.05556, 0.02593, 0.03333, 0.01481, 0.01852, 0.00741, 0.0037]

    # add border elements
    dx=gfx[2]-gfx[1]     # 5
    dy=minimum(gfy)/10   # 0.00037
    gfx=vcat([gfx[1]-dx],gfx,[gfx[end]+dx])
    gfy=vcat([dy],gfy,[dy])
    # Interpolation
    gfi = CubicSplineInterpolation(100:5:215, gfy)
    res=[]  # store the results of the interpolation
    for i in 100:1:215
        push!(res,gfi(i))
    end
    scaler=sum(res)  # scaler is used to sum(res)==1.0
    function bB_gf(d)
        reprate =  2.15
        if(d<100 || d>215)
            return 0.0
        end
        return reprate*gfi(d)/scaler
    end
    return bB_gf
end

# ========================================================================
# Tests
# ========================================================================

if(DEBUG==true)
    bB_am=amsel_bB()
    plot(bB_am,80:300)
    bB_el=elster_bB()
    plot!(bB_el,80:300)
    bB_kr=kraehe_bB()
    plot!(bB_kr,80:300)
    bB_ha=habicht_bB()
    plot!(bB_ha,80:300)
    bB_sp=sperling_bB()
    plot!(bB_sp,80:300)
    bB_gf=gruenfink_bB()
    plot!(bB_gf,80:300)
end

# Habicht: NB=1000*0.42 BP=1000*(1-0.42)=580
# Habicht: NB1=(1-0.43)*NB=239 BP1=(1-0.43)*BP=330
# Habicht: Juvenil=1.1*(BP+NB)=1100 Juvenil1=Juvenil*(1-0.62)=418
# Habicht: Gesamt1=239+330+418=987

# Elster: NB=1000*0.14=140 BP=1000(1-0.14)=860
# Elster: NB1=NB*(1-0.513)=68 BP1=BP*(1-0.513)=419
# Elster: Juvenil=1.2*(BP+NB)=1200 Juvenil1=Juvenil*(1-0.7)=360
# Elster: Gesamt1=68+419+360=847

