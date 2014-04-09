###
# Helper Functions for building the UC Models
###
using JuMP, Gurobi

    #VG 
    #There are a non-trivial number where 
    #ecomin = ecomax by hour for all hours (or some hours)
    #Some even where Ecomax = 0 for all hours (or some hours)
    #create this subclass within secondStage fcns
    #for gen-hrs where eco-min == eco-max, 
    # fix the on-off state, fix the production level (constant, non-adaptive)

    ######




function addFirstStage(m; HRS=24)
    @defVar(m, ys[1:HRS], Bin)  #stopping
    @defVar(m, zs[1:HRS], Bin)  #starting
    @defVar(m, xs[1:HRS], Bin)  #on
    
    for ix = 1:HRS
        @addConstraint(m, zs[ix] <= xs[ix])  #Start now -> On
        @addConstraint(m, ys[ix] <= 1-xs[ix])  #Stop now -> Off
    end
    
    #assume everyone starts offg
    @addConstraint(m, xs[1] <= zs[1])
    @addConstraint(m, ys[1] == 0)
    
    for ix = 2:HRS
        @addConstraint(m, zs[ix] <= 1-xs[ix-1]) #Start now -> Was off
        @addConstraint(m, ys[ix] <= xs[ix-1])  #Stop now -> Was on
    
        @addConstraint(m, xs[ix] - xs[ix-1] <= zs[ix])  #Was off, now on -> Start
        @addConstraint(m, xs[ix-1] - xs[ix] <= ys[ix])  #Was on, now off -> Start
    end
    return xs, zs, ys
end

#Neglects ending effects
#Interpret minUp <= 1 to be redundant.
function addMinUp(m, zs, xs, minUp::Int; HRS=24)
    if minUp <= 1
        return
    end
    for ix = 1:HRS
        for jx = ix+1:min(HRS, ix + minUp - 1)  #ix counts as 1 hour by itself
            @addConstraint(m, zs[ix] <= xs[jx])
        end
    end
end

function addMinDown(m, ys, xs, minDown::Int; HRS=24)
    if minDown <= 1
        return
    end
    for ix = 1:HRS
        for jx = ix+1:min(HRS, ix + minDown - 1)  #ix counts as 1 hour by itself
            @addConstraint(m, ys[ix] <= 1- xs[jx])
        end
    end
end

