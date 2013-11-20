""" Historical analysis of Nominal Problem
"""
import pdb, csv, numpy
import generator, buildUC2

# read in the data and build a single model
gen_dict = generator.doEverything()
load_ratio = 1.0

#thin the example
#gen_dict, load_ratio = generator.smallTestCase( gen_dict_, filt_percent = .1 )
TMSR_REQ = (1.245 * .5)
T10_REQ = (1.245)
T30_REQ = (1.245 + .5 * 1.237)

#Load up the averages
file_avg_loads = csv.reader(open("../load_means.csv", "rU") )
file_avg_loads.next() #burn one for the header
avg_loads_by_hr = []
for line in file_avg_loads:
    avg_loads_by_hr.append( float( line[1] ) * 1e-3 * load_ratio)
avg_loads_by_hr = numpy.array(avg_loads_by_hr)

print "Num Generators:\t", len(gen_dict)
for iType in ("Steam", "CT", "Diesel", "Hydro", "Nuclear", "FixedImport"):
    print iType, len( filter( lambda g: g.fuel_type == iType, gen_dict.values() ) )

buildUC2.calcResReqs(gen_dict, avg_loads_by_hr * load_ratio)

print "\n Capacity Requirements:"
print TMSR_REQ, T10_REQ, T30_REQ

print "\n Load Ratio:\t", load_ratio

# Build and solve a nominal model
on_vals, start_vals, fixed_cost_init, tot_cost_init, prod_by_hr, variable_costs = buildUC2.buildSolveNom(
        gen_dict, TMSR_REQ, T10_REQ, T30_REQ, avg_loads_by_hr, 
        includeIncDecs = False, useReserve = True, sparseRamps = False)

# dump some statistics for the run
file_out_costs = csv.writer(open("nominalReserve_backtest_costs.csv", "w"))
file_out_costs.writerow(["Date", "Fixed", "Predicted", "Fixed", "Variable"])

file_out_costs2 = csv.writer(open("hindsight_costs.csv", "w"))
file_out_costs2.writerow(["Date", "Fixed", "Predicted", "Fixed", "Variable"])


file_out = csv.writer(open("nominalReserve_backtest.csv", "w"))
file_out.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )

file_out2 = csv.writer(open("hindsight_backtest.csv", "w"))
file_out2.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )


#write out just once the planned production
file_out_planned = csv.writer(open("nominalReserve_backtest_planned.csv", "w"))
file_out_planned.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )
buildUC2.writeHourlySchedCap(file_out_planned, "Planned", on_vals, gen_dict)


#build the model for the second stage stuff
UC2obj =  buildUC2.__buildNomNoLoad(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, 
        includeIncDecs = False, useReserve=True, sparseRamps=False)


file_loads = csv.reader(open("../load_validate.csv", "rU") ) 
file_loads.next() # burn one for the header


# for each day in the validation set,
old_objs = []
for ix, line in enumerate(file_loads):
    if "NA" in line:
        raise ValueError()
    loads = [float(l) * 1e-3 * load_ratio for l in line[2:26] ]

    (on_vals, start_vals, fixedCost, tot_cost, prod_by_hr, variable_costs), old_objs = buildUC2.updateSolveSecondStage( 
            UC2obj, loads, old_objs, start_vals, on_vals, gen_dict )

    # dump the realized generation scheme including slacks....
    dt = line[1]
    file_out_costs.writerow([dt, fixed_cost_init, tot_cost_init - fixed_cost_init, fixedCost, tot_cost - fixedCost] )
    buildUC2.writeHourlyGens(file_out, dt, prod_by_hr)


    #solve it with perfect hindsight
    on_vals, start_vals, fixedCostHind, totcostHind, prod_by_hour, variable_costs = buildUC2.buildSolveNom(
            gen_dict, TMSR_REQ, T10_REQ, T30_REQ, loads)
    file_out_costs2.writerow([dt, fixed_cost_init, tot_cost_init - fixed_cost_init, fixedCostHind, totcostHind - fixedCostHind] )
    buildUC2.writeHourlyGens(file_out2, dt, prod_by_hr)
    
   
    
    
