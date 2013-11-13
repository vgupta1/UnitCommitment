""" Historical analysis of Nominal Problem
"""
import pdb, csv
import generator, buildUC2

# read in the data and build a single model
gen_dict = generator.doEverything()
file_loads = csv.reader(open("../load_validate_set.csv", "rU") ) 
file_loads.next() # burn one for the header
avg_loads_by_hr = file_loads.next()
avg_loads_by_hr = [float(l) * 1e-3 for l in avg_loads_by_hr[2:26]] #drop the date column and row indx
TMSR_REQ = 622.5 * 1e-3
T10_REQ  = 1245. * 1e-3
T30_REQ = 1883. * 1e-3

#thin the instance
gen_dict, avg_loads_by_hr, TMSR_REQ, T10_REQ, T30_REQ = generator.smallTestCase( 
        gen_dict, avg_loads_by_hr, TMSR_REQ, T10_REQ, T30_REQ)

print "Num Generators:\t", len(gen_dict)
for iType in ("Steam", "CT", "Diesel", "Hydro", "Nuclear", "FixedImport"):
    print iType, len( filter( lambda g: g.fuel_type == iType, gen_dict.values() ) )

# Build and solve a nominal model
on_vals, start_vals, fixed_cost, tot_cost, prod_by_hr, variable_costs = buildUC2.buildSolveNom(
        gen_dict, TMSR_REQ, T10_REQ, T30_REQ, avg_loads_by_hr, 
        includeIncDecs = False, useReserve = True, sparseRamps = True)

# dump some statistics for the run
file_out_costs = csv.writer(open("nominalReserve_backtest_costs.csv", "w"))
file_out_costs.writerow(["Date", "Fixed", "Variable"])
file_out_costs.writerow(["Planned", fixed_cost, tot_cost- fixed_cost ] )

file_out = csv.writer(open("nominalReserve_backtest.csv", "w"))
file_out.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )
buildUC2.writeHourlyGens(file_out, "Planned", prod_by_hr)

#build the model for the second stage stuff
UC2obj =  buildUC2.__buildNomNoLoad(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, 
        includeIncDecs = False, useReserve=True, sparseRamps=False)

# for each day in the validation set,
old_objs = []
for ix, line in enumerate(file_loads):
    if "NA" in line:
        continue
    loads = [float(l) for l in line[2:26]]

    if ix > 5:
        sys.exit()

    (on_vals, start_vals, fixedCost, tot_cost, prod_by_hr, variable_costs), old_objs = buildUC2.updateSolveSecondStage( 
            UC2obj, loads, old_objs, start_vals, on_vals, gen_dict )

    # dump the realized generation scheme including slacks....
    dt = line[1]
    file_out_costs.writerow([dt, fixedCost, tot_cost - fixedCost] )
    writeHourlyGens(file_out, dt, prod_by_hr)
    
    
    
    
    
