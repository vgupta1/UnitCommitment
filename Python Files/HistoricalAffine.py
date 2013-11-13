""" Historical analysis Affine Model
"""
import pdb, csv, numpy, math
import generator, buildUC2

# read in the data and build a single model
gen_dict = generator.doEverything()
file_loads = csv.reader(open("../load_validate_set.csv", "rU") ) 
file_loads.next() # burn one for the header
avg_loads_by_hr = file_loads.next()
avg_loads_by_hr = [float(l) for l in avg_loads_by_hr[2:26]] #drop the date column and row indx
TMSR_REQ = 622.5 * 1e-3
T10_REQ  = 1245. * 1e-3
T30_REQ = 1883. * 1e-3

#Thin it a bit
gen_dict, avg_loads_by_hr, TMSR_REQ, T10_REQ, T30_REQ = generator.smallTestCase( 
        gen_dict, avg_loads_by_hr, TMSR_REQ, T10_REQ, T30_REQ)

print "Num Generators:\t", len(gen_dict)
for iType in ("Steam", "CT", "Diesel", "Hydro", "Nuclear", "FixedImport"):
    print iType, len( filter( lambda g: g.fuel_type == iType, gen_dict.values() ) )

#pull out hte resiudlas
eps = .1
delta = .2
file_resids = csv.reader(open("../forestResid.csv", "rU") )
file_resids.next() #burn the headers
resids = [line[2:] for line in file_resids]  #burn two for the date and the row name
resids = numpy.array(resids, dtype=float) * 1e-3

#load prediction data by day
file_preds = csv.reader(open("../forestFit.csv", "rU"))
file_preds.next() #burn a header

# for each day in the validation set,
old_objs = []
file_out_costs = csv.writer(open("affine_backtest_costs.csv", "w"))
file_out_costs.writerow(["Date", "Fixed", "Predicted", "Fixed", "Variable"])

file_out = csv.writer(open("affine_backtest.csv", "w"))
file_out.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )

#build the model for the second stage stuff
UC2obj =  buildUC2.__buildNomNoLoad(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, 
        includeIncDecs = False, useReserve=True, sparseRamps=False)

#### Begin attempt to solve a nominal run for starting values for affine
# (on_vals_init, start_vals_init, fixedCost, totCost, prod_by_hour, variable_costs) = buildUC2.buildSolveNom(
#         gen_dict, TMSR_REQ, T10_REQ, T30_REQ, avg_loads_by_hr, includeIncDecs = False, 
#         useReserve = True, sparseRamps = True)

#### end solve nominal run for starting values for affine
#build one model for the affine stuff
(m, gen_dict, prod_vars, flex_loads, on_vars, start_vars, stop_vars, 
                fixed_cost_var, reserve_vars, res_cnsts, variable_cost_vars, slack, fvecs ,fvecs_norm, gprod_sys, uncSet) = (
                buildUC2.buildAffineModel(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, resids, eps, delta, 
                        includeIncDecs = False, sparseRamps = True, method="BSApprox", omitRamping=True) 
                         )

old_objs = []
for ix, (line_load, line_fit) in enumerate(zip(file_loads, file_preds)):
    if ix > 0:
        sys.exit()

    #Update the affine model and solve
    fit_load_by_hr = [float(l) * 1e-3 for l in line_fit[2:26] ]
    (on_vals, start_vals, fixedCost1, totCost1, prod_by_hr, variable_costs), old_objs = buildUC2.addSolveAffineLoadBalanceNaive(
            m, gen_dict, prod_vars, flex_loads, on_vars, start_vars, stop_vars, fixed_cost_var, reserve_vars, 
            res_cnsts, fit_load_by_hr, variable_cost_vars, old_objs, slack, fvecs, fvecs_norm, gprod_sys, uncSet)

    #extract the period 1 solutions, and solve stage 2
#    loads = [float(l) * 1e-3 for l in line_load[2:26]]
#    (on_vals, start_vals, fixedCost2, totCost2, prod_by_hr, variable_costs), old_objs = buildUC2.updateSolveSecondStage( 
#            UC2obj, loads, old_objs, start_vals, on_vals, gen_dict )

    # dump some statistics for the run
#    file_out_costs.writerow([ line_load[1], fixedCost1, totCost1 - fixedCost1, fixedCost2, totCost2 - fixedCost2] )
#    writeHourlyGens(file_out, line_load[1], prod_by_hr)    
