""" Historical analysis Affine Model
"""
import pdb, csv, numpy, math
import generator, buildUC2

# read in the data and build a single model
gen_dict_ = generator.doEverything()

#Load up the averages
file_avg_loads = csv.reader(open("../load_means.csv", "rU") )
file_avg_loads.next() #burn one for the header
avg_loads_by_hr = []
for line in file_avg_loads:
    avg_loads_by_hr.append( float( line[1] ) * 1e-3 )
avg_loads_by_hr = numpy.array(avg_loads_by_hr)

#Thin it a bit
gen_dict, load_ratio = generator.smallTestCase( gen_dict_, filt_percent = .1 )
#buildUC2.calcResReqs(gen_dict, avg_loads_by_hr * load_ratio)
TMSR_REQ = (0.683 * .5)/2.
T10_REQ = (0.683)/2.
T30_REQ = (0.683 + .5 * 0.272)/2.

#These represent the averages over the whole set
# TMSR_REQ = 622.5 * 1e-3
# T10_REQ  = 1245. * 1e-3
# T30_REQ = 1883. * 1e-3

TMSR_REQ = 0.
T10_REQ = 0.
T30_REQ = 0.


print "\n Num Generators:\t", len(gen_dict)
for iType in generator.FUEL_TYPES:
    print iType, len( filter( lambda g: g.fuel_type == iType, gen_dict.values() ) )

print "\n Capacity Requirements:"
print TMSR_REQ, T10_REQ, T30_REQ

print "\n Load Ratio: \t", load_ratio

#pull out hte residuals
eps = .1
delta = .2
file_resids = csv.reader(open("../forestResid.csv", "rU") )
file_resids.next() #burn the headers
resids = [line[2:] for line in file_resids]  #burn two for the date and the row name
resids = numpy.array(resids, dtype=float) * 1e-3 * load_ratio * 2.5
print "Mean Resids:\t", resids.shape
print numpy.mean(resids, axis=0)

#load prediction data by day
file_preds = csv.reader(open("../forestFit.csv", "rU"))
file_preds.next() #burn a header

file_loads = csv.reader(open("../load_validate.csv", "rU") ) 
file_loads.next() # burn one for the header

#output files
file_out_costs = csv.writer(open("affine_backtest_costs.csv", "w"))
file_out_costs.writerow(["Date", "Fixed", "Predicted", "Fixed", "Variable"])
file_out = csv.writer(open("affine_backtest.csv", "w"))
file_out.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )
file_out_planned = csv.writer(open("affine_backtest_planned.csv", "w"))
file_out_planned.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )

#build one model for the nominal second stage stuff
UC2obj =  buildUC2.__buildNomNoLoad(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, False, True, False)

#solve the nominal for a good starting solution

#build one model for the affine stuff
(m, gen_dict, prod_vars, flex_loads, on_vars, start_vars, stop_vars, 
                fixed_cost_var, reserve_vars, variable_cost_vars, slack, fvecs ,fvecs_norm, gprod_sys, uncSet) = (
                buildUC2.buildAffineModel(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, resids, eps, delta, 5, 
                                                                omitRamping=True) )

# for each day in the validation set,
old_objs_affine = []
old_objs_nom = []
for ix, (line_load, line_fit) in enumerate(zip(file_loads, file_preds)):
    if ix > 2:
        sys.exit()

    fit_load_by_hr = [float(l) * 1e-3  * load_ratio for l in line_fit[2:26] ]

    #solve the nominal for a mip start
#     on_vals_init, start_vals_init, fixedCost_init, totCost_init, prod_by_hour_init, variable_costs_init = (
#         buildUC2.buildSolveNom(gen_dict, TMSR_REQ, T10_REQ, T30_REQ, fit_load_by_hr) )
#     buildUC2.addMIPStart(m, on_vals_init, start_vals_init, on_vars, start_vars)

    #Update the affine model and solve
    (on_vals, start_vals, fixedCost1, totCost1, prod_by_hr, variable_costs), old_objs_affine = buildUC2.addSolveAffineLoadBalanceNaive(
            m, gen_dict, prod_vars, flex_loads, on_vars, start_vars, stop_vars, fixed_cost_var, reserve_vars, 
            fit_load_by_hr, variable_cost_vars, old_objs_affine, slack, fvecs, fvecs_norm, gprod_sys, uncSet)

    buildUC2.writeHourlySchedCap(file_out_planned, line_load[1], on_vals, gen_dict)

    #extract the period 1 solutions, and solve stage 2
    loads = [float(l) * 1e-3 * load_ratio for l in line_load[2:26]]

    ##VG HACK
    resids = [ (l - lp) * 2.5 for l, lp in zip(loads, fit_load_by_hr) ]
    loads= [l + r for l,r in zip(loads, resids)]

    (on_vals, start_vals, fixedCost2, totCost2, prod_by_hr, variable_costs), old_objs_nom = buildUC2.updateSolveSecondStage( 
           UC2obj, loads, old_objs_nom, start_vals, on_vals, gen_dict )

    # dump some statistics for the run
    file_out_costs.writerow([ line_load[1], fixedCost1, totCost1 - fixedCost1, fixedCost2, totCost2 - fixedCost2] )
    buildUC2.writeHourlyGens(file_out, line_load[1], prod_by_hr)    
