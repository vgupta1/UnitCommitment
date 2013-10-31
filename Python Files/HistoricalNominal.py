""" Historical analysis of Nominal Problem
"""
import csv, pdb
import gurobipy as grb
import generator, buildUC

def writeHourlyGens(file_out, label):
    file_out.writerow([label, "Total"] + 
        [ prod_by_type["TOTAL", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "Nuclear"] + 
        [ prod_by_type["Nuclear", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "Hydro"] + 
        [ prod_by_type["Hydro", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "Steam"] + 
        [ prod_by_type["Steam", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "CT"] + 
        [ prod_by_type["CT", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "Diesel"] + 
        [ prod_by_type["Diesel", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "FixedImport"] + 
        [ prod_by_type["FixedImport", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "Other"] + 
        [ prod_by_type["Other", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "INC"] + 
        [ prod_by_type["INC", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "DEC"] + 
        [ prod_by_type["DEC", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "FLEX"] + 
        [ prod_by_type["FLEX", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "FLEX"] + 
        [ prod_by_type["FLEX", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "LOAD"] + 
        [ prod_by_type["LOAD", hr] for hr in xrange(buildUC.HORIZON_LENGTH) ] )
    file_out.writerow([label, "Slack"] + 
        [ slack[hr].x for hr in xrange(buildUC.HORIZON_LENGTH) ] )

# read in the data and build a single model
gen_dict = generator.doEverything(True)
m = grb.Model("UCNominal")
on_vars, start_vars, stop_vars, cost_var = buildUC.genStage1Vars( m, gen_dict )
prod_vars, reserve_vars = buildUC.genStage2VarsNom(m, gen_dict, True )
variable_cost_vars = buildUC.addPiecewiseCosts(m, gen_dict, prod_vars )
buildUC.startStopConstraints(m, gen_dict, on_vars, start_vars, stop_vars)
buildUC.ecoMinMaxConsts(m, gen_dict, prod_vars, on_vars, reserve_vars)
buildUC.minUpConstraints(m, gen_dict, on_vars)
buildUC.minDownConstraints(m, gen_dict, on_vars)
flex_loads = buildUC.genFlexibleDemandVarsNom( m, gen_dict )
buildUC.rampingConsts(m, gen_dict, prod_vars, start_vars, stop_vars)
reserve_cnsts = buildUC.reserveRequirements(m, gen_dict, reserve_vars)
buildUC.reserveCapacity(m, gen_dict, reserve_vars, reserve_cnsts)

file_loads = csv.reader(open("ISO-data/load_validate_set.csv", "rU") ) 
avg_loads_by_hr = file_loads.next()
avg_loads_by_hr = [float(l) for l in avg_loads_by_hr]

# solve the nominal model, export the planned generation scheme
slack, balance_cnsts = buildUC.nominalLoadBalance(m, gen_dict, prod_vars, avg_loads_by_hr, 
                                                    flex_loads, {}, {})

m.params.mipgap = 1e-2
m.optimize()

file_out = csv.writer(open("nominalReserve_backtest.csv", "w"))
file_out.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )
file_out_costs = csv.writer(open("nominalReserve_backtest_costs.csv", "w"))

file_out_costs.writerow(["Date", "Fixed", "Variable"])
file_out_costs.writerow(["Planned", cost_var.x, m.objVal - cost_var.x ] )

prod_by_type = buildUC.computeGenByHr(gen_dict, prod_vars, 
                                                    flex_loads, {}, {})
writeHourlyGens(file_out, "Planned")

# for each day in the validation set,
for ix, line in enumerate(file_loads):
    if "NA" in line:
        continue
    loads = [float(l) for l in line]

    if ix == 0:
        slack, balance_cnsts, m = buildUC.solveSecondStage(m, balance_cnsts, gen_dict, prod_vars, loads, 
                                              flex_loads, {}, {}, on_vars, start_vars, stop_vars, reserve_cnsts, reserve_vars)
    else:  #only remove the reserve stuff once.
        slack, balance_cnsts, m = buildUC.solveSecondStage(m, balance_cnsts, gen_dict, prod_vars, loads, 
                                              flex_loads, {}, {}, on_vars, start_vars, stop_vars, {}, {})

    m.optimize()
    # dump the realized generation scheme including slacks....
    file_out_costs.writerow([str(ix), cost_var.x, m.objVal - cost_var.x] )
    prod_by_type = buildUC.computeGenByHr(gen_dict, prod_vars, 
                                                        flex_loads, {}, {})
    writeHourlyGens(file_out, str(ix))
    
    
    
    
    
