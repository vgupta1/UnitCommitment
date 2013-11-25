""" Historical analysis of Nominal Problem
"""
import pdb, csv, numpy, sys
import generator, buildNom, summarize

tag1 = sys.argv[1]
print tag1

gen_dict = generator.doEverything()
load_ratio = 1.0
TMSR_REQ = (1.245 * .5)
T10_REQ = (1.245)
T30_REQ = (1.245 + .5 * 1.237)
TMSR_REQ = T10_REQ = T30_REQ = 0.

#thin the example
#gen_dict, load_ratio = generator.smallTestCase( gen_dict_, filt_percent = .1 )

#Load up the averages
file_avg_loads = csv.reader(open("../load_means.csv", "rU") )
file_avg_loads.next() #burn one for the header
avg_loads_by_hr = []
for line in file_avg_loads:
    avg_loads_by_hr.append( float( line[1] ) * 1e-3 * load_ratio)
avg_loads_by_hr = numpy.array(avg_loads_by_hr)

# Set up the output files
file_out_costs = csv.writer(open(tag1 + "_costs.csv", "w"))
file_out_costs.writerow(["Date", "Fixed", "Predicted", "Fixed", "Variable"])

file_out_costs2 = csv.writer(open(tag1 + "_hind_costs.csv", "w"))
file_out_costs2.writerow(["Date", "Fixed", "Predicted", "Fixed", "Variable"])

file_out = csv.writer(open(tag1 + "_sched.csv", "w"))
file_out.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )

file_out2 = csv.writer(open(tag1 + "_hind_sched.csv", "w"))
file_out2.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )

print "Num Generators:\t", len(gen_dict)
for iType in ("Steam", "CT", "Diesel", "Hydro", "Nuclear", "FixedImport"):
    print iType, len( filter( lambda g: g.fuel_type == iType, gen_dict.values() ) )

print "\n Capacity Requirements:"
print TMSR_REQ, T10_REQ, T30_REQ

print "\n Load Ratio:\t", load_ratio

# Build and solve a nominal model
on_vals, start_vals, fixed_cost_init, tot_cost_init, prod_by_hr, variable_costs = buildNom.buildSolveNom(
        gen_dict, TMSR_REQ, T10_REQ, T30_REQ, avg_loads_by_hr)

summarize.calcOpLimits(on_vals, gen_dict)
pdb.set_trace()

# Try solving a similar problem with a lower profile
vFixedCosts = []
scaling_grid = numpy.linspace(.1, .4, num=10)
for scaling in scaling_grid:
    on_vals, start_vals, fixedCostHind, totcostHind, prod_by_hour, variable_costs = buildNom.buildSolveNom(
            gen_dict, TMSR_REQ, T10_REQ, T30_REQ, scaling * avg_loads_by_hr)
    vFixedCosts.append( fixedCostHind )

for s, val in zip(scaling_grid, vFixedCosts):
    print s, val, fixed_cost_init

sys.exit()

#write out just once the planned production
file_out_planned = csv.writer(open(tag1 + "_plan.csv", "w"))
file_out_planned.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )
summarize.writeHourlySchedCap(file_out_planned, "Planned", on_vals, gen_dict)

#build the model for the second stage stuff
UC2obj =  buildNom.__buildNomNoLoad(gen_dict, 0., 0., 0., False, False)

file_loads = csv.reader(open("../load_validate.csv", "rU") ) 
file_loads.next() # burn one for the header
# for each day in the validation set,
old_objs = []
for ix, line in enumerate(file_loads):
    if "NA" in line:
        raise ValueError()
    loads = [float(l) * 1e-3 * load_ratio for l in line[2:26] ]

    on_vals, start_vals, fixedCost, tot_cost, prod_by_hr, variable_costs = buildNom.updateSolveSecondStage(
            UC2obj, loads, gen_dict, on_vals, start_vals )

    # dump the realized generation scheme including slacks....
    dt = line[1]
    file_out_costs.writerow([dt, fixed_cost_init, tot_cost_init - fixed_cost_init, fixedCost, tot_cost - fixedCost] )
    summarize.writeHourlyGens(file_out, dt, prod_by_hr)

    #solve it with perfect hindsight
    on_vals, start_vals, fixedCostHind, totcostHind, prod_by_hour, variable_costs = buildNom.buildSolveNom(
            gen_dict, 0., 0., 0., loads)
    file_out_costs2.writerow([dt, fixed_cost_init, tot_cost_init - fixed_cost_init, fixedCostHind, totcostHind - fixedCostHind] )
    summarize.writeHourlyGens(file_out2, dt, prod_by_hr)
    
   
    
    
