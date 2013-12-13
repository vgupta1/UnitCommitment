""" Historical analysis Affine Model
"""
import pdb, csv, numpy, sys
import generator, buildAff, buildNom, summarize
import sparseAffineCG as cg
import gurobipy as grb

#VG use the skits version of bootstrapping to get a better value...

tag1 = sys.argv[1]
print tag1

# read in the data and thin it as necessary
genDict = generator.doEverything()
load_ratio = 1.0
TMSR_REQ = (1.245 * .5)
T10_REQ = (1.245)
T30_REQ = (1.245 + .5 * 1.237)

file_means = csv.reader(open("../load_means.csv", "rU"))
file_means.next()
avgLoads = numpy.array([ float(line[1]) * 1e-3 for line in file_means])
print "Num Loads:\t %d" % len(avgLoads )

#solve a small nominal problem to get the reserve requirements
genDict, load_ratio = generator.smallTestCase( genDict, filt_percent = .1 )
#out = buildNom.buildSolveNom(genDict, 0, 0, 0, avgLoads)
#print "Reserve Reqs:\t", summarize.calcResReqs(out[0], genDict)

TMSR_REQ = 0.4405
T10_REQ = 0.881
T30_REQ = 1.2225

TMSR_REQ = T10_REQ = T30_REQ = 0.


#########
# Set up the output files
#########
file_out_costs = csv.writer(open(tag1 + "_costs.csv", "w"))
file_out_costs.writerow(["Date", "Fixed", "PredVariable", "RealVariable"])

file_out = csv.writer(open(tag1 + "_sched.csv", "w"))
file_out.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )

file_out2 = csv.writer(open(tag1 + "_realized.csv", "w"))
file_out2.writerow([ "Date", "Type"] + ["H" + str(ix + 1) for ix in range(24) ] )


print "\n Num Generators:\t", len(genDict)
for iType in generator.FUEL_TYPES:
    print iType, len( filter( lambda g: g.fuel_type == iType, genDict.values() ) )

print "\n Capacity Requirements:"
print TMSR_REQ, T10_REQ, T30_REQ

print "\n Load Ratio: \t", load_ratio

#load the realized residuals
file_resids = csv.reader(open("../forestResid.csv", "rU") )
file_resids.next() #burn the headers
resids = [line[2:] for line in file_resids]  #burn two for the date and the row name
resids = numpy.array(resids, dtype=float) * 1e-3 * load_ratio
print "Mean Resids:\t", resids.shape
print numpy.mean(resids, axis=0)

#load prediction data by day
file_preds = csv.reader(open("../forestFit.csv", "rU"))
file_preds.next() #burn a header

file_loads = csv.reader(open("../load_validate.csv", "rU") ) 
file_loads.next() # burn one for the header

#model for the second stage stuff
UC2obj =  buildNom.__buildNomNoLoad(genDict, TMSR_REQ, T10_REQ, T30_REQ, False, False)

#the affine cut generator
eps = .1
delta = .2
m = grb.Model("Affine")
affCG = cg.SparseAffineCutGen(resids, eps, 5, m, .5 * delta, .5 * delta)

#build one model for the affine stuff
affModel = buildAff.__buildAffNoLoad(affCG, genDict, TMSR_REQ, T10_REQ, T30_REQ, False, False)

# for each day in the validation set,
for ix, (line_load, line_fit) in enumerate(zip(file_loads, file_preds)):
    if ix > 0:
        sys.exit()

    predLoads = [float(l) * 1e-3  * load_ratio for l in line_fit[2:26] ]

    #Update the affine model and solve
    onVals, startVals, fixedCostAff, totCostAff, prodByHrAff, varCostsAff = \
            buildAff.resolve(affModel, predLoads, genDict)
    summarize.writeHourlySchedCap(file_out, line_load[1], onVals, genDict)

    #extract the period 1 solutions, and solve stage 2
    loads = [float(l) * 1e-3 * load_ratio for l in line_load[2:26]]

    onValsCheck, startValsCheck, fixedCostCheck, totCostRealAff, prodByHrReal, varCostsReal = \
            buildNom.updateSolveSecondStage( UC2obj, loads, genDict, onVals, startVals)

    # dump some statistics for the run
    file_out_costs.writerow([ line_load[1], fixedCostAff, varCostsAff, varCostsReal] )
    summarize.writeHourlyGens(file_out2, line_load[1], prodByHrReal)    
