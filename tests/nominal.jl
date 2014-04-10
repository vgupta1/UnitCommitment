###
# Unit tests for Nominal solver functions
###
using Base.Test
include("../solvers/solver.jl")

function testFS()
	m = Model()
	ons, starts, stops = addFirstStage(m, HRS=4)
	@addConstraint(m, ons[2] == 1)
	@setObjective(m, Min, sum{ons[i], i = 1:4})
	solve(m)

	true_ons = true_starts = [0 1 0 0]
	true_stops = [0 0 1 0]
	for ix = 1:4
		@test_approx_eq getValue(ons[ix]) true_ons[ix]
		@test_approx_eq getValue(starts[ix]) true_starts[ix]
		@test_approx_eq getValue(stops[ix]) true_stops[ix]
	end
end
testFS()

function testUp()
	m = Model()
	ons, starts, stops = addFirstStage(m, HRS=4)
	@addConstraint(m, ons[2] == 1)
	@setObjective(m, Min, sum{ons[i], i = 1:4} + ons[1])
	addMinUp(m, ons, starts, 2) 
	solve(m)
	true_ons = [0 1 1 0]
	true_starts = [0 1 0 0]
	true_stops = [0 0 0 1]
	for ix = 1:4
		@test_approx_eq getValue(ons[ix]) true_ons[ix]
		@test_approx_eq getValue(starts[ix]) true_starts[ix]
		@test_approx_eq getValue(stops[ix]) true_stops[ix]
	end
end
testUp()

function testDown1()
	m = Model()
	ons, starts, stops = addFirstStage(m, HRS=4)
	@addConstraint(m, ons[2] == 0)
	@setObjective(m, Max, ons[1] + ons[4])
	addMinDown(m, ons, stops, 2)
	solve(m)
	true_ons = true_starts = [1 0 0 1]
	true_stops = [0 1 0 0]
	for ix = 1:4
		@test_approx_eq getValue(ons[ix]) true_ons[ix]
		@test_approx_eq getValue(starts[ix]) true_starts[ix]
		@test_approx_eq getValue(stops[ix]) true_stops[ix]
	end
end
testDown1()

function testDown2()
	m = Model()
	ons, starts, stops = addFirstStage(m, HRS=4)
	@addConstraint(m, ons[2] == 0)
	@setObjective(m, Max, ons[1] + ons[4])
	addMinDown(m, ons, stops, 3)
	solve(m)
	true_ons = [0 0 1 1]
	true_starts = [0 0 1 0]
	true_stops = [0 0 0 0]
	for ix = 1:4
		@test_approx_eq getValue(ons[ix]) true_ons[ix]
		@test_approx_eq getValue(starts[ix]) true_starts[ix]
		@test_approx_eq getValue(stops[ix]) true_stops[ix]
	end
end
testDown2()




