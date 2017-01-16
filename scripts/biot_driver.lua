------------------------------------------------------------------------------
--
--   Lua - Script for poroelacticity
--
--   Author: Arne Naegel
--          (derived from on solid_mechanics app by Raphael Prohl)
--
------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
--ug_load_script("solver_util.lua")
ug_load_script("generic.lua")
ug_load_script("cryer.lua")
ug_load_script("mandel.lua")

local kperm = 1e-0 -- m/s 1e-3
local poro = 0.2
local nu = 0.25

local EYoung = 2.0 * 1e+2 -- kPa 2.0 * 1e+4 
local Kmedium = EYoung/3.0*(1.0-2.0*nu)
local Kfluid = 2.2 * 1e+6 -- kPa -- 2.2 * 1e+6 --



print ("Kmedium = "..Kmedium)
print ("Kfluid  = "..Kfluid)

local problem = mandel3d --, mandel--, cryer3d

problem:init(kperm, poro, nu, 1.0/Kmedium, 1.0/Kfluid, 0.0)



local startTime  = util.GetParamNumber("-start", 0.0, "end time") 
local endTime    = util.GetParamNumber("-end", 1e+5, "end time") 
local dt       = util.GetParamNumber("-dt", 1e-3, "time step size")
local dtMin    = util.GetParamNumber("-dtmin", dt*1e-2, "minimal admissible time step size")
local dtMax    = util.GetParamNumber("-dtmax", dt*1e+3, "minimal admissible time step size")
local dtRed    = util.GetParamNumber("-dtred", 0.5, "time step size reduction factor on divergence")




local charTime = problem:get_char_time()
print("charTime="..charTime)
  
startTime = 0.0
endTime   = 1.0*charTime
  
dt  = 1e-0*charTime
dtMin = 1e-4*dt
dtMax = 0.1*endTime
  
  
 if (problem == mandel) then 
  local time = 1e-4
  while (time<=10.0) do
    problem:create_test_data(time*charTime)
    time = time *10.0;
  end
 end 


local doSteadyState = false
local doTransient = true

----------------------------------
----------------------------------
--  Settings
----------------------------------
----------------------------------


local dim = problem.dim
local cpu = problem.cpu or 3

-- ORDER OF ANSATZ-FUNCTIONS 
local porder = problem.porder or 1
local uorder = problem.uorder or (porder+1)


InitUG(dim, AlgebraType("CPU", cpu)); -- m=4 fpr block

-- REFINEMENT
-- choose number of pre-Refinements (before sending grid onto different processes)      
local numPreRefs = util.GetParamNumber("-numPreRefs", 1)
-- choose number of total Refinements (incl. pre-Refinements)
local numRefs = util.GetParamNumber("-numRefs", 1) --4



local epsilon    = util.GetParamNumber("-epsilon", 1e-4, "epsilon")
local timeTol    = util.GetParamNumber("-tol", 1e-2, "TOL")




-- OUTPUT-ASSISTANT FOR SEVERAL PROCESSES
GetLogAssistant():enable_file_output(true, "output_p_"..ProcRank()..
								"_Lev"..numRefs..".txt")
								
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Domain / ApproximationSpace setup
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

-- Create, Load, Refine and Distribute Domain

local gridName = problem.gridName
local mandatorySubsets = problem.mandatorySubsets



-- load grid into domain
--LoadDomain(dom, gridName)
--print("Loaded domain from " .. gridName)

local dom = problem:create_domain(numRefs, numPreRefs)


--quit();

--local refiner =  GlobalDomainRefiner(dom)
--refiner:refine();
--refiner:refine();
-----------------------------------------------------------------
--  Approximation Space
-----------------------------------------------------------------

print("Create ApproximationSpace... ")
local approxSpace = ApproximationSpace(dom) 
approxSpace:add_fct("p", "Lagrange", porder) 
approxSpace:add_fct("ux", "Lagrange", uorder)          
approxSpace:add_fct("uy", "Lagrange", uorder)   
if (dim==3) then approxSpace:add_fct("uz", "Lagrange", uorder) end

approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)              
print("... done!")


--------------------------------------------------------------------------------
-- Problem Setup
--------------------------------------------------------------------------------
print("FE discretization...") 
local bSteadyStateMechanics = true -- true
local domainDisc = DomainDiscretization(approxSpace)
problem:add_elem_discs(domainDisc, bSteadyStateMechanics)
problem:add_boundary_conditions(domainDisc, bSteadyStateMechanics)
print("done!")


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
local u = GridFunction(approxSpace)
local dbgSmooth=u:clone()
--------------------------------------------------------------------------------

----------------------------------------
-- create algebraic Preconditioner
----------------------------------------
local jac = Jacobi()
jac:set_damp(0.66)
local gs = GaussSeidel()
local sgs = SymmetricGaussSeidel()
local bgs = BackwardGaussSeidel() 
local ilu = ILU()
--ilu:set_beta(-0.5);
local ilut = ILUT()
ilut:set_threshold(1e-3)

--local egs_weights = u:clone();
--egs_weights:set(1.0);
--Interpolate(0.1, egs_weights, "p")

local egs = ElementGaussSeidel() -- patches per node
egs:select_schur_cmp({"p"}, 4.0)
egs:set_relax(0.125)

local cgs = ComponentGaussSeidel(1.0, {"p"}) -- patches per node
cgs:set_alpha(2.0)
cgs:set_beta(0.5) --- 0 > 0.25  (beta=0.0: no pressure change) -- 1.0: works
cgs:set_weights(true)
--cgs:set_damp(0.25)

gs = egs
-------------------------
-- create GMG
-------------------------

local	dbgWriter = GridFunctionDebugWriter(approxSpace)
	
	-- Base Solver
local	baseConvCheck = ConvCheck()
	baseConvCheck:set_maximum_steps(5000)
	baseConvCheck:set_reduction(1e-12)
	baseConvCheck:set_verbose(false)

local	base = BiCGStab()
	base:set_preconditioner(jac)
	base:set_convergence_check(baseConvCheck)
	
local	baseCG = CG()
	baseCG:set_preconditioner(jac)
	baseCG:set_convergence_check(baseConvCheck)
	
	-- exact base solver
local	baseLU = LU()
local	superLU = SuperLU()



	-- Geometric Multi Grid
local	gmg = GeometricMultiGrid(approxSpace)
gmg:set_discretization(domainDisc)
gmg:set_base_level(0)  -- was 1 in Cincyj
gmg:set_base_solver(baseLU)  -- was baseLU in Cincy
gmg:set_smoother(jac) --(jac)
gmg:set_cycle_type("F") -- 1:V, 2:W -- "F"
gmg:set_num_presmooth(2)
gmg:set_num_postsmooth(2)
gmg:set_rap(true)  -- mandatory, if set_stationary
--gmg:set_debug(dbgWriter)

--local transfer = StdTransfer()
--transfer:enable_p1_lagrange_optimization(false)
--gmg:set_transfer(transfer)

--------------------------------
-- debug solver /iter
--------------------------------
local p0 = problem.modelParameter.p0 or 1.0

local cmpConvCheck = CompositeConvCheck(approxSpace)
cmpConvCheck:set_component_check("ux", p0*1e-7, 1e-6)
cmpConvCheck:set_component_check("uy", p0*1e-7, 1e-6)
if (dim==3) then
cmpConvCheck:set_component_check("uz", p0*1e-7, 1e-6)
end
cmpConvCheck:set_component_check("p", p0*1e-9, 1e-6)
cmpConvCheck:set_maximum_steps(60)


local cmpConvCheck2 = CompositeConvCheck(approxSpace)
cmpConvCheck2:set_component_check("ux", p0*1e-7, 1e-6)
cmpConvCheck2:set_component_check("uy", p0*1e-7, 1e-6)
if (dim==3) then
cmpConvCheck2:set_component_check("uz", p0*1e-7, 1e-6)
end
cmpConvCheck2:set_component_check("p", p0*1e-9, 1e-6)
cmpConvCheck2:set_maximum_steps(20)

cmpConvCheck2 = ConvCheck(20, 1e-10, 1e-8)

local dbgSolver = LinearSolver()
dbgSolver:set_preconditioner(gs) -- cgs, gmg
dbgSolver:set_convergence_check(cmpConvCheck2)
--dbgSolver:set_debug(dbgWriter)
--dbgSolver:set_convergence_check(cmpConvCheck)

local dbgIter= DebugIterator()
dbgIter:set_preconditioner(gmg)  -- gmg is the 'real' preconditioner
dbgIter:set_solver(dbgSolver)
dbgIter:set_solution(dbgSmooth)
dbgIter:set_random_bounds(-5e-6, 5e-6)
dbgIter:set_debug(dbgWriter)  -- print t_0 anf t_N

--------------------------------
-- create and choose a Solver
--------------------------------

local convCheck = ConvCheck()
convCheck:set_maximum_steps(500)
convCheck:set_reduction(1e-8) 
convCheck:set_minimum_defect(1e-10)
--convCheck = cmpConvCheck

local iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilut)
iluSolver:set_convergence_check(convCheck)

local jacSolver = LinearSolver()
jacSolver:set_preconditioner(jac)
jacSolver:set_convergence_check(convCheck)

local gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(dbgIter) -- gmg, dbgIter
gmgSolver:set_convergence_check(convCheck)

local bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(dbgIter) --(gmg)
bicgstabSolver:set_convergence_check(convCheck)

local gmresSolver = GMRES(3)
gmresSolver:set_preconditioner(gmg) -- gmg, dbgIter
gmresSolver:set_convergence_check(convCheck)


local luSolver = SuperLU()
--local luSolver = LinearSolver()
--luSolver:set_preconditioner(LU())
--luSolver:set_convergence_check(convCheck)

-- choose a solver

local lsolver = luSolver
--solver = jacSolver
--lsolver = iluSolver
--lsolver = gmgSolver
lsolver = bicgstabSolver 

--lsolver = gmresSolver

--lsolver:set_compute_fresh_defect_when_finished(true)

local vtk=VTKOutput()
vtk:select_nodal("p", "PNodal")
if (dim == 2) then vtk:select({"ux", "uy"}, "uNodal") end
if (dim == 3) then vtk:select({"ux", "uy", "uz"}, "uNodal") end

--vtk:select_element( displacementEqDisc:displacement(), "DispElem")
--vtk:select_element( displacementEqDisc:divergence(), "DivElem")
--vtk:select_element( flowEqDisc:gradient(), "GradP")
--vtk:select(massLinker, "Mass")

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Solving Procedure (linear, steady state)
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------	


if (doSteadyState) then
local A = MatrixOperator()
u = GridFunction(approxSpace)
local b = GridFunction(approxSpace)
u:set(0.0)
b:set(0.0)
-- 1. assemble matrix and rhs
domainDisc:assemble_linear(A, b)

-- 2. set dirichlet values in start iterate
u:set(0.0)
domainDisc:adjust_solution(u)

-- 3. init solver for linear Operator
lsolver:init(A)

SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")

-- 4. apply solver
u:set_random(0.0, 1.0)
lsolver:apply_return_defect(u,b)
vtk:print("PoroElasticityStationary", u, 0, 0.0)


end



--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Solving Procedure (linear, transient)
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------	
if (doTransient) then

local lineSearch = StandardLineSearch();
lineSearch:set_maximum_steps(6)
--lineSearch:set_accept_best(true)

local newtonCheck = ConvCheck()
newtonCheck:set_maximum_steps(100)
newtonCheck:set_minimum_defect(1e-7)
newtonCheck:set_reduction(5e-6)
newtonCheck:set_verbose(true)

local newtonCheck2 = CompositeConvCheck(approxSpace)
newtonCheck2:set_component_check("ux",  p0*1e-7, 5e-6)
newtonCheck2:set_component_check("uy",p0*1e-7, 5e-6)
if (dim==3) then
newtonCheck2:set_component_check("uz",p0*1e-7, 5e-6)
end
newtonCheck2:set_component_check("p", p0*1e-9, 5e-6)
newtonCheck2:set_maximum_steps(2)

local newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(lsolver)
newtonSolver:set_convergence_check(newtonCheck)
--newtonSolver:set_line_search(lineSearch)
--newtonSolver:set_debug(dbgWriter)

local nlsolver = newtonSolver

print("Interpolation start values")
problem:interpolate_start_values(u, startTime)


function myStepCallback(u, step, time)
  problem:post_processing(u, step, time)
  vtk:print("Cryer.vtu", u, step, time)
end





print ("Integrating from 0.0 to "..endTime)
--util.SolveLinearTimeProblem(u, domainDisc, lsolver, myStepCallback, "PoroElasticityTransient",
--							   "ImplEuler", 1, startTime, endTime, dt, dtmin, dtred);
							   
dt =dt*1e-4 -- smaller => more complicated
print("dt="..dt/charTime)							   
util.SolveNonlinearTimeProblem(u, domainDisc, nlsolver, myStepCallback, "PoroElasticityTransient",
						   "ImplEuler", 1, startTime, endTime, dt, dtMin/2, dtRed); 
						   
local cAdaptiveStepInfo  ={
      ["TOLERANCE"] = 2e-2, 
      ["REDUCTION"] = 0.5, 
      ["INCREASE"]  = 1.2, 
      ["SAFETY"]    = 0.5,
      ["ESTIMATOR"] = GridFunctionEstimator("p", 2) -- compare p-error (2nd order quadrature)
}
				   
--util.SolveNonlinearProblemAdaptiveTimestep(u, domainDisc, nlsolver, myStepCallback, "PoroElasticityAdaptive", 
--   startTime, endTime, dt, dtMin, dtMax, cAdaptiveStepInfo)
						
						

						    
--util.SolveNonlinearProblemAdaptiveLimex(u, domainDisc, nlsolver, myStepCallback, "PoroElasticityAdaptiveLimex",
--	startTime, endTime, dt, dtMin, dtMax, cAdaptiveStepInfo)						   
						   
end
exit();
