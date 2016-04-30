------------------------------------------------------------------------------
--
--   Lua - Script for bio-/geomechanics
--
--   Author: Arne Naegel
--           (based on solid_mechanics app by Raphael Prohl)
--
------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("generic.lua")
ug_load_script("cryer.lua")
ug_load_script("mandel.lua")

local kperm = 1e-3 -- m/s
local poro = 0.2
local nu = 0.25

local EYoung = 2.0 * 1e+4 -- kPa
local Kmedium = EYoung/3.0*(1.0-2.0*nu)
local Kfluid = 2.2 * 1e+6 -- kPa

print ("Kmedium = "..Kmedium)
print ("Kfluid  = "..Kfluid)

local problem = cryer2d

problem:init(kperm, poro, nu, 1.0/Kmedium, 1.0/Kfluid, 0.0)



local startTime  = util.GetParamNumber("-start", 0.0, "end time") 
local endTime    = util.GetParamNumber("-end", 1e+5, "end time") 
local dt       = util.GetParamNumber("-dt", 1e-4, "time step size")
local dtMin    = util.GetParamNumber("-dtmin", dt*1e-4, "minimal admissible time step size")
local dtMax    = util.GetParamNumber("-dtmax", dt*1e+3, "minimal admissible time step size")
local dtRed    = util.GetParamNumber("-dtred", 0.5, "time step size reduction factor on divergence")




  local charTime = problem:get_char_time()
  print("charTime="..charTime)
  
  startTime = 0.0
  endTime   = 1.0*charTime
  
  dt  = 1e-4*charTime
  dtMin = 0.1*dt
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
local cpu = problem.cpu or 1

-- ORDER OF ANSATZ-FUNCTIONS 
local porder = problem.porder or 1
local uorder = problem.uorder or 2


InitUG(dim, AlgebraType("CPU", cpu)); -- m=4 fpr block

-- REFINEMENT
-- choose number of pre-Refinements (before sending grid onto different processes)      
local numPreRefs = util.GetParamNumber("-numPreRefs", 0)
-- choose number of total Refinements (incl. pre-Refinements)
local numRefs = util.GetParamNumber("-numRefs", 2) --4



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

local dom = Domain()

-- load grid into domain
LoadDomain(dom, gridName)
print("Loaded domain from " .. gridName)

dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, mandatorySubsets)

--local refiner =  GlobalDomainRefiner(dom)
--refiner:refine();

function createSphereProjector(dom)
     -- args: dom, center, radius, eps
     ProjectVerticesToSphere(dom, {0.0, 0.0}, 0.25, 0.05)  --
     ProjectVerticesToSphere(dom, {0.0, 0.0}, 0.5, 0.05)  --
     return SphericalFalloffProjector(dom, {0.0, 0.0}, 0.5, 0.05) -- For: grid, pos, center, inner_radius, outer_radius --
end

projector = createSphereProjector(dom)

local refProjector = DomainRefinementProjectionHandler(dom)
--refProjector:set_callback("INNER", projector)
SaveDomain(dom, "Sphere.ugx");

SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "CryerSpheres.ugx", 1)

--quit();


-----------------------------------------------------------------
--  Approximation Space
-----------------------------------------------------------------

print("Create ApproximationSpace")
local approxSpace = ApproximationSpace(dom) 
approxSpace:add_fct("p", "Lagrange", porder) 
approxSpace:add_fct("ux", "Lagrange", uorder)          
approxSpace:add_fct("uy", "Lagrange", uorder)   
if (dim==3) then approxSpace:add_fct("uz", "Lagrange", uorder) end


approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)              
print("end approx_init")


local domainDisc = DomainDiscretization(approxSpace)

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Problem Setup
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------





-- Biot
--local rho = 1.0
--local alpha = 1.0
--local F=1.0
--local M = 1.0;



-----------------------------------------------------------------
--  Boundary Conditions & Rhs
-----------------------------------------------------------------

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Setup FE Linear Element Discretization
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

problem:add_elem_discs(domainDisc)
problem:add_boundary_conditions(domainDisc)

-----------------------------------------------------------------
-- Domain discretization
-----------------------------------------------------------------

--[[
myParam = TwoMaterialProblem
myDirichletBnd = ConfinedCompressionBnd

flowEqDiscs = {}
displacementEqDiscs = {}


for i=1,#myParam do
	flowEqDiscs[i], displacementEqDiscs[i] = CreateElemDiscs(myParam[i])
	domainDisc:add(flowEqDiscs[i])
	domainDisc:add(displacementEqDiscs[i])
end
domainDisc:add(myDirichletBnd)

--]]
print("FE discretization-setup. done.")	

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

----------------------------------------
-- create algebraic Preconditioner
----------------------------------------
local jac = Jacobi()
jac:set_damp(0.6)
local gs = GaussSeidel()
local sgs = SymmetricGaussSeidel()
local bgs = BackwardGaussSeidel() 
local ilu = ILU()
--ilu:set_beta(-0.5);
local ilut = ILUT()
--vanka = ElementGaussSeidel(0.9, "vertex")
local vanka = ElementGaussSeidel(1.0, "element")

-------------------------
-- create GMG
-------------------------

local	dbgWriter = GridFunctionDebugWriter(approxSpace)
	
	-- Base Solver
local	baseConvCheck = ConvCheck()
	baseConvCheck:set_maximum_steps(5000)
	baseConvCheck:set_minimum_defect(1e-10)
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
	gmg:set_base_level(0)
	gmg:set_base_solver(superLU)
	gmg:set_smoother(vanka) --(jac)
	gmg:set_cycle_type(1) -- 1:V, 2:W
	gmg:set_num_presmooth(1)
	gmg:set_num_postsmooth(1)
	--gmg:set_debug(dbgWriter)


--------------------------------
-- create and choose a Solver
--------------------------------

local convCheck = ConvCheck()
convCheck:set_maximum_steps(500)
convCheck:set_reduction(1e-10) 
convCheck:set_minimum_defect(1e-12)

local iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilu)
iluSolver:set_convergence_check(convCheck)

local jacSolver = LinearSolver()
jacSolver:set_preconditioner(jac)
jacSolver:set_convergence_check(convCheck)

local gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(convCheck)

local cgSolver = CG()
cgSolver:set_preconditioner(gmg) --(jac)
cgSolver:set_convergence_check(convCheck)

local bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(gmg) --(gmg)
bicgstabSolver:set_convergence_check(convCheck)

local luSolver = SuperLU()

-- choose a solver

local lsolver = luSolver
--solver = jacSolver
--solver = iluSolver
--solver = bicgstabSolver 
--lsolver = gmgSolver
--solver = cgSolver


local u = GridFunction(approxSpace)

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

--newtonCheck = StandardConvCheck(100, 1e-10, 5e-6)
--newtonCheck:set_component_check("ux", 1e-10, 5e-6)
--newtonCheck:set_component_check("uy", 1e-10, 5e-6)
--newtonCheck:set_component_check("p", 1e-10, 5e-6)

local newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(lsolver)
newtonSolver:set_convergence_check(newtonCheck)
newtonSolver:set_line_search(lineSearch)
--newtonSolver:set_debug(dbgWriter)

local nlsolver = newtonSolver

print("Interpolation start values")
problem:interpolate_start_values(u, startTime)


function myPostprocessingCallback(u, step, time)
  problem:post_processing(u, step, time)
end

--step = myPostprocessingCallback
step = vtk



print ("Integrating from 0.0 to "..endTime)
util.SolveLinearTimeProblem(u, domainDisc, lsolver, step, "PoroElasticityTransient",
							   "ImplEuler", 1, startTime, endTime, dt, dtmin, dtred);
							   
--util.SolveNonlinearTimeProblem(u, domainDisc, nlsolver, vtk, "PoroElasticityTransient",
--						   "ImplEuler", 1, startTime, endTime, dt, dtMin, dtRed); 
						   
				   
--util.SolveNonlinearProblemAdaptiveTimestep(u, domainDisc, nlsolver, vtk, "PoroElasticityAdaptive",
--	startTime, endTime, dt, dtMin, dtMax, dtRed, timeTol)
						    
--util.SolveNonlinearProblemAdaptiveLimex(u, domainDisc, nlsolver, vtk, "PoroElasticityAdaptiveLimex",
--	startTime, endTime, dt, dtMin, dtMax, dtRed, timeTol)						   
						   
end
exit();
