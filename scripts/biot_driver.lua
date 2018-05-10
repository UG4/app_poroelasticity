------------------------------------------------------------------------------
--
--   Lua - Script for poroelacticity
--
--   Author: Arne Naegel
--          (derived from on solid_mechanics app by Raphael Prohl)
--
------------------------------------------------------------------------------

ug_load_script("ug_util.lua")
ug_load_script("plugins/Limex/limex_util.lua")
ug_load_script("generic.lua")
ug_load_script("cryer.lua")
-- ug_load_script("mandel.lua")


-- TIMES AND TIME-STEPPING
local startTime  = util.GetParamNumber("--start", 0.0, "end time") 
local endTime    = util.GetParamNumber("--end", 1e+5, "end time") 
local dtFrac     = util.GetParamNumber("--dtFrac", 1e-5, "time step size")
local dtMinFrac  = util.GetParamNumber("--dtminFrac", 1e-2, "minimal admissible time step size")
local dtMaxFrac  = util.GetParamNumber("--dtmaxFrac", 0.1, "minimal admissible time step size (as fraction of tend)")
local dtRed      = util.GetParamNumber("--dtred", 0.5, "time step size reduction factor on divergence")


-- REFINEMENT
local numPreRefs   = util.GetParamNumber("--numPreRefs", 0, "number of pre-Refinements (before distributing grid)")
local numRefs      = util.GetParamNumber("--num-refs", 3, "total number of refinements (incl. pre-Refinements)") --4 -- 



local params = {
  problemID = util.GetParam("--problem-id", "deleeuw2d"), -- cryer3d‚
  solverID =  util.GetParam("--solver-id", "UzawaMGKrylov"),  --  "FixedStressEX", "UzawaMG", "UzawaSmoother","UzawaMGKrylov"

 
  MGCycleType = util.GetParam("--mg-cycle-type", "W", "V,F,W"),
  MGBaseLevel = util.GetParamNumber("--mg-base-level", 0, "some non-negative integer"),  
  MGNumSmooth = util.GetParamNumber("--mg-num-smooth", 3, "some positive integer"), 
  MGSmootherType =  util.GetParam("--mg-smoother-type", "uzawa", "uzawa,cgs"),
  
  -- LIMEX
  LimexTOL     = util.GetParamNumber("--limex-tol", 1e-3, "TOL"),
  LimexNStages = util.GetParamNumber("--limex-num-stages", 4, "number of LIMEX stages q"),
}



-- Set parameters
local kperm   = 1e-0 -- m/s 1e-3
local poro    = 0.2
local nu      = 0.25

local EYoung  = 2.0 * 1e+2                  -- kPa 2.0 * 1e+4   
local Kmedium = EYoung/(3.0*(1.0-2.0*nu))   
local Kfluid  = 2.2 * 1e+6                  -- kPa -- 2.2 * 1e+6 --

print ("Kmedium = "..Kmedium)
print ("Kfluid  = "..Kfluid)

--deleeuw2d--deleeuw2d -- cryer3d --cryer2d -- mandel3d --, mandel--, cryer3d
local problemList = {
  ["deleeuw2d"] = deleeuw2d,
  ["deleeuw3d"] = deleeuw3d,
  ["deleeuw3dTet"] = deleeuw3dTet,
  ["cryer3d"] = cryer3d,
  ["cryer3dTet"] = cryer3dTet,
}

local problem = problemList[params.problemID]
problem:init(kperm, poro, nu, 1.0/Kmedium, 1.0/Kfluid, 0.0)

local charTime = problem:get_char_time()
print("charTime="..charTime)
  
startTime = 0.0
endTime   = 2.0*charTime
  
local dt  = dtFrac*charTime
local dtMin = dtMinFrac
local dtMax = endTime
  
  
 if (problem == mandel) then 
  local time = 1e-4
  while (time<=10.0) do
    problem:create_test_data(time*charTime)
    time = time *10.0;
  end
 end 


local doSteadyState = true
local doTransient = true

----------------------------------
----------------------------------
--  Settings
----------------------------------
----------------------------------

local dim = problem.dim
local cpu = problem.cpu or dim+1    -- default: block

-- Order for Ansatz functions.
local porder = problem.porder or 1
local uorder = problem.uorder or (porder+1)

InitUG(dim, AlgebraType("CPU", cpu));



-- OUTPUT-ASSISTANT FOR SEVERAL PROCESSES
GetLogAssistant():enable_file_output(true, "output_p_"..ProcRank()..
								"_Lev"..numRefs..".txt")
GetLogAssistant():set_debug_level("SchurDebug", 7);					
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Domain / ApproximationSpace setup
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

-- Create, Load, Refine and Distribute Domain

local gridName = problem.gridName
local mandatorySubsets = problem.mandatorySubsets

local dom = problem:create_domain(numRefs, numPreRefs)


--local refiner =  GlobalDomainRefiner(dom)
--refiner:refine();
--refiner:refine();
-----------------------------------------------------------------
--  Approximation Space
-----------------------------------------------------------------

print("Create ApproximationSpace... ")
local approxSpace = ApproximationSpace(dom) 
approxSpace:add_fct("p", "Lagrange", porder) 

if false then 
  -- Does not work due to registration issues in SmallStrain mechanics.
  uorder=1
  approxSpace:add_fct("ux", "mini", 1)          
  approxSpace:add_fct("uy", "mini", 1)  
else
  local utype = "Lagrange" --"Lagrange"  --"mini" -- "Lagrange"
  approxSpace:add_fct("ux", utype, uorder)          
  approxSpace:add_fct("uy", utype, uorder)   
  if (dim==3) then approxSpace:add_fct("uz", utype, uorder) end
end

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

local uzawaSchurUpdateDisc = DomainDiscretization(approxSpace)
problem:add_uzawa_discs(uzawaSchurUpdateDisc)
print("done!")


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
local u = GridFunction(approxSpace)
local dbgVector=u:clone()
--------------------------------------------------------------------------------

----------------------------------------
-- create algebraic Preconditioner
----------------------------------------
local jac = Jacobi()
jac:set_damp(0.66)
local gs = GaussSeidel()
local sgs = SymmetricGaussSeidel()
local bgs = BackwardGaussSeidel() 
bgs:enable_overlap(true)
gs:enable_overlap(true)
sgs:enable_overlap(true)
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
cgs:set_alpha(1.0)
cgs:set_beta(1.0) --- 0 > 0.25  (beta=0.0: no pressure change) -- 1.0: works
cgs:set_weights(true)


local dbgWriter = GridFunctionDebugWriter(approxSpace)
local uzawaSchurUpdateOp = AssembledLinearOperator()
uzawaSchurUpdateOp:set_discretization(uzawaSchurUpdateDisc)



--- Factory for Uzawa iteration.
-- @function createUzawaIteration
-- @param #string sSchurCmp  Schur complement will be built for this unknown.
-- @param aiForward Approximate Inverse (forward problem)
-- @param aiSchur Approximate Inverse (Schur complement)
-- @param aiBackward Approximate Inverse (backward problem)
function createUzawaIteration(sSchurCmp, aiForward, aiSchur, aiBackward, uzawaSchurUpdateOp, uzawaSchurWeight)

  local uzawa = UzawaBase(sSchurCmp)              
  local weight = uzawaSchurWeight or 1.0
  if (aiForward) then uzawa:set_forward_iter(aiForward)  end
  if (aiSchur) then uzawa:set_schur_iter(aiSchur) end
  if (aiBackward) then uzawa:set_backward_iter(aiBackward)  end

  uzawa:set_schur_operator_update(uzawaSchurUpdateOp, weight)
  -- uzawa:set_debug(dbgWriter)
  
  return uzawa
end

local uzawaWeight = 1.0
local uzawaForward = createUzawaIteration("p", gs, Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
local uzawaBackward = createUzawaIteration("p", nil, Jacobi(0.66), bgs, uzawaSchurUpdateOp, uzawaWeight)

--local uzawaForward = createUzawaIteration("p", Jacobi(0.66), Jacobi(0.66), nil, uzawaSchurUpdateOp, uzawaWeight)
--local uzawaBackward = createUzawaIteration("p", nil, Jacobi(0.66), Jacobi(0.66), uzawaSchurUpdateOp, uzawaWeight)
local uzawaSym = createUzawaIteration("p", gs, sgs, bgs, uzawaSchurUpdateOp, uzawaWeight)
--local uzawaBackward = createUzawaIteration("p", nil, Jacobi(0.5), Jacobi(0.66), uzawaSchurOp, uzawaWeight)
local uzawa = uzawaForward



local preSmoother
local postSmoother

if (params.MGSmootherType == "uzawa") then
  preSmoother = uzawaForward
  postSmoother = uzawaBackward
elseif (params.MGSmootherType == "cgs") then
  preSmoother = cgs
  postSmoother = cgs
else
  quit()

end

-------------------------
-- create GMG
-------------------------

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
gmg:set_base_level(params.MGBaseLevel)  -- was 1 in Cincy
gmg:set_base_solver(superLU)  -- was baseLU in Cincy
gmg:set_presmoother(preSmoother) --(jac)
gmg:set_postsmoother(postSmoother) 
gmg:set_cycle_type(params.MGCycleType) -- 1:V, 2:W -- "F"
gmg:set_num_presmooth(params.MGNumSmooth)
gmg:set_num_postsmooth(params.MGNumSmooth)
gmg:set_rap(true)  -- mandatory, if set_stationary
--gmg:set_debug(dbgWriter)

local gmgP = GeometricMultiGrid(approxSpace)
gmgP:set_discretization(domainDisc)
gmgP:set_base_level(numPreRefs)  -- was 1 in Cincyj
gmgP:set_base_solver(baseLU)  -- was baseLU in Cincy
gmgP:set_presmoother(sgs) 
gmgP:set_postsmoother(sgs) 
gmgP:set_cycle_type("V") -- 1:V, 2:W -- "F"
gmgP:set_num_presmooth(3)
gmgP:set_num_postsmooth(3)
gmgP:set_rap(true)  -- mandatory, if set_stationary
--gmg:set_debug(dbgWriter)

local uzawaTotal       = createUzawaIteration("p", ILUT(1e-8), ILUT(1e-8), nil, uzawaSchurUpdateOp, 1.0)      -- ???
 
local fixedStressLU = createUzawaIteration("p", nil, ILUT(1e-12), ILUT(1e-12), uzawaSchurUpdateOp, 1.0)
local fixedStressMG    = createUzawaIteration("p", nil, ILUT(1e-12), ILUT(1e-12), uzawaSchurUpdateOp, 1.0)


--local transfer = StdTransfer()
--transfer:enable_p1_lagrange_optimization(false)
--gmg:set_transfer(transfer)

--------------------------------
-- debug solver /iter
--------------------------------u
local p0 = problem.modelParameter.p0 or 1.0

local cmpConvCheck = CompositeConvCheck(approxSpace)
cmpConvCheck:set_component_check("ux", p0*1e-14, 1e-6)
cmpConvCheck:set_component_check("uy", p0*1e-14, 1e-6)
if (dim==3) then
cmpConvCheck:set_component_check("uz", p0*1e-14, 1e-6)
end
cmpConvCheck:set_component_check("p", p0*1e-11, 1e-6)
cmpConvCheck:set_maximum_steps(60)


local cmpConvCheck2 = CompositeConvCheck(approxSpace)
  cmpConvCheck2:set_component_check("ux", p0*1e-12, 1e-6)
  cmpConvCheck2:set_component_check("uy", p0*1e-12, 1e-6)
if (dim==3) then
  cmpConvCheck2:set_component_check("uz", p0*1e-12, 1e-6)
end
cmpConvCheck2:set_component_check("p", p0*1e-12, 1e-6)
cmpConvCheck2:set_maximum_steps(50)

cmpConvCheck2 = ConvCheck(2, 1e-10, 1e-8)

local dbgSolver = LinearSolver()
dbgSolver:set_preconditioner(uzawa) -- cgs, gmg, uzawa
dbgSolver:set_convergence_check(cmpConvCheck2)
--dbgSolver:set_debug(dbgWriter)
--dbgSolver:set_convergence_check(cmpConvCheck)

local dbgIter= DebugIterator()
dbgIter:set_preconditioner(gmg)  -- gmg is the 'real' preconditioner
dbgIter:set_solver(dbgSolver)
dbgIter:set_solution(dbgVector)
dbgIter:set_random_bounds(-5e-6, 5e-6)
dbgIter:set_debug(dbgWriter)  -- print t_0 anf t_N

--------------------------------
-- create and choose a Solver
--------------------------------

local solver = {}

local convCheck = ConvCheck()
convCheck:set_maximum_steps(500)
convCheck:set_reduction(1e-8) 
convCheck:set_minimum_defect(1e-16)
--convCheck = cmpConvCheck

local iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilut)
iluSolver:set_convergence_check(convCheck)

local jacSolver = LinearSolver()
jacSolver:set_preconditioner(jac)
jacSolver:set_convergence_check(convCheck)

solver["UzawaSmoother"] = LinearSolver()
solver["UzawaSmoother"]:set_preconditioner(uzawaForward)
solver["UzawaSmoother"]:set_convergence_check(convCheck)

solver["UzawaMG"] = LinearSolver()
solver["UzawaMG"]:set_preconditioner(gmg) -- gmg, dbgIter
solver["UzawaMG"]:set_convergence_check(convCheck) -- cmpConvCheck

solver["FixedStressEX"] = BiCGStab()
solver["FixedStressEX"]:set_preconditioner(fixedStressLU)
solver["FixedStressEX"]:set_convergence_check(convCheck)

solver["UzawaMGKrylov"] = BiCGStab()
solver["UzawaMGKrylov"]:set_preconditioner(gmg) -- gmg, dbgIter
solver["UzawaMGKrylov"]:set_convergence_check(convCheck) -- cmpConvCheck

local bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(dbgIter) --(gmg)
bicgstabSolver:set_convergence_check(convCheck)

local cgSolver = CG()
cgSolver:set_preconditioner(dbgIter) --(gmg)
cgSolver:set_convergence_check(convCheck)


local gmresSolver = GMRES(3)
gmresSolver:set_preconditioner(gmg) -- gmg, dbgIter
gmresSolver:set_convergence_check(convCheck)

local sluSolver = SuperLU()

local luSolver = LinearSolver()
luSolver:set_preconditioner(LU())
luSolver:set_convergence_check(convCheck)

-- Select solver.
local lsolver = solver[params.solverID]
--solver = jacSolver
--lsolver = iluSolver
--lsolver = gmgSolver
--lsolver = cgSolver
--lsolver = bicgstabSolver 
--lsolver = gmresSolver
--lsolver:set_compute_fresh_defect_when_finished(true)
-- lsolver = sluSolver

local vtk=VTKOutput()
vtk:select_nodal("p", "PNodal")
if (dim == 2) then vtk:select({"ux", "uy"}, "uNodal") end
if (dim == 3) then vtk:select({"ux", "uy", "uz"}, "uNodal") end

--vtk:select_element( displacementEqDisc:displacement(), "DispElem")
--vtk:select_element( displacementEqDisc:divergence(), "DivElem")
--vtk:select_element( flowEqDisc:gradient()dbgSolver, "GradP")
--vtk:select(massLinker, "Mass")



-- Init error estimator.
local biotErrorEst = ScaledGridFunctionEstimator()
biotErrorEst:add(L2ComponentSpace("p", 2))        -- L2 norm for p, 2nd order quadrature
--[[ 
biotErrorEst:add(H1ComponentSpace("ux", 4))  
biotErrorEst:add(H1ComponentSpace("uy", 4)) 
if (dim==3) then biotErrorEst:add(H1ComponentSpacer("uz", 4)) end
--]]


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Solve transient (linear) problem
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------	
if (doTransient) then

local lineSearch = StandardLineSearch();
lineSearch:set_maximum_steps(6)
lineSearch:set_accept_best(true)

local newtonCheck = ConvCheck()
newtonCheck:set_maximum_steps(10)
newtonCheck:set_minimum_defect(1e-14)
newtonCheck:set_reduction(5e-6)
newtonCheck:set_verbose(true)

local newtonCheck2 = CompositeConvCheck(approxSpace)
newtonCheck2:set_component_check("ux", p0*1e-7, 5e-6)
newtonCheck2:set_component_check("uy", p0*1e-7, 5e-6)
if (dim==3) then
newtonCheck2:set_component_check("uz", p0*1e-7, 5e-6)
end
newtonCheck2:set_component_check("p", p0*1e-9, 5e-6)
newtonCheck2:set_maximum_steps(2)

local newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(lsolver)
newtonSolver:set_convergence_check(newtonCheck)
--newtonSolver:set_line_search(lineSearch)
--newtonSolver:set_debug(dbgWriter)

local nlsolver = newtonSolver

print(lsolver:config_string())

print("Interpolation start values")
problem:interpolate_start_values(u, startTime)


function myStepCallback(u, step, time)
  problem:post_processing(u, step, time)
  vtk:print("Cryer.vtu", u, step, time)
end

print ("Integrating from 0.0 to "..endTime)
--util.SolveLinearTimeProblem(u, domainDisc, lsolver, myStepCallback, "PoroElasticityTransient",
--							   "ImplEuler", 1, startTime, endTime, dt, dtmin, dtred);
							   
--dt =dt*1e-4*problem:get_char_time() -- smaller => more complicated

dt = 1e-4*problem:get_char_time()
dtMin = 1e-2*dt

if (false) then

print("dt="..dt/charTime)							   
util.SolveNonlinearTimeProblem(u, domainDisc, nlsolver, myStepCallback, "PoroElasticityTransient",
						   "ImplEuler", 1, startTime, endTime, dt, dtMin/2, dtRed); 
end
						   
local cAdaptiveStepInfo  ={
      ["TOLERANCE"] = 2e-2, 
      ["REDUCTION"] = 0.5, 
      ["INCREASE"]  = 1.2, 
      ["SAFETY"]    = 0.5,
      ["ESTIMATOR"] = biotErrorEst
}
				   

if (true) then

-- Adjust NEWTON for LIMEX.
newtonCheck:set_maximum_steps(1)
newtonCheck:set_supress_unsuccessful(true) 

-- Create & configure LIMEX
-- LIMEX descriptor.
local limexDesc = {
  nstages = params.LimexNStages,
  steps = {1,2,3,4,5,6,7,8,9,10},
  domainDisc=domainDisc,
 
  nonlinSolver = nlsolver,
  tol = params.LimexTOL,
  dt = dt,
  dtmin = dtMin,
  dtmax = dtMax,
  
   -- gammaDiscOPT= gammaTensorDisc,  -- no gamma for linear problem
}


-- Call factory.
local limex = util.limex.CreateIntegrator(limexDesc)
limex:add_error_estimator(biotErrorEst)
-- limex:set_tolerance(0.001)
-- limex:set_time_step(dt)
-- limex:set_dt_min(dtMin)
-- limex:set_dt_max(dtMax)
limex:set_stepsize_safety_factor(0.8)
limex:set_stepsize_greedy_order_factor(0.0)

limex:disable_matrix_cache()        -- This problem is linear
--limex:set_time_derivative(udot)   -- 

-- Create observers.
local vtkFull = VTKOutput()
local vtkobserver = VTKOutputObserver("PoroElasticityTransient.vtk", vtkFull)

local luaobserver = LuaCallbackObserver()
function myLuaPostProcess(step, time, currdt)
  local usol=luaobserver:get_current_solution()
  problem:post_processing(usol, step, time)
  return 0;
end

luaobserver:set_callback("myLuaPostProcess")


-- Attach observers.
-- if (limex_output==1) then
   limex:attach_observer(luaobserver)
   limex:attach_observer(vtkobserver)
--end


-- Solve problem with LIMEX.
local myclock = CuckooClock()
myclock:tic()
limex:apply(u, endTime, u, startTime)
print("CDELTA="..myclock:toc())
end

end -- doTransient

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Solve linear, steady state problem (easy!)
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
  vtk:print("PoroElasticitySteadyState", u, 0, 0.0)


end  --  doSteadyState

				   

