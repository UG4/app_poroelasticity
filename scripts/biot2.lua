------------------------------------------------------------------------------
--
--   Lua - Script for bio-/geomechanics
--
--   Author: Arne Naegel
--           (based on solid_mechanics app by Raphael Prohl)
--
------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = 3

-----------------------------------------------------------------
--  Settings
-----------------------------------------------------------------

-- ORDER OF ANSATZ-FUNCTIONS 
porder = 1
uorder = 1
InitUG(dim, AlgebraType("CPU", 4));

-- REFINEMENT
-- choose number of pre-Refinements (before sending grid onto different processes)      
numPreRefs = util.GetParamNumber("-numPreRefs", 1)
-- choose number of total Refinements (incl. pre-Refinements)
numRefs = util.GetParamNumber("-numRefs", 2) --4

startTime  = util.GetParamNumber("-start", 0.0, "end time") 
endTime    = util.GetParamNumber("-end", 100.0, "end time") 
dt 		   = util.GetParamNumber("-dt", 1.0, "time step size")
dtMin	   = util.GetParamNumber("-dtmin", 1e-4, "minimal admissible time step size")
dtMax	   = util.GetParamNumber("-dtmax", 1e+2, "minimal admissible time step size")
dtRed	   = util.GetParamNumber("-dtred", 0.5, "time step size reduction factor on divergence")
timeTol    = 1e-4



-- OUTPUT-ASSISTANT FOR SEVERAL PROCESSES
GetLogAssistant():enable_file_output(true, "output_p_"..ProcRank()..
								"_Lev"..numRefs..".txt")
								
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Domain / ApproximationSpace setup
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

-- Create, Load, Refine and Distribute Domain

-- grid:
--gridName = util.GetParam("-grid", "../grids/plate_with_hole_100.ugx");
--neededSubsets = {"Inner", "Top", "SymY", "SymX"};

--gridName = util.GetParam("-grid", "../grids/mandel.ugx");
--gridName = util.GetParam("-grid", "../grids/tunnel.ugx");
--gridName = util.GetParam("-grid", "../grids/mandel.ugx");
gridName = util.GetParam("-grid", "../grids/cube3.ugx");
neededSubsets = {"INNER", "TOP", "BOTTOM", "LEFT", "RIGHT", "SRC"};
neededSubsets = {"INNER", "TOP", "BOTTOM", "LEFT", "RIGHT", "SRC"};

dom = Domain()

-- load grid into domain
LoadDomain(dom, gridName)
print("Loaded domain from " .. gridName)

dom = util.CreateAndDistributeDomain(gridName, numRefs, numPreRefs, neededSubsets)

-----------------------------------------------------------------
--  Approximation Space
-----------------------------------------------------------------

print("Create ApproximationSpace")
approxSpace = ApproximationSpace(dom) 
approxSpace:add_fct("p", "Lagrange", porder) 
approxSpace:add_fct("ux", "Lagrange", uorder)          
approxSpace:add_fct("uy", "Lagrange", uorder)   
if (dim==3) then approxSpace:add_fct("uz", "Lagrange", uorder) end


approxSpace:init_levels()
approxSpace:init_top_surface()
approxSpace:print_statistic()
approxSpace:print_local_dof_statistic(2)              
print("end approx_init")


--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
-- Problem Setup
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


lambda = 1.0
mu = 1.0
Lambda = 1.0
epsilon = 1.0
F=1.0
alpha = 1.0


rho = 1.0
M = 1.0;


-----------------------------------------------------------------
--  Boundary Conditions & Rhs
-----------------------------------------------------------------

if (dim==2) then
-- zero dirichlet boundary conditions


	function Dirichlet0(x, y, t) return true, 0.0 end      -- p  (left, right)
	function UxDirichletTop(x, y, t) return true, 0.0 end   -- ux (bottom)
	--function UyDirichletTop(x, y, t) return true, (1.0-math.exp(-t*0.1))*(x*x-25.0)*0.04 end   -- uy (bottom)
	function UyDirichletTop(x, y, t) return true, 1.0 end   -- uy (bottom)
	function UyDirichletBottom(x, y, t) return true, 0.0 end   -- uy (bottom)
	function PointSource(x, y, t) return true, 0.0 end      -- p  (left, right)
	
	function ThreeRegionFlowPerm(x, y, t, si) 
	if math.abs(y)<2.5 then return Lambda*epsilon, 0, 0, Lambda*epsilon
	else return Lambda, 0, 0, Lambda end
	
	
	function ThreeRegionElastLambda(x, y, t, si) 
	if x>2.5 then return 1.0/epsilon
	else return 1.0 end
end
end

end

if  (dim==3) then
	function Dirichlet0(x, y, z, t) return true, 0.0 end  
	function PointSource(x, y, z, t) return true, 1.0 end     
	
	function ThreeRegionFlowPerm(x, y, z, t, si) 
	if math.abs(y)<0.5 then 
	
	return Lambda*epsilon, 0, 0, 
			0, Lambda*epsilon, 0, 
			0,  0, Lambda*epsilon
	else return Lambda, 0, 0,  
				0, Lambda, 0, 
				0, 0, Lambda 
	end
	
	function ThreeRegionElastLambda(x, y, z, t, si) 
	if x>0.5 then return 1.0/epsilon
	else return 1.0 end
end
end
end 





--function ourVolumeForceField2d(x, y, t)
--	return 0.0, 0.0 
--end

--volumeForceField = LuaUserVector("ourVolumeForceField2d")

--dirichletZero = LuaCondUserNumber("DirichletBndZero" .. dim .. "d")
--neumannBndVal = LuaCondUserNumber("NeumannBnd" .. dim .. "d")

mandelDirichletBnd = DirichletBoundary()
mandelDirichletBnd:add(LuaCondUserNumber("Dirichlet0"), "p", "TOP")
--mandelDirichletBnd:add(LuaCondUserNumber("PointSource"), "p", "SRC")
--mandelDirichletBnd:add(LuaCondUserNumber("UxDirichletTop"), "ux", "TOP")
--mandelDirichletBnd:add(LuaCondUserNumber("UyDirichletTop"), "uy", "TOP")
mandelDirichletBnd:add(LuaCondUserNumber("Dirichlet0"), "ux", "BOTTOM")
mandelDirichletBnd:add(LuaCondUserNumber("Dirichlet0"), "uy", "BOTTOM")
if (dim==3) then mandelDirichletBnd:add(LuaCondUserNumber("Dirichlet0"), "uz", "BOTTOM") end 
--mandelDirichletBnd:add(LuaCondUserNumber("Dirichlet0"), "p", "BOTTOM")


dirichletBnd = mandelDirichletBnd
--neumannBnd = NeumannBoundaryFE("ux")
--neumannBnd:add(neumannBndVal, "Top", "Inner")

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Setup FE Linear Element Discretization
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------


-- define eqns for pressure,  displacement
flowEqDisc = ConvectionDiffusion("p", "INNER, INNER2", "fe")
if (dim==2) then displacementEqDisc = SmallStrainMechanics2d("ux,uy", "INNER, INNER2") end
if (dim==3) then displacementEqDisc = SmallStrainMechanics3d("ux,uy,uz", "INNER, INNER2") end

--matLaw:set_hooke_elasticity_tensor_E_nu(E, nu)
matLaw = HookeLaw()
matLaw:set_hooke_elasticity_tensor(lambda, mu)
displacementEqDisc:set_material_law(matLaw)


-- specify eqn (for displacement)
forceLinker = ScaleAddLinkerVector()
forceLinker:add(-alpha, flowEqDisc:gradient())
displacementEqDisc:set_volume_forces(forceLinker)

--scalarLinker = ScaleAddLinkerNumber()
--scalarLinker:add(alpha, flowEqDisc:value())
--displacementEqDisc:set_pressure(scalarLinker)

--mechOut = MechOutputWriter()
--displacementEqDisc:set_output_writer(mechOut)
--displacementEqDisc:displacement()


-- specify flow eqn (for pressure)
compressionLinker = ScaleAddLinkerNumber()
compressionLinker:add(alpha, displacementEqDisc:divergence())
flowEqDisc:set_mass_scale(1.0/M);
flowEqDisc:set_mass(compressionLinker);
flowEqDisc:set_diffusion(0.01);

if (porder==1) then flowEqDisc:set_quad_order(2) end
if (uorder==1) then displacementEqDisc:set_quad_order(2) end
-- print info

--print(flowEqDisc:config_string())


if dim == 3 then
	--if (order == 1) then
	--	displacementEqDisc:set_quad_order(3)  
		--3:#ip`s: 6, 2:#ip`s: 8 
	--	end
	if (uorder == 2) then
		displacementEqDisc:set_quad_order(7) 
		flowEqDisc:set_quad_order(7)
		--#ip`s: 31
		end
	--if (order == 3) then	
	--	elemDisc:set_quad_order(11) 
		--#ip`s: 90
	--	end
	--if (order > 3) then	
	--	elemDisc:set_quad_order(11) 
		--#ip`s: 90
	--	end
end
--print(flowEqDisc:config_string())
print(displacementEqDisc:config_string())

massLinker = ScaleAddLinkerNumber()
massLinker:add(rho/M, flowEqDisc:value())
massLinker:add(rho*alpha, displacementEqDisc:divergence())

-----------------------------------------------------------------
-- Domain discretization
-----------------------------------------------------------------

domainDisc = DomainDiscretization(approxSpace)
domainDisc:add(displacementEqDisc)
domainDisc:add(flowEqDisc)
--domainDisc:add(neumannBnd)
domainDisc:add(dirichletBnd)

print("FE discretization-setup. done.")	

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Algebra
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

----------------------------------------
-- create algebraic Preconditioner
----------------------------------------
jac = Jacobi()
jac:set_damp(0.6)
gs = GaussSeidel()
sgs = SymmetricGaussSeidel()
bgs = BackwardGaussSeidel() 
ilu = ILU()
--ilu:set_beta(-0.5);
ilut = ILUT()

-------------------------
-- create GMG
-------------------------

	dbgWriter = GridFunctionDebugWriter(approxSpace)
	
	-- Base Solver
	baseConvCheck = ConvCheck()
	baseConvCheck:set_maximum_steps(5000)
	baseConvCheck:set_minimum_defect(1e-10)
	baseConvCheck:set_reduction(1e-12)
	baseConvCheck:set_verbose(false)

	base = BiCGStab()
	base:set_preconditioner(jac)
	base:set_convergence_check(baseConvCheck)
	
	baseCG = CG()
	baseCG:set_preconditioner(jac)
	baseCG:set_convergence_check(baseConvCheck)
	
	-- exact base solver
	baseLU = LU()

	-- Geometric Multi Grid
	gmg = GeometricMultiGrid(approxSpace)
	gmg:set_discretization(domainDisc)
	gmg:set_base_level(0)
	gmg:set_base_solver(baseLU)
	gmg:set_smoother(sgs) --(jac)
	gmg:set_cycle_type(1) -- 1:V, 2:W
	gmg:set_num_presmooth(2)
	gmg:set_num_postsmooth(2)
	--gmg:set_debug(dbgWriter)

-- Exact LU solver
luSolver = LU()

--------------------------------
-- create and choose a Solver
--------------------------------

convCheck = ConvCheck()
convCheck:set_maximum_steps(500)
convCheck:set_reduction(1e-10) 
convCheck:set_minimum_defect(1e-12)

iluSolver = LinearSolver()
iluSolver:set_preconditioner(ilu)
iluSolver:set_convergence_check(convCheck)

jacSolver = LinearSolver()
jacSolver:set_preconditioner(jac)
jacSolver:set_convergence_check(convCheck)

gmgSolver = LinearSolver()
gmgSolver:set_preconditioner(gmg)
gmgSolver:set_convergence_check(convCheck)

cgSolver = CG()
cgSolver:set_preconditioner(ilu) --(jac)
cgSolver:set_convergence_check(convCheck)

bicgstabSolver = BiCGStab()
bicgstabSolver:set_preconditioner(ilu) --(gmg)
bicgstabSolver:set_convergence_check(convCheck)

luSolver = LU()

-- choose a solver

--lsolver = luSolver
--solver = jacSolver
--solver = iluSolver
--solver = bicgstabSolver 
lsolver = gmgSolver
--solver = cgSolver


u = GridFunction(approxSpace)

vtk=VTKOutput()
vtk:select_nodal("p", "PNodal")
if (dim == 2) then vtk:select({"ux", "uy"}, "DispNodal") end
if (dim == 3) then vtk:select({"ux", "uy", "uz"}, "DispNodal") end

vtk:select_element( displacementEqDisc:displacement(), "DispElem")
vtk:select_element( displacementEqDisc:divergence(), "DivElem")
vtk:select_element( flowEqDisc:gradient(), "GradP")
vtk:select(massLinker, "Mass")
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Solving Procedure (linear, steady state)
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------	
if (true) then
A = MatrixOperator()
u = GridFunction(approxSpace)
b = GridFunction(approxSpace)
u:set(0.0)
b:set(0.0)
-- 1. assemble matrix and rhs
domainDisc:assemble_linear(A, b)

-- 2. set dirichlet values in start iterate
u:set(0.0)
domainDisc:adjust_solution(u)

-- 3. init solver for linear Operator
lsolver:init(A)

--SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")

-- 4. apply solver
lsolver:apply_return_defect(u,b)
vtk:print("MandelStationary", u, 0, 0.0)
end




lineSearch = StandardLineSearch();
lineSearch:set_maximum_steps(6)
--lineSearch:set_accept_best(true)

newtonCheck = ConvCheck()
newtonCheck:set_maximum_steps(100)
newtonCheck:set_minimum_defect(1e-7)
newtonCheck:set_reduction(5e-6)
newtonCheck:set_verbose(true)

--newtonCheck = StandardConvCheck(100, 1e-10, 5e-6)
--newtonCheck:set_component_check("ux", 1e-10, 5e-6)
--newtonCheck:set_component_check("uy", 1e-10, 5e-6)
--newtonCheck:set_component_check("p", 1e-10, 5e-6)

newtonSolver = NewtonSolver()
newtonSolver:set_linear_solver(lsolver)
newtonSolver:set_convergence_check(newtonCheck)
newtonSolver:set_line_search(lineSearch)
--newtonSolver:set_debug(dbgWriter)

nlsolver = newtonSolver

print("Interpolation start values")

if (dim==2) then
function Pres0(x, y, t) if x*x+y*y < 4 then return 1.0 else return 0.0 end end
function VelX0(x, y, t) return 0.0 end
function VelY0(x, y, t) return 0.0 end

Interpolate("Pres0", u, "p", startTime)
Interpolate("VelX0", u, "ux", startTime)
Interpolate("VelY0", u, "uy", startTime)
end

if (dim==3) then
function Pres0(x, y, z, t)
 if ((x-0.5)*(x-0.5)+y*y+(z+0.5)*(z+0.5)) < 0.0625 
 	then return 1.0  
 elseif ((x+0.5)*(x+0.5)+(y-0.5)*(y-0.5)+(z+0.5)*(z+0.5)) < 0.0625 
 	then return 1.0  
 elseif ((x+0.5)*(x+0.5)+(y+0.5)*(y+0.5)+(z+0.5)*(z+0.5)) < 0.0625 
 	then return 1.0 
 else return 0.0 end
 end 
 
 
function VelX0(x, y, z, t) return 0.0 end
function VelY0(x, y, z, t) return 0.0 end
function VelZ0(x, y, z, t) return 0.0 end

Interpolate("Pres0", u, "p", startTime)
Interpolate("VelX0", u, "ux", startTime)
Interpolate("VelY0", u, "uy", startTime)
Interpolate("VelZ0", u, "uz", startTime)
end


--util.SolveLinearTimeProblem(u, domainDisc, lsolver, vtk, "MandelTransient",
--							   "ImplEuler", 1, startTime, endTime, dt, dtmin, dtred);
							   
--util.SolveNonlinearTimeProblem(u, domainDisc, nlsolver, vtk, "PoroElasticityTransient",
--						   "ImplEuler", 1, startTime, endTime, dt, dtMin, dtRed); 
						   

						   
util.SolveNonlinearProblemAdaptiveTimestep(u, domainDisc, nlsolver, vtk, "PoroElasticityAdaptive",
	startTime, endTime, dt, dtMin, dtMax, dtRed, timeTol)
						    
						   
						   

exit();
