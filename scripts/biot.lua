------------------------------------------------------------------------------
--
--   Lua - Script for bio-/geomechanics
--
--   Author: Arne Naegel
--           (based on solid_mechanics app by Raphael Prohl)
--
------------------------------------------------------------------------------

ug_load_script("ug_util.lua")

dim = 2

-----------------------------------------------------------------
--  Settings
-----------------------------------------------------------------

-- ORDER OF ANSATZ-FUNCTIONS 
porder = 1
uorder = 1 
InitUG(dim, AlgebraType("CPU", 1));

-- REFINEMENT
-- choose number of pre-Refinements (before sending grid onto different processes)      
numPreRefs = util.GetParamNumber("-numPreRefs", 1)
-- choose number of total Refinements (incl. pre-Refinements)
numRefs = util.GetParamNumber("-numRefs", 3)


startTime  = util.GetParamNumber("-start", 0.0, "start time")
endTime    = util.GetParamNumber("-end", 50000.0, "end time") 
dt 		   = util.GetParamNumber("-dt", 1.0, "time step size")
dtmin	   = util.GetParamNumber("-dtmin", 1.0, "minimal admissible time step size")
dtred	   = util.GetParamNumber("-dtred", 0.5, "time step size reduction factor on divergence")



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

gridName = util.GetParam("-grid", "../grids/mandel.ugx");
neededSubsets = {"INNER", "TOP", "BOTTOM", "LEFT", "RIGHT"};

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

-- find first n roots of tan x = alpha x
function findRoots(n, alpha, x1)



x=x1
for k=1,n,1 do 

--e=x  -- reset error
--F=100.0

F = math.tan(x) - alpha*x
while (math.abs(F)>1e-4) do
	y=x
	print(k..":"..x.." "..F)
--dF = (4.0*math.cos(x)*math.cos(x))/((math.cos(2.0*x)+1.0)*(math.cos(2.0*x)+1.0)) - alpha
	dF = 1.0/(math.cos(x)*math.cos(x)) - alpha
	s=-F/dF
	--y=math.atan(alpha*x)
    --e=math.abs(x-y)
	--print(k..":"..x.." "..y.." "..e.." "..math.tan(x)-alpha*x)
	
	Fold=F
	lambda=1.0
	repeat 
		lambda = 0.5*lambda
		x = y + lambda * s
		F = math.tan(x) - alpha*x
	until (F<Fold)
	print(" :"..x.." "..F)

end


-- determine next guess
x = x + math.pi
end

end

--findRoots(4, 4, 1.4);
--quit();
-----------------------------------------------------------------
--  Boundary Conditions & Rhs
-----------------------------------------------------------------

-- zero dirichlet boundary conditions
function PDirichletSide(x, y, t) return true, 0.0 end       -- p  (left, right)
function UxDirichletPlate(x, y, t) return true, 0.0 end   -- ux (bottom)
function UyDirichletPlate(x, y, t) return true, (-1.0)*y*math.exp(-0.01*t) end   -- uy (bottom)

--function ourVolumeForceField2d(x, y, t)
--	return 0.0, 0.0 
--end

--volumeForceField = LuaUserVector("ourVolumeForceField2d")

--dirichletZero = LuaCondUserNumber("DirichletBndZero" .. dim .. "d")
--neumannBndVal = LuaCondUserNumber("NeumannBnd" .. dim .. "d")

mandelDirichletBnd = DirichletBoundary()
mandelDirichletBnd:add(LuaCondUserNumber("PDirichletSide"), "p", "LEFT")
mandelDirichletBnd:add(LuaCondUserNumber("PDirichletSide"), "p", "RIGHT")
mandelDirichletBnd:add(LuaCondUserNumber("UxDirichletPlate"), "ux", "TOP")
mandelDirichletBnd:add(LuaCondUserNumber("UyDirichletPlate"), "uy", "TOP")
mandelDirichletBnd:add(LuaCondUserNumber("UxDirichletPlate"), "ux", "BOTTOM")
mandelDirichletBnd:add(LuaCondUserNumber("UyDirichletPlate"), "uy", "BOTTOM")

dirichletBnd = mandelDirichletBnd
--neumannBnd = NeumannBoundaryFE("ux")
--neumannBnd:add(neumannBndVal, "Top", "Inner")

--------------------------------------------------------------------------------
--------------------------------------------------------------------------------
--  Setup FE Linear Element Discretization
--------------------------------------------------------------------------------
--------------------------------------------------------------------------------

matLaw = HookeLaw()
E      = 5.94*1e+9    -- Young's modulus [Pa]
nu     = 0.2          -- Poissons's ratio 

alpha  = 1.0;         -- Biot's constant
M      = 1.65*10e+10;  -- Biot's modulus
--M=1.0

k      = 1          -- permeability
mu     = 1.0          -- fluid viscosity
Lambda = k/mu;

B      = 0.83333      -- Skempton coeff
rho0   = 1.0;


matLaw:set_hooke_elasticity_tensor_E_nu(E, nu)
--matLaw:set_hooke_elasticity_tensor(110743.82,80193.80)

-- define eqns for displacement, pressure
displacementEqDisc = SmallStrainMechanics("ux, uy", "INNER")
flowEqDisc = ConvectionDiffusion("p", "INNER", "fe")

-- specify eqn (for displacement)
forceLinker = ScaleAddLinkerVector()
forceLinker:add(alpha, flowEqDisc:gradient())

displacementEqDisc:set_material_law(matLaw)
--displacementEqDisc:set_mass_scale(0.0)
displacementEqDisc:set_volume_forces(forceLinker)
--mechOut = MechOutputWriter()
--displacementEqDisc:set_output_writer(mechOut)
--displacementEqDisc:displacement()


-- specify flow eqn (for pressure)
massLinker = ScaleAddLinkerNumber()
massLinker:add(alpha, displacementEqDisc:divergence())

flowEqDisc:set_mass_scale(1.0/M);
flowEqDisc:set_mass(massLinker);
flowEqDisc:set_diffusion(Lambda);

-- print info
print(displacementEqDisc:config_string())
--print(flowEqDisc:config_string())
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
	gmg:set_base_level(1)
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

lsolver = luSolver
--solver = jacSolver
--solver = iluSolver
--solver = bicgstabSolver 
--lsolver = gmgSolver
--solver = cgSolver


u = GridFunction(approxSpace)

vtk=VTKOutput()
vtk:select_nodal("p", "PNodal")
vtk:select({"ux", "uy"}, "DispNodal")
vtk:select_element( displacementEqDisc:displacement(), "DispElem")
vtk:select_element( displacementEqDisc:divergence(), "DivElem")

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

SaveMatrixForConnectionViewer(u, A, "Stiffness.mat")

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
newtonSolver:set_debug(dbgWriter)

nlsolver = newtonSolver

print("Interpolation start values")

function Pres0(x, y, t) return 0.0 end
function VelX0(x, y, t) return 0.0 end
function VelY0(x, y, t) return 0.0 end


Interpolate("Pres0", u, "p", startTime)
Interpolate("VelX0", u, "ux", startTime)
Interpolate("VelY0", u, "uy", startTime)

--util.SolveLinearTimeProblem(u, domainDisc, lsolver, vtk, "MandelTransient",
--							   "ImplEuler", 1, startTime, endTime, dt, dtmin, dtred);
							   
util.SolveNonlinearTimeProblem(u, domainDisc, nlsolver, vtk, "MandelTransient",
						   "ImplEuler", 1, startTime, endTime, dt, dtmin, dtred);  

exit();
