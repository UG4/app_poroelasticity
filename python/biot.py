
import sys
sys.path.append("/home/anaegel/Software/ug4-git/")
print(sys.path)

print("Imports:")
import ug4py.pyugcore as ug4
import ug4py.pylimex as limex
import ug4py.pyconvectiondiffusion as cd
import ug4py.pysmallstrainmechanics as mech
import ug4py.pyporoelasticity as biot
import ug4py.pyutil as util


print("Modules loaded:")
print(ug4)
print(util)
print(limex)
print(cd) 
print(mech)
print(biot)


# Configure disc.
config = biot.BiotDiscConfig("ux,uy",1, "p", 1, 1.0/12.0)
problem = biot.BarryMercerProblem2dCPU1(config)
print(problem)



# Create domain, two refinements.
dom=util.CreateDomain2d(problem.get_gridname(), 0)

# Approximation space.
approxSpace = ug4.ApproximationSpace2d(dom)
approxSpace.add_fct(config.pcmp(), "Lagrange", config.porder())
approxSpace.add_fct("ux", "Lagrange", config.uorder())
approxSpace.add_fct("uy", "Lagrange", config.uorder())
approxSpace.init_levels()
approxSpace.init_surfaces()
approxSpace.init_top_surface()
approxSpace.print_statistic()
approxSpace.print_layout_statistic()
approxSpace.print_local_dof_statistic(2)  



# Domain disc(s).
bSteadyStateMechanics = True
domainDisc = ug4.DomainDiscretization2dCPU1(approxSpace)
problem.add_elem_discs(domainDisc, bSteadyStateMechanics) 
problem.add_boundary_conditions(domainDisc, bSteadyStateMechanics)  

uzawaSchurUpdateDisc = ug4.DomainDiscretization2dCPU1(approxSpace)
problem.add_uzawa_discs(uzawaSchurUpdateDisc, bSteadyStateMechanics)



# Time stuff

# Create a LIMEX time disc.
limexDefaultDesc = {
        "nstages": 2,
        "lsolver" : ug4.LUCPU1(), 
        "TOL": 1e-3,
        "metricSpace": None
}

# Get default LIMEX time disc.
def GetLimexDefaultDesc():
    return limexDefaultDesc

# Create a LIMEX time disc.
def CreateLimexIntegrator(domainDisc, limexDesc=limexDefaultDesc, dt=1.0, dtmin=None, dtmax=None):

    nstages = limexDesc["nstages"] 
    lsolver = limexDesc["lsolver"]
    TOL     = limexDesc["TOL"]  
    metricSpace = limexDesc["metricSpace"]

    if (dtmin is None):
        dtmin=dt*1e-3

    if (dtmax is None):
        dtmax=dt*1e+2   

    if (nstages<2):
        print("Using implicit Euler (nstages=1).")
        return None
    
    # LIMEX config.
    timeInt = limex.LimexTimeIntegrator2dCPU1(nstages)
    nlsolver = limex.LimexNewtonSolverCPU1()
    nlsolver.set_linear_solver(lsolver)
    for i in range(nstages):
        timeInt.add_stage(i+1, nlsolver, domainDisc)

    # Time stepping config.
    timeInt.set_time_step(dt)
    timeInt.set_dt_min(dtmin)
    timeInt.set_dt_max(dtmax)
    timeInt.set_increase_factor(1.5) # max. Faktor, um den dt erhöht werden darf
    timeInt.disable_matrix_cache()

    # LIMEX w/ error estimation
    timeInt.set_tolerance(TOL)

    # Definition of error estimator.
    errorEst = None
    if metricSpace is not None:
        print("Using custom metric space for error estimation.")
        errorEst = limex.CompositeGridFunctionEstimator2dCPU1()
        errorEst.add(metricSpace)
    else:
        print("Using  euclidean norm for error estimation.")
        errorEst = limex.Norm2EstimatorCPU1() # Euclidean norm.
    
    timeInt.add_error_estimator(errorEst)
    return timeInt
    


# Time data.
startTime = problem.start_time()
charTime = problem.get_char_time()
endTime = 2.0*charTime

# Grid function(s)
print("Interpolation start values")
u = ug4.GridFunction2dCPU1(approxSpace)
problem.interpolate_start_values(u, startTime)
ug4.Interpolate(0.0, u, "p")


# Test solver for robustness
lsolver = ug4.LUCPU1()
timeDisc=ug4.ThetaTimeStepCPU1(domainDisc, 1.0) 
timeInt=limex.ConstStepLinearTimeIntegrator2dCPU1(timeDisc)
timeInt.set_linear_solver(lsolver)

dt = 1e-3*charTime
timeInt.set_time_step(dt)

# timeInt.apply(u, dt, u, 0.0)


# Solve time dependent problem.
limexDesc=GetLimexDefaultDesc()
limex = CreateLimexIntegrator(domainDisc, limexDesc, charTime*1e-3, charTime*1e-8, charTime/2)
limex.set_dt_min(charTime*1e-8)
limex.set_dt_max(charTime/2)



print("Solve problem")
limex.apply(u, endTime, u, startTime)



