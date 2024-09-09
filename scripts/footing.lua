-- Sample problem from N. Castelleto et al., J Comp Phys 327 (2016) 894-918
print ("Loading 'footing2D/3D'...")

FOOTING_DEFAULT_VALUES = {
-- Castelletto et al, JCP, 2016
  EYOUNG = 1e+6, -- Pa
  NU = 0.2,
  
  ALPHA  = 1.0,
  KAPPA  = 1e-9, -- m^2/(Pa s)
  
  FORCE = 200 -- N/m
}

FOOTING_DEFAULT_VALUES.LAMBDA = (FOOTING_DEFAULT_VALUES.NU*FOOTING_DEFAULT_VALUES.EYOUNG)/((1.0+FOOTING_DEFAULT_VALUES.NU)*(1.0-2.0*FOOTING_DEFAULT_VALUES.NU))
FOOTING_DEFAULT_VALUES.MU = 0.5*FOOTING_DEFAULT_VALUES.EYOUNG/(1.0+FOOTING_DEFAULT_VALUES.NU)

footing2D = {
 
  gridName = "../grids/footing2D_checker.ugx",
  dim = 2,
  cpu = 1,
  
  porder = 1, uorder = 2, vStab = 0.0, 
--  porder = 1, uorder = 1, vStab = 1.0/12.0, 
  
  mandatorySubsets ={"INNER", "IMPERMEABLE", "DRAINAGE", "FOOT"},
  
  modelParameter = { },
  elemDiscParams = { },

  bRAP = true,
  bAdjustTransfers = true,

}

function GenericParseCmdArgs(self)
  self.bRAP    = self.bRAP or util.HasParamOption("--use-rap", "Using Galerkin product")
  self.porder  = util.GetParamNumber("--orderP", 1, "Order for pressure") 
  self.uorder  = util.GetParamNumber("--orderU", 2, "Order for displacement") 
  self.vStab   = util.GetParamNumber("--stab", 0.0, "Stabilization") 
  
  if (self.bRAP) then 
    print ("bAdjustTransfers=false") 
    self.bAdjustTransfers = true
  else 
    print ("bAdjustTransfers=true") 
    self.bAdjustTransfers = true
  end
  
  print ("pOrder="..self.porder)
  print ("uOrder="..self.uorder)
  print ("vStab0="..self.vStab)
end

function GenericAddStabilization(self, domainDisc)
  -- CommonAddBiotStabDiscs(self, domainDisc) 
  -- print ("vStab1="..self.vStab) 
  if (self.vStab and self.porder==self.uorder) then 
  
    local stab = self.vStab or 0.0
    print ("stab="..stab)
    
    self.stabDisc = {}  
     for i=1,#self.elemDiscParams do
        local _parami = self.elemDiscParams[i]
        local _gammai = (_parami.LAMBDA+2*_parami.MU)
         print ("GenericAddStabilization: "..stab/_gammai.." on ".._parami["VOLUME"])
         self.stabDisc[i] = ConvectionDiffusionStabFE("p", _parami["VOLUME"], stab/_gammai)
         domainDisc:add(self.stabDisc[i])
     end -- for
  end --if
end

function GenericFootingInit(self)

  local E = 1e+6        -- Young's elasticity modulus [Pa]
  local E2 = 1e+8 
  local nu = 0.2        -- Poisson"s ratio  [1]
  local kappa = 1e-12   --  permeability [m*m]  
  local mu = 1e-3       -- Pa*s    => Diff Coeff 1e-9
  local alpha = 1.0
  local Kcomp = E/(3*(1-2*nu))  -- compression (or bulk) modulus)
  
  local Kdim = {}
  local Kv = 2.0*E/(1+nu)*(1.0-nu)/(1.0-2.0*nu)                                -- uni-axial drained bulk modulus 
  
  local _LAME_LAMBDA=E*nu/((1.0+nu)*(1.0-2.0*nu))
  local _LAME_MU=0.5*E/(1+nu)
  
  local _LAME_LAMBDA2=E2*nu/((1.0+nu)*(1.0-2.0*nu))
  local _LAME_MU2=0.5*E2/(1+nu)
   
  Kdim[1] = Kv
  Kdim[2] = Kv/(2.0-2.0*nu)
  Kdim[3] = Kcomp*2.0 -- 2.0 = Wheeler& Mikelic convergence
  
  self.elemDiscParams[1] = { 
     VOLUME = "INNER",
     KAPPA = kappa/mu, 
     LAMBDA= _LAME_LAMBDA, --E*nu/((1.0+nu)*(1.0-2.0*nu)), 
     MU = _LAME_MU, -- 0.5*E/(1+nu), 
     ALPHA=alpha, 
     PHI= 0.0, 
    --  THETA=(alpha*alpha)/Kdim[self.dim] 
     THETA=0.5*(alpha*alpha)/(2.0*_LAME_MU/self.dim + _LAME_LAMBDA)  -- Wheeler & Mikelic
  }
  
  self.elemDiscParams[2] = { 
     VOLUME = "INNER2",
     KAPPA = kappa/mu, 
     LAMBDA= _LAME_LAMBDA2, --E*nu/((1.0+nu)*(1.0-2.0*nu)), 
     MU = _LAME_MU2, -- 0.5*E/(1+nu), 
     ALPHA=alpha, 
     PHI= 0.0, 
     THETA=0.5*(alpha*alpha)/(2.0*_LAME_MU2/self.dim + _LAME_LAMBDA2)  -- Wheeler & Mikelic
  }
  
  print ("beta_stab= "..self.elemDiscParams[1].THETA)
  print ("beta_stab= "..self.elemDiscParams[2].THETA)
end


function GenericFootingCharTime(self)
  local tchar = 0.0
 for i=1,#self.elemDiscParams do
  local consolidation  = self.elemDiscParams[i].KAPPA* (self.elemDiscParams[i].LAMBDA + 2*self.elemDiscParams[i].MU)
  print ("Characteristic time: " ..  (1.0*1.0)/consolidation)
  tchar = math.max(tchar, (1.0*1.0)/consolidation)
 end -- for

 return tchar  -- seconds
end
-- Read parameters from command line.
function footing2D:parse_cmd_args()

  self.bRAP    = self.bRAP or util.HasParamOption("--use-rap", "Using Galerkin product")
  self.porder  = util.GetParamNumber("--orderP", self.porder, "Order for pressure") 
  self.uorder  = util.GetParamNumber("--orderU", self.uorder, "Order for displacement") 
  self.vStab   = util.GetParamNumber("--stab", self.vStab, "Stabilization") 
  
  if (self.bRAP) then 
    print ("bAdjustTransfers=false") 
    self.bAdjustTransfers = true
  else 
    print ("bAdjustTransfers=true") 
    self.bAdjustTransfers = true
  end
  
  print ("pOrder="..self.porder)
  print ("uOrder="..self.uorder)
  print ("vStab0="..self.vStab)
  
end

function footing2D:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end


function footing2D:add_elem_discs(domainDisc, bStationary)

  -- Add standard element discs
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
  
  -- Add stabilization.
  GenericAddStabilization(self, domainDisc)
end


function footing2D:add_uzawa_discs(domainDisc)
  CommonAddMassMatrixDiscs(self, domainDisc)
end

--- Boundary conditions.
function footing2D:add_boundary_conditions(domainDisc, bStationary)
 
 local doStationary = bStationary or false

 local PaUnit=1.0

 -- Footing zone.
 local tz = 1e+3*PaUnit  -- 1kPa  (negative sign due to discretisation!)
 
 local neumannZ = NeumannBoundaryFE("uy")
 neumannZ:add(tz, "FOOT", "INNER")
 
 if doStationary then 
    neumannZ:set_stationary() 
 end
 domainDisc:add(neumannZ)

  -- Drainage zone.
  -- local dirichlet = DirichletBoundary(false,self.bAdjustTransfers)
  local dirichlet = DirichletBoundary()
  dirichlet:add(0.1, "p", "DRAINAGE")
  dirichlet:add(0.0, "ux", "IMPERMEABLE,CORNERS")
  dirichlet:add(0.0, "uy", "BOTTOM,CORNERS")
  
  domainDisc:add(dirichlet)
 
end

-- Initialize all variables
function footing2D:init(kperm, nporo, nu, cmedium, cfluid, csolid, volumetricweight)
  GenericFootingInit(self) 
end

function footing2D:get_char_time()
  return GenericFootingCharTime(self)
end

-- Initial values.
function footing2D:interpolate_start_values(u, startTime)
  Interpolate(0.1, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
end


-- Create error estimator.
function footing2D:error_estimator()
  local p = self.elemDiscParams[1]
  local gamma2 = (p.LAMBDA+2*p.MU)*(p.LAMBDA+2*p.MU)
  print("footing2D:error_estimator: gamma="..gamma2)
  
  local biotErrorEst = CompositeGridFunctionEstimator()
        
  biotErrorEst:add(H1SemiComponentSpace("ux", 4, gamma2, p.VOLUME))  
  biotErrorEst:add(H1SemiComponentSpace("uy", 4, gamma2, p.VOLUME))
  biotErrorEst:add(L2ComponentSpace("p", 2)) 
  
  return biotErrorEst
end


-- Post processing (after each step)
function footing2D:post_processing(u, step, time)
 -- print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..(Integral(u, "p", "CENTER")/self.modelParameter.p0).."\t"..Cryer.ComputePressure(self, time))
end


footing2D_tri = {
 
  gridName = "../grids/footing2D-tri.ugx",
  dim = 2,
  cpu = 1,
  
  porder = 1,
  uorder = 2,
  -- vStab = 1.0/12.0, -- *h^2 * \triangle p 
  
  mandatorySubsets ={"INNER", "IMPERMEABLE", "DRAINAGE", "FOOT"},
  
  modelParameter = { },
  elemDiscParams = { },

  bRAP = false,
  bAdjustTransfers = true,

}

-- Read parameters from command line.
function footing2D_tri:parse_cmd_args() 
  GenericParseCmdArgs(self)
end

function footing2D_tri:create_domain(numRefs, numPreRefs)
  return util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
end

function footing2D_tri:add_elem_discs(domainDisc, bStationary)
  -- Add standard element discs
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
  GenericAddStabilization(self, domainDisc)
end


function footing2D_tri:add_uzawa_discs(domainDisc)
  CommonAddMassMatrixDiscs(self, domainDisc)
end

--- Boundary conditions.
function footing2D_tri:add_boundary_conditions(domainDisc, bStationary)
 
 local doStationary = bStationary or false

 local PaUnit=1.0

 -- Footing zone.
 local tz = 1e+3*PaUnit  -- 1kPa  (negative sign due to discretisation!)
 
 local neumannZ = NeumannBoundaryFE("uy")
 neumannZ:add(tz, "FOOT", "INNER")
 
 if doStationary then 
    neumannZ:set_stationary() 
 end
 domainDisc:add(neumannZ)

  -- Drainage zone.
 -- local dirichlet = DirichletBoundary(false, true) --true = truncate (default) false, true
  local dirichlet = DirichletBoundary()
  dirichlet:add(0.1, "p", "DRAINAGE")
  dirichlet:add(0.0, "ux", "IMPERMEABLE,CORNERS")
  dirichlet:add(0.0, "uy", "BOTTOM,CORNERS")
  
  domainDisc:add(dirichlet)
 
end

-- Initialize all variables
function footing2D_tri:init(kperm, nporo, nu, cmedium, cfluid, csolid, volumetricweight)
  print ("WARNING: Ignoring all parameters") 
  GenericFootingInit(self)
end

function footing2D_tri:get_char_time() 
  return GenericFootingCharTime(self)
end

-- Initial values.
function footing2D_tri:interpolate_start_values(u, startTime)
  Interpolate(0.1, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
end


-- Create error estimator.
function footing2D_tri:error_estimator()
  local p = self.elemDiscParams[1]
  local gamma2 = (p.LAMBDA+2*p.MU)*(p.LAMBDA+2*p.MU)
  print("footing2D:error_estimator: gamma="..gamma2)
  
  local biotErrorEst = CompositeGridFunctionEstimator()
        
  biotErrorEst:add(H1SemiComponentSpace("ux", 4, gamma2, p.VOLUME))  
  biotErrorEst:add(H1SemiComponentSpace("uy", 4, gamma2, p.VOLUME))
  biotErrorEst:add(L2ComponentSpace("p", 2)) 
  
  return biotErrorEst
end


-- Post processing (after each step)
function footing2D_tri:post_processing(u, step, time)
 -- print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..(Integral(u, "p", "CENTER")/self.modelParameter.p0).."\t"..Cryer.ComputePressure(self, time))
end


--[[ 
3D PROBLEM
--]]

footing3D = {
 
  gridName = "../grids/footing3D-tet.ugx",
  -- gridName = "../grids/footing3D.ugx",
  
  dim = 3,
  cpu = 1,
  
  porder = 1,
  uorder = 1,
  mandatorySubsets ={"INNER", "IMPERMEABLE", "DRAINAGE", "FOOT", "EDGES"},
  
  modelParameter = { },
  elemDiscParams = { },
  
  bRAP = false,
  bAdjustTransfers = true,

}

-- Read parameters from command line.
function footing3D:parse_cmd_args()
  GenericParseCmdArgs(self)
end

function footing3D:create_domain(numRefs, numPreRefs)
  print("footing3D:create_domain")
  return util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
end


function footing3D:add_elem_discs(domainDisc, bStationary)
  print("footing3D:add_elem_discs")
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
  GenericAddStabilization(self, domainDisc)
end


function footing3D:add_uzawa_discs(domainDisc)
  print("footing3D:add_uzawa_discs")
  CommonAddMassMatrixDiscs(self, domainDisc)
end

--- Boundary conditions.
function footing3D:add_boundary_conditions(domainDisc, bStationary)
 print("footing3D:add_boundary_conditions")
 local doStationary = bStationary or false

 local PaUnit=1.0

 -- Footing zone.
 local tz = 1000*PaUnit  -- 1kPa  (negative sign due to discretisation!)
 
 local neumannZ = NeumannBoundaryFE("uz")
 neumannZ:add(tz, "FOOT", "INNER")
 
 if doStationary then 
    neumannZ:set_stationary() 
 end
 domainDisc:add(neumannZ)

  -- Drainage zone.
  --local dirichlet = DirichletBoundary(false, true)
  local dirichlet = DirichletBoundary()
  dirichlet:add(0.1, "p", "DRAINAGE")
  dirichlet:add(0.0, "ux", "IMPERMEABLE,EDGES,BOTTOM_CAGE")
  dirichlet:add(0.0, "uy", "IMPERMEABLE,EDGES,BOTTOM_CAGE")
  dirichlet:add(0.0, "uz", "BOTTOM,BOTTOM_CAGE")
  
  domainDisc:add(dirichlet)
 
end

-- Initialize all variables
function footing3D:init(kperm, nporo, nu, cmedium, cfluid, csolid, volumetricweight)
  print("footing3D:init")
  print ("WARNING: Ignoring all parameters")
  GenericFootingInit(self)
end

function footing3D:get_char_time()
  return GenericFootingCharTime(self)
end

-- Initial values.
function footing3D:interpolate_start_values(u, startTime)
  Interpolate(0.1, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
  Interpolate(0.0, u, "uz", startTime)
end

-- Post processing (after each step)
function footing3D:post_processing(u, step, time)
 -- print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..(Integral(u, "p", "CENTER")/self.modelParameter.p0).."\t"..Cryer.ComputePressure(self, time))
end

print ("    ... done!")