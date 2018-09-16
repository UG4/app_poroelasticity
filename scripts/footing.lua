-- Sample problem from N. Castelleto et al., J Comp Phys 327 (2016) 894-918
print ("Loading 'footing2D/3D'...")

footing2D = {
 
  gridName = "../grids/footing2D.ugx",
  dim = 2,
  cpu = 1,
  
  porder = 1,
  uorder = 2,
 -- vStab = 0.0, 
  
  mandatorySubsets ={"INNER", "IMPERMEABLE", "DRAINAGE", "FOOT"},
  
  modelParameter = { },
  elemDiscParams = { },

  bRAP = false,

}

-- Read parameters from command line.
function footing2D:parse_cmd_args()
end

function footing2D:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end


function footing2D:add_elem_discs(domainDisc, bStationary)

  -- Add standard element discs
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
  
  -- Add stabilization.
  if (self.vStab and self.porder==self.uorder) then 
    CommonAddBiotStabDiscs(self, domainDisc)   
  end
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
  local dirichlet = DirichletBoundary(false,true)
  dirichlet:add(0.1, "p", "DRAINAGE")
  dirichlet:add(0.0, "ux", "IMPERMEABLE,CORNERS")
  dirichlet:add(0.0, "uy", "BOTTOM,CORNERS")
  
  domainDisc:add(dirichlet)
 
end

-- Initialize all variables
function footing2D:init(kperm, nporo, nu, cmedium, cfluid, csolid, volumetricweight)
  print ("WARNING: Ignoring all parameters")
  
  
  local E = 1e+6        -- Young's elasticity modulus [Pa]
  local nu = 0.2        -- Poisson"s ratio  [1]
  local kappa = 1e-12   --  permeability [m*m]  
  local mu = 1e-3       -- Pa*s    => Diff Coeff 1e-9
  local alpha = 1.0
  local Kcomp = E/(3*(1-2*nu))  -- compression (or bulk) modulus)
  
  local Kdim = {}
  local Kv = 2.0*E/(1+nu)*(1.0-nu)/(1.0-2.0*nu)                                -- uni-axial drained bulk modulus 
  
  Kdim[1] = Kv
  Kdim[2] = Kv/(2.0-2.0*nu)
  Kdim[3] = Kcomp
  
  self.elemDiscParams[1] = { 
    VOLUME = "INNER",
     KAPPA = kappa/mu, 
     LAMBDA=E*nu/(1.0+nu)*(1.0-2.0*nu), 
     MU = E/(1+nu), 
     ALPHA=alpha, 
     PHI= 0, 
     THETA=(alpha*alpha)/Kdim[self.dim] }
  
  print ("theta_stab= "..self.elemDiscParams[1].THETA)
  
  
end

function footing2D:get_char_time()
local consolidation = self.elemDiscParams[1].KAPPA* (self.elemDiscParams[1].LAMBDA + 2*self.elemDiscParams[1].MU)
 print ("Characteristic time: " ..  (1.0*1.0)/consolidation)
 return (1.0*1.0)/consolidation -- seconds
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

}

-- Read parameters from command line.
function footing2D_tri:parse_cmd_args()
end

function footing2D_tri:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end


function footing2D_tri:add_elem_discs(domainDisc, bStationary)

  -- Add standard element discs
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
  
  -- Add stabilization.
  -- CommonAddBiotStabDiscs(self, domainDisc)   
  if (self.vStab and self.porder==self.uorder) then 
  
    local stab = self.vStab or 0.0
     self.stabDisc = {}
     
     for i=1,#self.elemDiscParams do
        local _parami = self.elemDiscParams[i]
        local _gammai = (_parami.LAMBDA+2*_parami.MU)
        self.stabDisc[i] = ConvectionDiffusionStabFE("p", _parami["VOLUME"], stab/_gammai)
        domainDisc:add(self.stabDisc[i])
     end -- for
  end --if
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
  local dirichlet = DirichletBoundary(false, true)
  dirichlet:add(0.1, "p", "DRAINAGE")
  dirichlet:add(0.0, "ux", "IMPERMEABLE,CORNERS")
  dirichlet:add(0.0, "uy", "BOTTOM,CORNERS")
  
  domainDisc:add(dirichlet)
 
end

-- Initialize all variables
function footing2D_tri:init(kperm, nporo, nu, cmedium, cfluid, csolid, volumetricweight)
  print ("WARNING: Ignoring all parameters")
  
  
  local E = 1e+6        -- Young's elasticity modulus [Pa]
  local nu = 0.2        -- Poisson"s ratio  [1]
  local kappa = 1e-12   --  permeability [m*m]  
  local mu = 1e-3       -- Pa*s    => Diff Coeff 1e-9
  local alpha = 1.0
  local Kcomp = E/(3*(1-2*nu))  -- compression (or bulk) modulus)
  
  local Kdim = {}
  local Kv = 2.0*E/(1+nu)*(1.0-nu)/(1.0-2.0*nu)                                -- uni-axial drained bulk modulus 
  
  Kdim[1] = Kv
  Kdim[2] = Kv/(2.0-2.0*nu)
  Kdim[3] = Kcomp
  
  self.elemDiscParams[1] = { 
    VOLUME = "INNER",
     KAPPA = kappa/mu, 
     LAMBDA=E*nu/(1.0+nu)*(1.0-2.0*nu), 
     MU = E/(1+nu), 
     ALPHA=alpha, 
     PHI= 0, 
     THETA=(alpha*alpha)/Kdim[self.dim] }
  
  print ("theta_stab= "..self.elemDiscParams[1].THETA)
  
  
end

function footing2D_tri:get_char_time()
local consolidation = self.elemDiscParams[1].KAPPA* (self.elemDiscParams[1].LAMBDA + 2*self.elemDiscParams[1].MU)
 print ("Characteristic time: " ..  (1.0*1.0)/consolidation)
 return (1.0*1.0)/consolidation -- seconds
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
 
  gridName = "../grids/footing3D.ugx",
  dim = 3,
  cpu = 1,
  
  porder = 1,
  uorder = 1,
  mandatorySubsets ={"INNER", "IMPERMEABLE", "DRAINAGE", "FOOT", "EDGES"},
  
  modelParameter = { },
  elemDiscParams = { },
  
  bRAP = false

}

-- Read parameters from command line.
function footing3D:parse_cmd_args()
end

function footing3D:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end


function footing3D:add_elem_discs(domainDisc, bStationary)
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
end


function footing3D:add_uzawa_discs(domainDisc)
  CommonAddMassMatrixDiscs(self, domainDisc)
end

--- Boundary conditions.
function footing3D:add_boundary_conditions(domainDisc, bStationary)
 
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
  local dirichlet = DirichletBoundary()
  dirichlet:add(0.1, "p", "DRAINAGE")
  dirichlet:add(0.0, "ux", "IMPERMEABLE,EDGES")
  dirichlet:add(0.0, "uy", "IMPERMEABLE,EDGES")
  dirichlet:add(0.0, "uz", "BOTTOM,EDGES")
  
  domainDisc:add(dirichlet)
 
end

-- Initialize all variables
function footing3D:init(kperm, nporo, nu, cmedium, cfluid, csolid, volumetricweight)
  print ("WARNING: Ignoring all parameters")
  
  
  local E = 1e+6        -- Young's elasticity modulus [Pa]
  local nu = 0.2        -- Poisson"s ratio  [1]
  local kappa = 1e-12   --  permeability [m*m]  
  local mu = 1e-3       -- Pa*s    => Diff Coeff 1e-9
  local alpha = 1.0
  local Kcomp = E/(3*(1-2*nu))  -- compression (or bulk) modulus)
  
  local Kdim = {}
  local Kv = 2.0*E/(1+nu)*(1.0-nu)/(1.0-2.0*nu)                                -- uni-axial drained bulk modulus 
  
  Kdim[1] = Kv
  Kdim[2] = Kv/(2.0-2.0*nu)
  Kdim[3] = Kcomp
  
  self.elemDiscParams[1] = { 
    VOLUME = "INNER",
     KAPPA = kappa/mu, 
     LAMBDA=E*nu/(1.0+nu)*(1.0-2.0*nu), 
     MU = E/(1+nu), 
     ALPHA=alpha, 
     PHI= 0, 
     THETA=(alpha*alpha)/Kdim[self.dim] }
  
  print ("theta_stab= "..self.elemDiscParams[1].THETA)
  
  
end

function footing3D:get_char_time()
 local tchar = self.elemDiscParams[1].KAPPA* (self.elemDiscParams[1].LAMBDA + 2*self.elemDiscParams[1].MU)
 print ("Characteristic time: " .. tchar)
 return tchar -- seconds
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