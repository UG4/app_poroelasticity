generic = {}

function CreateMassMatrixDiscs(param, dim)

  print ("CreateMassMatrixDiscs:"..param.VOLUME)
  
  local theta = param.THETA or 0.0   -- this is: alpha^2/K
  print ("theta= "..theta)  
  
  local massEqDisc = ConvectionDiffusion("p", param["VOLUME"], "fv1")
  massEqDisc:set_mass_scale(theta)
  
  return massEqDisc;
end


function CreatePoissonMatrixDiscs(param, dim)

  print ("CreatePoissonMatrixDiscs:"..param["VOLUME"])
  
  local schurEqDisc = ConvectionDiffusion("p", param["VOLUME"], "fv1")
  schurEqDisc:set_diffusion(param["THETA"])
  schurEqDisc:set_stationary()

  return schurEqDisc;  
end

--[[
  \partial_t (PHI * p + ALPHA div(u) ) + \nabla \cdot [-KAPPA \nabla p] = f
  \nabla \cdot [SIGMA-ALPHA*p*Id] = 0,        SIGMA  := func(LAMBDA, MU)
  
  
  param["LAMBDA"] 
  param["MU"]
  param["ALPHA"]
  param["PHI"]
  param["KAPPA"]
  
  returns pair of flow and displacement disc
--]]



function DEPRECATED_CreateBiotElemDiscs(param, dim, bSteadyState)

  local doSteadyState = bSteadyState or false
  -- define eqns for pressure,  displacement
  print ("CreateBiotElemDiscs: "..param["VOLUME"])
  
  -- elasticity
  print("=> const_lambda = "..param["LAMBDA"])
  print("=> const_mu     = "..param["MU"])
  
  -- quasi-static Bio
  print("=> const_alpha  = "..param["ALPHA"])
  print("=> const_Phi    = "..param["PHI"])
  print("=> const_kappa  = "..param["KAPPA"])
  
  
  local flowEqDisc = ConvectionDiffusion("p", param["VOLUME"], "fe")
  
  local displacementEqDisc
  if (dim==2) then displacementEqDisc = SmallStrainMechanics("ux,uy", param["VOLUME"]) end
  if (dim==3) then displacementEqDisc = SmallStrainMechanics("ux,uy,uz", param["VOLUME"]) end
 
  if doSteadyState then displacementEqDisc:set_stationary() end  -- do not scale with tau?
  
  -- specify displacement eq (for u)
  local matLaw = HookeLaw()
  matLaw:set_hooke_elasticity_tensor(param["LAMBDA"], param["MU"])  -- corresponds to plane strain in 2D
  displacementEqDisc:set_material_law(matLaw)
  displacementEqDisc:set_mass_scale(0.0)

  -- if (param["ALPHA"]>0) then
    if (false) then
      -- tested
      local forceLinker = ScaleAddLinkerVector()
      forceLinker:add(-param["ALPHA"], flowEqDisc:gradient())
      displacementEqDisc:set_volume_forces(forceLinker)
    else
      -- more natural
      local divLinker = ScaleAddLinkerNumber()
      divLinker:add(param["ALPHA"], flowEqDisc:value())
      displacementEqDisc:set_div_factor(divLinker)
    end


    -- specify flow eq (for pressure)
    local compressionLinker = ScaleAddLinkerNumber()
    compressionLinker:add(param["ALPHA"], displacementEqDisc:divergence())
    flowEqDisc:set_mass(compressionLinker);
  
 -- end -- alpha >0
  
  
  flowEqDisc:set_mass_scale(param["PHI"]);  -- Storativity 1.0/M = S ‰
  flowEqDisc:set_diffusion(param["KAPPA"]);


  -- adjust quadrature order
  if (dim==2) then
  
    if (porder==1) then flowEqDisc:set_quad_order(4) end
    if (uorder==2) then displacementEqDisc:set_quad_order(4) end
  --  displacementEqDisc:set_quad_order(5)


  elseif (dim == 3) then
    --if (order == 1) then
    --  displacementEqDisc:set_quad_order(2)  
    --3:#ip`s: 6, 2:#ip`s: 8 
    --end
    if (uorder == 2) then
    displacementEqDisc:set_quad_order(5)  -- 7
    flowEqDisc:set_quad_order(3) --7
    --#ip`s: 31
    end
    --if (order == 3) then  
    --  elemDisc:set_quad_order(11) 
    --#ip`s: 90
    --  end
    --if (order > 3) then 
    --  elemDisc:set_quad_order(11) 
      --#ip`s: 90
    --  end
  end
  
    -- print info
  --print(flowEqDisc:config_string())
  
  print(displacementEqDisc:config_string())

return flowEqDisc, displacementEqDisc 

end


-- Create elem discs for consistency
function DEPRECATED_CreateConsistencyElemDiscs(param, dim, bSteadyState)

  local doSteadyState = bSteadyState or false
  -- define eqns for pressure,  displacement
  print ("...for subset "..param["VOLUME"])
  
  -- elasticity
  print("=> const_lambda = "..param["LAMBDA"])
  print("=> const_mu     = "..param["MU"])
  
  
  print("=> const_alpha  = "..param["ALPHA"])
  print("=> const_Phi    = "..param["PHI"])
  print("=> const_kappa  = "..param["KAPPA"])
  
  
  local flowEqDisc = ConvectionDiffusion("p", param["VOLUME"], "fe")
  
  local displacementEqDisc
  if (dim==2) then displacementEqDisc = SmallStrainMechanics("ux,uy", param["VOLUME"]) end
  if (dim==3) then displacementEqDisc = SmallStrainMechanics("ux,uy,uz", param["VOLUME"]) end
 
  if doSteadyState then displacementEqDisc:set_stationary() end  -- do not scale with tau?
  
  -- specify displacement eq (for u)
  local matLaw = HookeLaw()
  matLaw:set_hooke_elasticity_tensor(param["LAMBDA"], param["MU"])  -- corresponds to plane strain in 2D
  displacementEqDisc:set_material_law(matLaw)
  displacementEqDisc:set_mass_scale(0.0)

    if (param["ALPHA"]>0) then
    if (false) then
      -- tested
      local forceLinker = ScaleAddLinkerVector()
      forceLinker:add(-param["ALPHA"], flowEqDisc:gradient())
      displacementEqDisc:set_volume_forces(forceLinker)
    else
      -- more natural
      local divLinker = ScaleAddLinkerNumber()
      divLinker:add(param["ALPHA"], flowEqDisc:value())
      displacementEqDisc:set_div_factor(divLinker)
    end


  end -- alpha >0
  
   -- specify flow eq (for pressure)
  flowEqDisc:set_mass_scale(param["PHI"]);  -- Storativity 1.0/M = S ‰
  flowEqDisc:set_diffusion(0.0);


  -- adjust quadrature order
  if (dim==2) then
  
 --   if (porder==1) then flowEqDisc:set_quad_order(4) end
  --  if (uorder==2) then displacementEqDisc:set_quad_order(4) end
  -- displacementEqDisc:set_quad_order(5)


  elseif (dim == 3) then
    if (uorder == 1) then
      displacementEqDisc:set_quad_order(2)  
    --3:#ip`s: 6, 2:#ip`s: 8 
    end
    if (uorder == 2) then
    displacementEqDisc:set_quad_order(5)  -- 7
    flowEqDisc:set_quad_order(3) --7
    --#ip`s: 31
    end
    --if (order == 3) then  
    --  elemDisc:set_quad_order(11) 
    --#ip`s: 90
    --  end
    --if (order > 3) then 
    --  elemDisc:set_quad_order(11) 
      --#ip`s: 90
    --  end
  end
  
    -- print info
  --print(flowEqDisc:config_string())
  
  print(displacementEqDisc:config_string())

return flowEqDisc, displacementEqDisc 

end


-- Creates discretization for stiffness matrix.
function CommonAddBiotElemDiscs(self, domainDisc,  bStationary)
  
  self.flowDisc = {}
  self.dispDisc = {}
  

  -- NEW C++ style
  local ucmps="ux,uy"
  if (dim==3) then ucmps = ucmps + ",uz" end
  local factory = BiotElemDiscFactory(ucmps, self.uorder, "p", self.porder, bStationary)
  
  -- Create element discretizations.
  for i=1,#self.elemDiscParams do
     local param = self.elemDiscParams[i]
    
    if (false) then
      -- Old LUA style (i.e. forward call to 'generic.lua')
      self.flowDisc[i], self.dispDisc[i] = DEPRECATED_CreateBiotElemDiscs(param, self.dim, bStationary)
    else
      -- New C++ style 
      -- TODO: -> JSON-style init
      local biotparams = BiotSubsetParameters(param["VOLUME"], param["ALPHA"], param["KAPPA"], param["PHI"], param["LAMBDA"], param["MU"], param["THETA"])
      local jstringTest = util.json.encode(param)
      print(jstringTest)

      local biot_disc = factory:create_elem_discs(biotparams) 
      self.dispDisc[i] = biot_disc:displacement_disc()
      self.flowDisc[i] = biot_disc:pressure_disc()
    end
   
    -- Add to domain disc.
    domainDisc:add(self.flowDisc[i])
    domainDisc:add(self.dispDisc[i])
    
  end
end

-- Creates discretization for stabilization matrix.
function CommonAddBiotStabDiscs(self, domainDisc, bStationary)
  
  self.stabDisc = {}
  
  local stab = self.vStab or 0.0;
  if (stab==0.0) then return end
  
  -- forward call to 'generic.lua'
  if (bStationary == false) then print ("ERROR: Only implemented for Mass matrix!") return quit; end;
     
  for i=1,#self.elemDiscParams do
    local _parami = self.elemDiscParams[i]
    self.stabDisc[i] = ConvectionDiffusionStabFE("p", _parami["VOLUME"], stab)
    domainDisc:add(self.stabDisc[i])
  end
end


-- Creates discretization for mass matrix.
function CommonAddMassMatrixDiscs(self, domainDisc)

  self.massDisc = {}

  -- use generic.lua
  for i=1,#self.elemDiscParams do
    self.massDisc[i] = CreateMassMatrixDiscs(self.elemDiscParams[i], self.dim)
    domainDisc:add(self.massDisc[i])
  end
end

-- Add discretizations for flow (and Schur complement!)
function AddElemDiscsP(self, domainDiscP)

  -- use generic.lua
  for i=1,#self.elemDiscParams do
    domainDiscP:add(self.flowDisc[i])
    domainDiscP:add(self.massDisc[i])
    if (self.stabDisc[i]) then domainDiscP:add(self.stabDisc[i]) end
  end
end

-- Add discretizations for deformations
function AddElemDiscsU(self, domainDiscU)

  -- use generic.lua
  for i=1,#self.elemDiscParams do
    domainDiscU:add(self.dispDisc[i])
  end
end
