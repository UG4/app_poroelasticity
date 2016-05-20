function CreateElemDiscs(param, dim, bSteadyState)

--param["ALPHA"] = 0.0  --debugging
  local doSteadyState = bSteadyState or false
  -- define eqns for pressure,  displacement
  print (" "..param["VOLUME"])
  
  local flowEqDisc = ConvectionDiffusion("p", param["VOLUME"], "fe")
  
  local displacementEqDisc
  if (dim==2) then displacementEqDisc = SmallStrainMechanics("ux,uy", param["VOLUME"]) end
  if (dim==3) then displacementEqDisc = SmallStrainMechanics("ux,uy,uz", param["VOLUME"]) end
 

  if doSteadyState then displacementEqDisc:set_stationary() end  -- do not scale with tau?

  
  -- specify eqn (for displacement)
  local matLaw = HookeLaw()
  matLaw:set_hooke_elasticity_tensor(param["LAMBDA"], param["MU"])  -- corresponds to plane strain in 2D
  displacementEqDisc:set_material_law(matLaw)
  displacementEqDisc:set_mass_scale(0.0)

  local forceLinker = ScaleAddLinkerVector()
  forceLinker:add(-param["ALPHA"], flowEqDisc:gradient())
  displacementEqDisc:set_volume_forces(forceLinker)
 
--[[ 
  local presLinker = ScaleAddLinkerNumber()
  presLinker:add(param["ALPHA"], flowEqDisc:value()) 
  displacementEqDisc:set_div_factor(presLinker)
--]]

  --mechOut = MechOutputWriter()
  --displacementEqDisc:set_output_writer(mechOut)
  --displacementEqDisc:displacement()

  -- specify flow eqn (for pressure)
  local compressionLinker = ScaleAddLinkerNumber()
  compressionLinker:add(param["ALPHA"], displacementEqDisc:divergence())
 
  flowEqDisc:set_mass(compressionLinker);
  flowEqDisc:set_mass_scale(param["PHI"]); -- 1.0/M = S 
  flowEqDisc:set_diffusion(param["KAPPA"]);

  if (porder==1) then flowEqDisc:set_quad_order(2) end
  if (uorder==1) then displacementEqDisc:set_quad_order(2) end

  
  -- print info




  if dim == 3 then
    --if (order == 1) then
    --  displacementEqDisc:set_quad_order(2)  
    --3:#ip`s: 6, 2:#ip`s: 8 
    --end
    if (uorder == 2) then
    displacementEqDisc:set_quad_order(7)  -- 7
    flowEqDisc:set_quad_order(7) --7
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
  
  --print(flowEqDisc:config_string())
  print(displacementEqDisc:config_string())

  --massLinker = ScaleAddLinkerNumber()
  --massLinker:add(rho/M, flowEqDisc:value())
  --massLinker:add(rho*alpha, displacementEqDisc:divergence())

return flowEqDisc, displacementEqDisc 

end