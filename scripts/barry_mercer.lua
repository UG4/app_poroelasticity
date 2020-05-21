

--[[
 Barry-Mercer problem as reported by Phillips, UT Austin, 2005
--]]

BARRY_MERCER_DATA = {
  NAPPROX = 512, --hires=64
  X0 = 0.25,
  Y0 = 0.25,
}

BarryMercerProblem = {}

function BarryMercerFourierCoeff_P(n, q, t_norm)
  local beta = BARRY_MERCER_DATA.BETA 
  local x0=BARRY_MERCER_DATA.X0
  local y0=BARRY_MERCER_DATA.Y0
  
  if (n%4==0) or (q%4==0) then return 0.0 end
  
  local lambda_n = n*math.pi
  local lambda_q = q*math.pi
  local _lambda_nq = lambda_n*lambda_n + lambda_q*lambda_q
   
  local val1 = -2.0 * math.sin(lambda_n*x0) * math.sin(lambda_q*y0)
  local val2 = (_lambda_nq*math.sin(t_norm) - math.cos(t_norm) + math.exp(-_lambda_nq*t_norm))
  local val3 = (1 + _lambda_nq*_lambda_nq)
  -- s[[
  --print (t_norm)
  -- print ("val2: (n="..n..", q="..q.."):".. _lambda_nq*math.sin(t_norm) .."-"..math.cos(t_norm).."+".. math.exp(-_lambda_nq*t_norm))
 --  print ("coeff: (n="..n..", q="..q..", m="..n+q.."):".. val1 .."*"..val2.."/"..val3.."="..(val1*val2)/val3 )
  --]]
  return  (val1*val2)/val3
end


function BarryMercerPressure2D(x, y, t)

  local N = BARRY_MERCER_DATA.NAPPROX;
  local t_norm = BARRY_MERCER_DATA.BETA*t
  
  local beta = BARRY_MERCER_DATA.BETA 
  local x0=BARRY_MERCER_DATA.X0
  local y0=BARRY_MERCER_DATA.Y0
 
   local sinbt = math.sin(t_norm) 
   local cosbt = math.cos(t_norm)
  
  
  local pp = 0.0
--  local sumCoeff2 = 0.0
local m=0
    for m=2,N do
      for k=1,m-1 do
  --     for k=1,N do
        local n = k
        local sinnx = math.sin(math.pi*n*x)
        local lambda_n = n*math.pi
 --   for q=1,N do
       local q = m-k
       local lambda_q = q*math.pi
       local sinqy = math.sin(math.pi*q*y)
 
       local _coeff_nq = BarryMercerFourierCoeff_P(n,q, t_norm) 
       -- print ("x= ("..x..", "..y.."\tn="..n.."q="..q..",m="..m.."), c=".._coeff_nq.."\t"..sinnx.."\t"..sinqy.."=> \tpp="..pp.."\tupdate=".._coeff_nq * sinnx * sinqy)
        pp = pp + _coeff_nq * sinnx * sinqy
   end
   -- print("break: "..sinnx)
  end
  -- print ("x= ("..t..","..t_norm..","..x..","..y..")"..pp)
  return -4.0*(BARRY_MERCER_DATA.LAMBDA+2.0*BARRY_MERCER_DATA.MU)*pp;
end

function BarryMercerVelX2D(x, y, t)

  local N = BARRY_MERCER_DATA.NAPPROX;
  local t_norm = BARRY_MERCER_DATA.BETA*t
 
  local u = 0.0
  local m=0
    for m=2,N do
      for k=1,m-1 do
  --for k=1,N do
  
    local n=k
    local _lambda_n = math.pi*n
    local _cosnx = math.cos(_lambda_n*x)*_lambda_n
    
    -- for q=1,N do
      local q = m-k
        local _lambda_nq = (math.pi*math.pi)*(n*n+ q*q)
        local _coeff_nq = BarryMercerFourierCoeff_P(n,q, t_norm) 
     -- < print ("x= ("..x..","..y..","..n..","..q.."), coeff="..(_coeff_nq/_lambda_nq)..", u="..u)
        u = u + _coeff_nq * _cosnx * math.sin(math.pi*q*y)/_lambda_nq
    end
  end
  
  return 4.0*u;
end

function BarryMercerVelY2D(x, y, t)

  local N = BARRY_MERCER_DATA.NAPPROX;
  local t_norm = BARRY_MERCER_DATA.BETA*t
 
  local u = 0.0
   local m=0
    for m=2,N do
      for k=1,m-1 do
  -- for q=1,N do
  local q = m-k
    local cosqx = math.cos(math.pi*q*y)*math.pi*q
    -- for n=1,N do
      local n=k
        local _lambda_nq = (math.pi*math.pi)*(n*n+ q*q)
        local _coeff_nq = BarryMercerFourierCoeff_P(n,q, t_norm) 
      --  print ("x= ("..x..","..y..")"..pp..", c=".._coeff_nq)
        u = u + _coeff_nq * math.sin(math.pi*n*x)*cosqx/_lambda_nq
    end
  end
  
  return 4.0*u;
end



-- This defines the point singularity (Note: time should be selected according to char_time!)
function BarryMercerSource2D(x,y,t,si)
  local delta = BARRY_MERCER_DATA.DELTA
  local beta = BARRY_MERCER_DATA.BETA
  local t_hat = beta*t
  -- NOTE: t_hat= beta*t
   if ( math.abs(x-0.25) < 0.5*delta and  math.abs(y-0.25) < 0.5*delta) then 
      return 2.0*beta*math.sin(t_hat)/(delta*delta)   
   else 
     return 0.0
   end
end


-- This defines the point singularity (Note: time should be selected according to char_time!)
function BarryMercerDiracSource2D(x,y,t,si)

  local beta = BARRY_MERCER_DATA.BETA
  local t_hat = beta*t  -- rescaling time to [0,2*PI]
   
  return 2.0*beta*math.sin(t_hat)
   
end

-- problem definition
barrymercer2D_tri = {
 
 gridName = "../grids/barrymercer2D-tri.ugx",
 -- gridName = "../grids/barrymercer2D.ugx",
  dim = 2,
  cpu = 1,
  
  porder = 1,
  uorder = 2,
  -- vStab = 1.0/12.0, -- *h^2 * \triangle p 
  
  mandatorySubsets ={"INNER", "HORIZONTAL", "VERTICAL", "CORNERS"},
  
  modelParameter = { },
  elemDiscParams = { },

  bRAP = true,

}


-- Initialize all variables
function barrymercer2D_tri:init(kperm, nporo, nu, cmedium, cfluid, csolid, volumetricweight)
  print ("WARNING: Ignoring all parameters")
  
  
  local E = 1e+5        -- Young's elasticity modulus [Pa]
  local nu = 0.4        -- Poisson"s ratio  [1]
  local kappa = 1e-5   --  permeability [m*m]  
  local muf = 1e-3       -- Pa*s    => Diff Coeff 1e-9
  local alpha = 1.0
  local Kcomp = E/(3*(1-2*nu))  -- compression (or bulk) modulus)
  
  local Kdim = {}
  local Kv = 2.0*E/(1+nu)*(1.0-nu)/(1.0-2.0*nu)                                -- uni-axial drained bulk modulus 
  
  Kdim[1] = Kv
  Kdim[2] = Kv/(2.0-2.0*nu)
  Kdim[3] = Kcomp
  
  self.elemDiscParams[1] = { 
     VOLUME = "INNER",
     KAPPA = kappa/muf, 
     LAMBDA=(E*nu)/((1.0+nu)*(1.0-2.0*nu)), 
     MU = 0.5*E/(1+nu), 
     ALPHA=alpha, 
     PHI= 0, 
     THETA=(alpha*alpha)/Kdim[self.dim] }
  
  print ("theta_stab= "..self.elemDiscParams[1].THETA)
  
  BARRY_MERCER_DATA.KAPPA =  self.elemDiscParams[1].KAPPA
  BARRY_MERCER_DATA.LAMBDA =  self.elemDiscParams[1].LAMBDA
  BARRY_MERCER_DATA.MU =  self.elemDiscParams[1].MU
 
  local beta = self.elemDiscParams[1].KAPPA*(self.elemDiscParams[1].LAMBDA + 2*self.elemDiscParams[1].MU)
  self.elemDiscParams[1].BETA =  beta
  BARRY_MERCER_DATA.BETA = beta
   print ("beta= "..beta)
  
  -- Reset stabilization
  if (self.vStab and self.vStab ~= 0.0) then
      print ("vStab(preReset)="..self.vStab)
      local lambda = self.elemDiscParams[1].LAMBDA
      local G = self.elemDiscParams[1].MU
      local factor = lambda+2.0*G
      self.vStab = self.vStab/factor
      print ("vStab(postReset)="..self.vStab)
  end
  
end

-- Read parameters from command line.
function barrymercer2D_tri:parse_cmd_args()

  local numRefs      = util.GetParamNumber("--num-refs", 3, "total number of refinements (incl. pre-Refinements)")
  BARRY_MERCER_DATA.DELTA = 0.25*math.pow(0.5,numRefs)
  print("BARRY_MERCER_DATA.GRID_SIZE="..BARRY_MERCER_DATA.DELTA)
  
  -- Number of terms in approximation.
  BARRY_MERCER_DATA.NAPPROX =  util.GetParamNumber("--bm-napprox", 256, "frequency bound (for analytic solution) ") 
  
  
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


function barrymercer2D_tri:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end


function barrymercer2D_tri:add_elem_discs(domainDisc, bStationary)

  -- Add standard element discs
  CommonAddBiotElemDiscs(self, domainDisc, bStationary, self.uorder, self.porder)
  
  -- Add a singular source 
  -- self.flowDisc[1]:set_source("BarryMercerSource2D")
  --local pointSourceDisc = ConvectionDiffusionFV1("p", "INNER")
  --pointSourceDisc:set_source("BarryMercerSource2D")
  
  local pointSourceDisc = DiracSourceDisc("p", "SINGULARITY")
  pointSourceDisc:add_source("BarryMercerDiracSource2D", Vec2d(0.25, 0.25))
  domainDisc:add(pointSourceDisc)
  
  
  
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


function barrymercer2D_tri:add_uzawa_discs(domainDisc)
  CommonAddMassMatrixDiscs(self, domainDisc)
end

--- Boundary conditions.
function barrymercer2D_tri:add_boundary_conditions(domainDisc, bStationary)
 
 local doStationary = bStationary or false

  -- Drainage zone.
  -- local dirichlet = DirichletBoundary(false, true)
  local dirichlet = DirichletBoundary(false)
  dirichlet:add(0.0, "p", "VERTICAL,HORIZONTAL,CORNERS")
  dirichlet:add(0.0, "ux", "HORIZONTAL,CORNERS")
  dirichlet:add(0.0, "uy", "VERTICAL,CORNERS")
  
  domainDisc:add(dirichlet)
 
end


function barrymercer2D_tri:get_char_time()
  local consolidation = self.elemDiscParams[1].KAPPA* (self.elemDiscParams[1].LAMBDA + 2*self.elemDiscParams[1].MU)
  local tchar = (1.0*1.0)/consolidation*math.pi
    print ("Characteristic time: " ..  tchar)
  return tchar -- seconds
end

-- Initial values.
function barrymercer2D_tri:interpolate_start_values(u, startTime)
  Interpolate(0.0, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
end


-- Create error estimator.
function barrymercer2D_tri:error_estimator()

  local p = self.elemDiscParams[1]
  local gamma2 = (p.LAMBDA+2*p.MU)*(p.LAMBDA+2*p.MU)
  print("barrymercer2D_tri:error_estimator: gamma^2="..gamma2)
  
  local biotErrorEst = CompositeGridFunctionEstimator()
  
  biotErrorEst:add(L2ComponentSpace("p", 2))       
  biotErrorEst:add(H1SemiComponentSpace("ux", 4, gamma2, p.VOLUME))  
  biotErrorEst:add(H1SemiComponentSpace("uy", 4, gamma2, p.VOLUME))
  
  return biotErrorEst
end



-- Computes various norms
function ComputeNorms(uref)

local normDesc={}

local porder = 2
local uorder = 4

normDesc["l2norm-p"] = L2Norm(uref, "p", porder)
normDesc["l2norm-ux"] = L2Norm(uref, "ux", uorder)
normDesc["l2norm-uy"] = L2Norm(uref, "uy", uorder)

normDesc["h1semi-ux"] = H1SemiNorm(uref, "ux", uorder)
normDesc["h1semi-uy"] = H1SemiNorm(uref, "uy", uorder)

return normDesc

end

function PrintNorms(normDesc)
for key, val in pairs(normDesc) do  print(key.."\t"..val) end
end

function CompareNorms(normDesc, errDesc, refDesc)
for key, val in pairs(normDesc) do  print(key.."\t"..val.."\t"..errDesc[key].."\t("..errDesc[key]/refDesc[key]..")") end
end


-- Post processing (after each step)
function barrymercer2D_tri:post_processing(u, step, time)
 
 local napprox = BARRY_MERCER_DATA.NAPPROX
 if (napprox<=0) then return end
  
  
  local normRef = {} -- norm of reference
  local normSol = {} -- norm of solution
  local normErr = {} -- norm of error
  
  local vtk = VTKOutput()
  local uref= u:clone()
 --  local norm = L2Norm(u, "p", 2)
 -- local err = L2Error("BarryMercerPressure2D", u, "p", time, 2)
  print ("NAPROX ="..napprox)
  
  
  uref:set(0.0)
  Interpolate("BarryMercerPressure2D", uref, "p", "INNER,SINGULARITY,CORNERS,HORIZONTAL,VERTICAL", time) -- "SINGULARITY"
  Interpolate("BarryMercerVelX2D", uref, "ux", "INNER,SINGULARITY,CORNERS,HORIZONTAL,VERTICAL", time) -- "SINGULARITY"
  Interpolate("BarryMercerVelY2D", uref, "uy", "INNER,SINGULARITY,CORNERS,HORIZONTAL,VERTICAL", time) -- "SINGULARITY"
 -- Interpolate("BarryMercerVelX2D", uref, "ux") -- "SINGULARITY"
 -- Interpolate("BarryMercerVelY2D", uref, "uy") -- "SINGULARITY"
  vtk:print("BarryMercer2D_Ref.vtu", uref, step, time)
  normSol = ComputeNorms(u)
  normRef = ComputeNorms(uref)
  print ("REFERENCE:")
  PrintNorms(normRef)

  
  VecScaleAdd2(uref, 1.0, uref, -1.0, u)
  vtk:print("BarryMercer2D_Err.vtu", uref, step, time)
  normErr = ComputeNorms(uref)
  print ("SOLUTION/ERROR:")
  CompareNorms(normSol, normErr, normRef)
  
  -- grep'able output
  local charTime =self:get_char_time()
  
  print ("deltaP:\t"..time.."\t"..time/charTime.."\t"..
    normErr["l2norm-p"].."\t"..normSol["l2norm-p"].."\t"..normRef["l2norm-p"])
    
  print ("deltaU1A:\t"..time.."\t"..time/charTime.."\t"..
    normErr["h1semi-ux"].."\t"..normSol["h1semi-ux"].."\t"..normRef["h1semi-ux"])
    
  print ("deltaU2A:\t"..time.."\t"..time/charTime.."\t"..
    normErr["h1semi-uy"].."\t"..normSol["h1semi-uy"].."\t"..normRef["h1semi-uy"])
    
  print ("deltaU1B:\t"..time.."\t"..time/charTime.."\t"..
    normErr["l2norm-ux"].."\t"..normSol["l2norm-ux"].."\t"..normRef["l2norm-ux"])
 
  print ("deltaU2B:\t"..time.."\t"..time/charTime.."\t"..
    normErr["l2norm-uy"].."\t"..normSol["l2norm-uy"].."\t"..normRef["l2norm-uy"])

end


-- This is a check evaluating the accuracy of approximation to the analytical solution
function barrymercer2D_tri:check(u)

u:set(0)

local time = self:get_char_time()*0.25 -- Evaluate for peak at PI/2
local vtk = VTKOutput()

local uref = u:clone()
local error = u:clone()
local normRef = {}
local normErr = {}
local normU ={}

local porder = 2
local uorder = 4


-- Most accurate solution
local origOrder = BARRY_MERCER_DATA.NAPPROX -- e.g., 256-4 for 3 numrefs 

if (origOrder<=0) then return end

local listOfTimes = { 
  self:get_char_time()*0.01,
  self:get_char_time()*0.25
}

for key,time in pairs(listOfTimes) do

print("TIME="..time)
BARRY_MERCER_DATA.NAPPROX = origOrder

local kmax = 4

local filenameSol = "BarryMercerCheck2DSol_".. BARRY_MERCER_DATA.NAPPROX .. "key"..key..".vtu"
local filenameErr = "BarryMercerCheck2DErr_".. BARRY_MERCER_DATA.NAPPROX .. "key"..key..".vtu"

-- Reference solution.
Interpolate("BarryMercerPressure2D", uref, "p", "INNER,SINGULARITY", time) -- "SINGULARITY"
Interpolate("BarryMercerVelX2D", uref, "ux", "INNER,SINGULARITY", time) -- "SINGULARITY"
Interpolate("BarryMercerVelY2D", uref, "uy", "INNER,SINGULARITY", time) -- "SINGULARITY"
vtk:print(filenameSol, uref, kmax, time)

normRef = ComputeNorms(uref)
print ("Reference for k="..kmax.." using m='"..BARRY_MERCER_DATA.NAPPROX)
-- print (kmax..": norm(uref,p)="..normRef["p"]..","..normRef["ux"]..","..normRef["uy"])
PrintNorms(normRef)

-- Comparison with less accurate terms.
for k=kmax-1,0,-1 do


BARRY_MERCER_DATA.NAPPROX = BARRY_MERCER_DATA.NAPPROX/2

print ("Check for k="..k.." using m='"..BARRY_MERCER_DATA.NAPPROX)

Interpolate("BarryMercerPressure2D", u, "p", "INNER,SINGULARITY", time) -- "SINGULARITY"
Interpolate("BarryMercerVelX2D", u, "ux", "INNER,SINGULARITY", time) -- "SINGULARITY"
Interpolate("BarryMercerVelY2D", u, "uy", "INNER,SINGULARITY", time) -- "SINGULARITY"
vtk:print(filenameSol, u, k, time)
normU = ComputeNorms(u)


VecScaleAdd2(error, 1.0, uref, -1.0, u)
vtk:print(filenameErr, error, k, time)

normErr = ComputeNorms(error)
CompareNorms(normU, normErr, normRef) -- Evaluate

end

end -- listOfTImes
-- restore
BARRY_MERCER_DATA.NAPPROX = origOrder 
-- quit() -- debugging
end

