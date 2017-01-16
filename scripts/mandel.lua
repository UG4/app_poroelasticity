
-- A new inheritsFrom() function
--
function inheritsFrom( baseClass )

    local new_class = {}
    local class_mt = { __index = new_class }

    function new_class:create()
        local newinst = {}
        setmetatable( newinst, class_mt )
        return newinst
    end

    if nil ~= baseClass then
        setmetatable( new_class, { __index = baseClass } )
    end

    -- Implementation of additional OO properties starts here --

    -- Return the class object of the instance
    function new_class:class()
        return new_class
    end

    -- Return the super class object of the instance
    function new_class:superClass()
        return baseClass
    end

    -- Return true if the caller is an instance of theClass
    function new_class:isa( theClass )
        local b_isa = false

        local cur_class = new_class

        while ( nil ~= cur_class ) and ( false == b_isa ) do
            if cur_class == theClass then
                b_isa = true
            else
                cur_class = cur_class:superClass()
            end
        end

        return b_isa
    end

    return new_class
end

---------------------------
-- problem definition
---------------------------
mandel2d = {
  
  gridName = "../grids/mandel-aniso.ugx",
  dim = 2,
  cpu = 1,
  
  porder = 1,
  uorder = 2,
  mandatorySubsets ={"INNER"},
  
  a = 100.0,
  b = 10.0,
  qForce = 1.0,
  
  n_approx = 200,
  
  modelParameter = {},
  elemDiscParams = {}
  
}


mandel3d = {
  
  gridName = "../grids/mandel3d-aniso.ugx",
  dim = 3,
  cpu = 1,
  
  porder = 1,
  uorder = 1,
  mandatorySubsets ={"INNER"},
  
  a = 100.0,
  b = 10.0,
  qForce = 1.0,
  
  n_approx = 200,
  
  modelParameter = {},
  elemDiscParams = {}
  
}



local function MandelSetModelParameters_(self, nu_poisson, nporo, cmedium, cfluid, csolid, D)

  -- lower case: non dimensional
  
  local nu = nu_poisson or 0.2
  local poro = nporo or 0.0
  
  local CMedium =  cmedium or 1.0  
  local CFluid =  cfluid or 0.0           -- default: incompressible
  local CSolid =  csolid or 0.0           -- default: incompressible => alpha =1.0

  local cscm = CSolid/CMedium;
  local cfcm = CFluid/CMedium;
   
  local alpha = 1.0-cscm;                -- non-dimensional Biot coeff

  self.modelParameter.alpha = alpha;
  self.modelParameter.poro = poro; 
  self.modelParameter.S = poro * CFluid + (alpha-poro) * CSolid;         -- Storativity S
 
 
  self.modelParameter.nu0 = nu;                                          -- Poisson's ratio
  self.modelParameter.K = 1.0/CMedium;                                   -- compression (bulk) modulus of the medium
  local lambda = (3.0*self.modelParameter.K*nu)/(1.0+nu)                 -- 1st Lame parameter
  local mu = (3.0*self.modelParameter.K*(1.0-2.0*nu))/(2.0+2.0*nu)       -- 2nd Lame parameter (also: shear modulus G)
  self.modelParameter.mu = mu
   
  self.modelParameter.cscm = cscm;
  self.modelParameter.cfcm = cfcm;

  local KS_0  = poro * cfcm + (alpha-poro)*cscm;
  self.modelParameter.KS_0  = KS_0;
    
  local KSG_0 = (3.0*KS_0)/(2.0+2.0*nu);    -- 
  self.modelParameter.KSG_0 = KSG_0;
  
  -- non-dimensional parameter (required for roots)
  local eta = (1.0-nu)*(alpha*alpha+KSG_0)/(alpha*alpha*(1.0-2.0*nu));
  self.modelParameter.eta = eta; 
  self.modelParameter.nu1 = nu+(1.0-nu)/(2.0*eta);  
  
  
  self.elemDiscParams = {}
  self.elemDiscParams[1] =
  {VOLUME = "INNER",  KAPPA = D, LAMBDA=lambda, MU = mu, ALPHA=alpha, PHI=self.modelParameter.S}
  
  print ("Phi=S="..self.modelParameter.S)
  print ("Alpha="..alpha)
  print ("Kappa="..D)
  
end

-- return consolidation coefficient (from modelParameter)
local function MandelGetCharTime_(self)

  local var =  (self.modelParameter.K + 4.0*self.modelParameter.mu/3.0)
  local cnu = (self.elemDiscParams[1].KAPPA) * var;
  cnu = cnu/(self.modelParameter.alpha*self.modelParameter.alpha + var*self.modelParameter.S)
  return cnu/(self.a*self.a);
end

--Verruijt: Poroelasticity (lecture notes)
local function MandelInitRoots_(self)

  self.ROOT = {}
  local N = self.n_approx;
  local eps=0.000001;
  
  local PI=math.pi;
  local eta = self.modelParameter.eta;

  for i=1,N do

    local a1=(i-1)*PI+PI/4;
    local a2=a1+PI/2-eps;
  local am;
    for j=1,40 do
      local y1=math.tan(a1)-2*eta*a1;
      local y2=math.tan(a2)-2*eta*a2;
      am=(a1+a2)/2;
      local ym=math.tan(am)-2*eta*am;
      if (ym*y1>0) then a1=am else a2=am end
      if (math.abs(y2)<eps) then am=a2 end
    end
  self.ROOT[i]=am
  --print (am);
  end
end -- function 'initRoots_'


local function MandelAddElemDiscs(self, domainDisc, bStationary)
  
  self.flowDisc = {}
  self.dispDisc = {}
  
 for i=1,#self.elemDiscParams do
    self.flowDisc[i], self.dispDisc[i] = CreateElemDiscs(self.elemDiscParams[i], self.dim, bStationary)
    domainDisc:add(self.flowDisc[i])
    domainDisc:add(self.dispDisc[i])
  end
end


-- compute normalized pore pressure
local function MandelComputeP(self, xhat, yhat, t)

 local eta = self.modelParameter.eta; 
 local N = self.n_approx;
 
 local pp0 = 0.0;
 for i=1,N do
    local si  = math.sin(self.ROOT[i]);
    local ci  = math.cos(self.ROOT[i]);
    local cxi = math.cos(xhat*self.ROOT[i]);
    local aai  = 2.0*eta*self.ROOT[i]*si-(2.0*eta-1.0)*ci
    local tti  = self.ROOT[i]*self.ROOT[i]*t;
    local coeffi =  2.0*eta/aai
     -- add to sum (abort, if contribution is small)
    if (tti<20) then  pp0 = pp0 + 2.0*eta*(cxi-ci)*math.exp(-tti)/aai -- 2.0*eta/aai = Coeff_i(Gai)
    else i=N+1 end
  end
 
 return pp0;

end -- function 'computeP'


local function MandelComputeUx(self, x, y, t)
  local xhat = x/self.a
  local yhat = x/self.b
  
  local eta = self.modelParameter.eta;
  local N = self.n_approx; 
  
  local nu0 = self.modelParameter.nu0;
  local nu1 = self.modelParameter.nu1;
  local F = self.qForce;
  local G = self.modelParameter.mu;

  local ux0 = 0.5*(F*nu0)/(self.a*G);-- 0.5(Fv/Ga)
  local ux1 = 0.0
  
  local fac0 =  (F*nu1)/(self.a*G)  -- (F*vu/Ga)
  local fac1 =  F/G; -- (F/G)

  for i=1,N do
    
    local si  = math.sin(self.ROOT[i]);
    local ci  = math.cos(self.ROOT[i]);
    local cxi = math.cos(xhat*self.ROOT[i]);
    local sxi = math.sin(xhat*self.ROOT[i]);
    
    local aai  = 2.0*eta*self.ROOT[i]*si-(2.0*eta-1.0)*ci
    local tti  = self.ROOT[i]*self.ROOT[i]*t;
    local coeffi = 2.0*eta/aai -- 2.0*eta/aai = Coeff_i(Gai)
    
    if (tti<20) then 
      ux0 = ux0 + fac0*coeffi*ci*math.exp(-tti)  
      ux1 = ux1 + fac1*coeffi*math.exp(-tti)*(2.0*eta*self.ROOT[i])*sxi --
    else i=N+1 end 
      
  end
  
 -- print(ux0.."*"..xhat.."+"..ux1.."="..ux0*xhat+ux1)
  return ux0*xhat+ux1
end

local function MandelComputeUy(self, x, y, t)
 
 local yhat = y/self.b
 local eta = self.modelParameter.eta; 
 local N = self.n_approx;
 
 local nu0 = self.modelParameter.nu0;
 local nu1 = self.modelParameter.nu1;
 local F = self.qForce;
 local Ga = self.modelParameter.mu*self.a;
 
 
 local uy = -0.5*F*(1.0-nu0)/Ga; --  
 local fac = F*(1.0-nu1)/Ga;
 
 for i=1,N do
    local si  = math.sin(self.ROOT[i]);
    local ci  = math.cos(self.ROOT[i]);
    local aai  = 2.0*eta*self.ROOT[i]*si-(2.0*eta-1.0)*ci
    local tti  = self.ROOT[i]*self.ROOT[i]*t;
    local coeffi = 2.0*eta/aai  -- 2.0*eta/aai = Coeff_i(Gai)
    
    -- add to sum (abort, if contribution is small)
    if (tti<20) then uy = uy + fac*coeffi*ci*math.exp(-tti)  
    else i=N+1 end 
  end
 
 return uy*y;
end -- function 'compute_uy'

function mandel2d:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end

function mandel2d:set_load(Fy)
  self.qForce = Fy
end

function mandel2d:init(kperm, poro, nu, cmedium, cfluid, csolid, volumetricweight)
  volumetricweight = volumetricweight or 1.0 -- kg/m^3
  MandelSetModelParameters_(self, nu, poro, cmedium, cfluid, csolid, kperm*volumetricweight)
  MandelInitRoots_(self)
end

-- return consolidation coefficient
function mandel2d:get_char_time()
  return MandelGetCharTime_(self)
end

function mandel2d:compute_p(xhat, yhat, t)
  return MandelComputeP(self, xhat, yhat, t)
end -- function 'compute_p'

function mandel2d:compute_ux(x, y, t)
  return MandelComputeUx(self, x, y, t)
end

function mandel2d:compute_uy(x, y, t)
 return MandelComputeUy(self, x, y, t)
end 


function Mandel2dPressure(x,y,t)
  return mandel2d:compute_p(x,y,t)
end

function Mandel2dHorizDisp(x,y,t)
  return mandel2d:compute_ux(x,y,t)
end

function Mandel2dVertDisp(x,y,t)
  return mandel2d:compute_uy(x,y,t)
end

function mandel2d:add_boundary_conditions(domainDisc, bStationary)

  local neumann = NeumannBoundaryFE("ux")
  neumann:add(Mandel2dPressure, "TOP", "INNER")
  neumann:add(Mandel2dPressure, "BOT", "INNER")
  --neumann:add(self.modelParameter.alpha*self.flowDisc[1]:value(), "BOT", "INNER")
  domainDisc:add(neumann)
  
 
  local dirichlet = DirichletBoundary()
  dirichlet:add(Mandel2dHorizDisp, "ux", "TOP, BOT")
  dirichlet:add(Mandel2dVertDisp, "uy", "TOP, BOT")
  dirichlet:add(0.0, "p", "RIGHT")
  dirichlet:add(0.0, "p", "LEFT")
  domainDisc:add(dirichlet)
end

function mandel2d:add_elem_discs(domainDisc, bStationary)
  MandelAddElemDiscs(self, domainDisc, bStationary)
end

function mandel2d:interpolate_start_values(u, startTime)
  Interpolate(0.5*(self.qForce/self.a), u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
end

-- post processing (after each step)
function mandel2d:post_processing(u, step, time)
  local testResults = {}
  for k=0,20 do
    testResults[k] = self:compute_p(k/20, 0.0, time)
    print ((self.a*k/20) .."\t"..testResults[k])
  end

--  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..Integral(u, "p", "CENTER").."\t"..self:compute_pressure(0.0,0.0,time))
end

function mandel2d:create_test_data(time)
local testResults = {}
  print ("-----------------------")
  print ("time="..time)
  print ("-----------------------")
  testResults[time] = {}
  local k;  
  for k=0,20 do
    testResults[time][k] = self:compute_p(k/20, 0.0, time)
    print ((self.a*k/20) .."\t"..testResults[time][k])
  end
  print ("-----------------------")  
end






------------
--
--
--------

function mandel3d:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end

function mandel3d:set_load(Fy)
  self.qForce = Fy
end

function mandel3d:init(kperm, poro, nu, cmedium, cfluid, csolid, volumetricweight)
  volumetricweight = volumetricweight or 1.0 -- kg/m^3
  MandelSetModelParameters_(self, nu, poro, cmedium, cfluid, csolid, kperm*volumetricweight)
  MandelInitRoots_(self)
end

-- return consolidation coefficient
function mandel3d:get_char_time()
  return MandelGetCharTime_(self)
end

function mandel3d:compute_p(xhat, yhat, t)
  return MandelComputeP(self, xhat, yhat, t)
end -- function 'compute_p'

function mandel3d:compute_ux(x, y, t)
  return MandelComputeUx(self, x, y, t)
end

function mandel3d:compute_uy(x, y, t)
 return MandelComputeUy(self, x, y, t)
end 

function Mandel3dPressure(x,y,z,t)
  return mandel3d:compute_p(x,y,t)
end

function Mandel3dHorizDisp(x,y,z,t)
  return mandel3d:compute_ux(x,y,t)
end

function Mandel3dVertDisp(x,y,z,t)
  return mandel3d:compute_uy(x,y,t)
end

function mandel3d:add_boundary_conditions(domainDisc, bStationary)

  local neumann = NeumannBoundaryFE("ux")
  neumann:add(Mandel3dPressure, "TOP", "INNER")
  neumann:add(Mandel3dPressure, "BOT", "INNER")
  --neumann:add(self.modelParameter.alpha*self.flowDisc[1]:value(), "BOT", "INNER")
  domainDisc:add(neumann)
  
 
  local dirichlet = DirichletBoundary()
  dirichlet:add(Mandel3dHorizDisp, "ux", "TOP, BOT")
  dirichlet:add(0.0, "uy", "TOP, BOT")
  dirichlet:add(0.0, "uy", "FRONT, BACK")
  dirichlet:add(Mandel3dVertDisp, "uz", "TOP, BOT")
  dirichlet:add(0.0, "p", "RIGHT, LEFT")
  domainDisc:add(dirichlet)
end

function mandel3d:add_elem_discs(domainDisc, bStationary)
  MandelAddElemDiscs(self, domainDisc, bStationary)
end

function mandel3d:interpolate_start_values(u, startTime)
  Interpolate(0.5*(self.qForce/self.a), u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
  Interpolate(0.0, u, "uz", startTime)
end

-- post processing (after each step)
function mandel3d:post_processing(u, step, time)
  local testResults = {}
  for k=0,20 do
    testResults[k] = self:compute_p(k/20, 0.0, time)
    print ((self.a*k/20) .."\t"..testResults[k])
  end

--  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..Integral(u, "p", "CENTER").."\t"..self:compute_pressure(0.0,0.0,time))
end

function mandel3d:create_test_data(time)
local testResults = {}
  print ("-----------------------")
  print ("time="..time)
  print ("-----------------------")
  testResults[time] = {}
  local k;  
  for k=0,20 do
    testResults[time][k] = self:compute_p(k/20, 0.0, time)
    print ((self.a*k/20) .."\t"..testResults[time][k])
  end
  print ("-----------------------")  
end





CompositeProblem ={
  gridName = "../grids/cube_composite.ugx",
  mandatorySubsets = {"INNER", "INNER2", "TOP", "BOTTOM", "LEFT", "RIGHT"},
}


-- problem definitions
TwoMaterialProblem = {
  {VOLUME = "INNER", MU = 1.0, LAMBDA=1.0, KAPPA = 1.0, PHI =0.0},
  {VOLUME = "INNER2", MU = 1.0, LAMBDA=1.0, KAPPA = epsilon, PHI =0.0}, 
}

--[[
Lambda =1.0;
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
  function Dirichlet1(x, y, z, t) return true, 1.0 end  
  function DirichletNeg1(x, y, z, t) return true, -1.0 end  
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
--]]

ThreeRegionProblem = {
}


ConfinedCompression ={
}


function ConfinedCompression:add_dirichlet_boundary(domainDisc)

local ConfinedCompressionBnd = DirichletBoundary()

ConfinedCompressionBnd:add(0.0, "p", "TOP")
ConfinedCompressionBnd:add(-1.0, "uz", "TOP") 

ConfinedCompressionBnd:add(0.0, "ux", "FRONT")
ConfinedCompressionBnd:add(0.0, "uy", "FRONT")

ConfinedCompressionBnd:add(0.0, "ux", "BACK")
ConfinedCompressionBnd:add(0.0, "uy", "BACK")

ConfinedCompressionBnd:add(0.0, "ux", "RIGHT")
ConfinedCompressionBnd:add(0.0, "uy", "RIGHT")

ConfinedCompressionBnd:add(0.0, "ux", "LEFT")
ConfinedCompressionBnd:add(0.0, "uy", "LEFT")

ConfinedCompressionBnd:add(0.0, "ux", "BOTTOM")
ConfinedCompressionBnd:add(0.0, "uy", "BOTTOM")
ConfinedCompressionBnd:add(0.0, "uz", "BOTTOM")
domainDisc:add(ConfinedCompressionBnd);
end

