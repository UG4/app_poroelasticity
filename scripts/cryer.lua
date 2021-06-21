
-- Neumann boundary (2D) 
-- opposite sign since we discretise
-- div[-sigma + alpha p*Id] = -g


function Sphere2dNormalX(x,y,t)
  return 1.0*x/math.sqrt(x*x+y*y)
end

function Sphere2dNormalY(x,y,t)
  return 1.0*y/math.sqrt(x*x+y*y)
end


-- Neumann boundary (3D)
function Cylinder3dNormalX(x,y,z,t)
  return 1.0*x/math.sqrt(x*x+y*y)
end

function Cylinder3dNormalY(x,y,z,t)
  return 1.0*y/math.sqrt(x*x+y*y)
end


-- Neumann boundary (3D)
function Sphere3dNormalX(x,y,z,t)
  return 1.0*x/math.sqrt(x*x+y*y+z*z)
end

function Sphere3dNormalY(x,y,z,t)
  return 1.0*y/math.sqrt(x*x+y*y+z*z)
end

function Sphere3dNormalZ(x,y,z,t)
  return 1.0*z/math.sqrt(x*x+y*y+z*z)
end




--[[ NOTE: Cryer in 2D should be de Leeuw! --]]
cryer2d = {

  gridName = "../grids/cryer2d-with-projector.ugx",
  dim = 2,
  cpu = 1,
  
  porder = 1,
  uorder = 2, -- should be 2
  mandatorySubsets ={"INNER"},
  
  a = 0.5, 
  qForce = 1.0,
  
  n_approx = 50,
  
  modelParameter = {},
  elemDiscParams = {}

}




local DeLeeuw = {}


-- Verruijt, pp. 77
function DeLeeuw.J0(u)

  local a,t,f,x,xx,c1,c2,c3,c4,c5,c6,c7;
  x=math.abs(u);
  
  if (x<3) then
  
    xx=x*x/9.0;
    c1=-2.2499997;
    c2=1.2656208; 
    c3=-0.3163866;
    c4=0.0444479;
    c5=-0.0039444;
    c6=0.0002100;
    f=1+xx*(c1+xx*(c2+xx*(c3+xx*(c4+xx*(c5+xx*c6)))));
  else
    xx=3.0/x;
    c1=0.79788456;
    c2=-0.00000077;
    c3=-0.00552740;
    c4=-0.00009512;
    c5=0.00137237;
    c6=-0.00072805;
    c7=0.00014476;
    a=c1+xx*(c2+xx*(c3+xx*(c4+xx*(c5+xx*(c6+xx*c7)))));
    
    c1=-0.78539816;
    c2=-0.04166397;
    c3=-0.00003954;
    c4=0.00262573;
    c5=-0.00054125;
    c6=-0.00029333;
    c7=0.00013558;
    t=x+c1+xx*(c2+xx*(c3+xx*(c4+xx*(c5+xx*(c6+xx*c7)))));
    
    f=a*math.cos(t)/math.sqrt(x);
  end
  return f;
end 

function DeLeeuw.J1(u)

  local a,t,f,x,xx,c1,c2,c3,c4,c5,c6,c7;
  x= math.abs(u);
  if (x<3) then
    xx=x*x/9.0
    c1=-0.56249985;
    c2=0.21093573;
    c3=-0.03954289;
    c4=0.00443319;
    c5=-0.00031761;
    c6=0.00001109;
    a=0.5+xx*(c1+xx*(c2+xx*(c3+xx*(c4+xx*(c5+xx*c6)))));
    f=a*x;
  else
    xx=3.0/x;
    c1=0.79788456;
    c2=0.00000156;
    c3=0.01659667;
    c4=0.00017105;
    c5=-0.00249511;
    c6=0.00113653;
    c7=-0.00020033;
    a=c1+xx*(c2+xx*(c3+xx*(c4+xx*(c5+xx*(c6+xx*c7)))));
    
    c1=-2.35619449;
    c2=0.12499612;
    c3=0.00005650;
    c4=-0.00637879;
    c5=0.00074348;
    c6=0.00079824;
    c7=-0.00029166;
    t=x+c1+xx*(c2+xx*(c3+xx*(c4+xx*(c5+xx*(c6+xx*c7)))));
    f=a*math.cos(t)/math.sqrt(x);
end

    if (u<0) then return -f 
    else return f end
end 

function DeLeeuw.SetModelParameters(self, nporo, nu_poisson, dCmedium, dCfluid, dCsolid, kappa, dim)

  -- lower case: non dimensional
  local nu      = nu_poisson or 0.25
  local nporo   = nporo or 0.0
  
  local CMedium =  dCmedium or 1.0  
  local CFluid  =  dCfluid or 0.0         -- 1/Pa (default: 0.0, i.e., incompressible)
  local CSolid  =  dCsolid or 0.0         -- default: incompressible => alpha =1.0
  

  local cscm    = CSolid/CMedium;
  local cfcm    = CFluid/CMedium;
  local BSkempton = (CMedium-CSolid)/((CMedium-CSolid)+ nporo*(CFluid-CSolid))
   
  local alpha = 1.0-cscm;                                 -- non-dimensional Biot coeff
  -- alpha = 0 -- for DEBUGGING
  local KS  = nporo * cfcm + (alpha-nporo)*cscm;          -- K*S
  
  self.modelParameter.cscm = cscm;
  self.modelParameter.cfcm = cfcm;
  self.modelParameter.nporo = nporo; 
  
  self.modelParameter.KS  = KS;
  self.modelParameter.alpha = alpha;

  -- non-dimensional parameter (required for roots)
  self.modelParameter.MM= (1-nu)*(1.0+3.0*KS/(2.0*(1.0+nu)*alpha*alpha))/(2.0*(1.0-2.0*nu));
  
  -- storativity
  self.modelParameter.S     = nporo * CFluid + (alpha-nporo) * CSolid;    -- Storativity S
  
  
  -- diffuson
  self.modelParameter.kappa = kappa
  
  -- elasticity parameters
  self.modelParameter.nu0   = nu;                     -- Poisson's ratio
  self.modelParameter.K  = 1.0/CMedium;               -- compression (bulk) modulus of the medium [Pa]
  self.modelParameter.Ku = 1.0/CMedium + (alpha*alpha)/self.modelParameter.S;
  
  self.modelParameter.EYoung = 3.0*self.modelParameter.K*(1.0-2.0*nu)   -- Young's modulus

  -- Express (lambda, mu=G) as (K,nu)
  local lambda = (3.0*self.modelParameter.K*nu)/(1.0+nu)                -- Lame's 1st parameter
  local G = (3.0*self.modelParameter.K*(1.0-2.0*nu))/(2.0+2.0*nu)       -- Lame's 2nd parameter (shear modulus)
  self.modelParameter.lambda = lambda
  self.modelParameter.G = G

  local Kv = 2.0*G*(1.0-nu)/(1.0-2.0*nu)                                -- uni-axial drained bulk modulus 
  self.modelParameter.Kv = Kv
  
 -- self.modelParameter.eps0 = self.qForce/(alpha + KS/alpha)
  self.modelParameter.p0 = self.qForce/(alpha + KS/alpha)  -- initial pore pressure
  
  local Kdim = {}
   Kdim[1] = Kv
   Kdim[2] = Kv/(2.0-2.0*nu)
   Kdim[3] = 1.0/CMedium
  
  print("Kd="..Kdim[1]..", "..Kdim[2]..", "..Kdim[3])
 -- alpha = 0 -- for DEBUGGING
  self.elemDiscParams = {}
  self.elemDiscParams[1] =
 -- { VOLUME = "INNER",  KAPPA = D, LAMBDA=lambda, MU = G, ALPHA=alpha, PHI=self.modelParameter.S, THETA=(alpha*alpha)*CMedium }

  { VOLUME = "INNER",  KAPPA = kappa, LAMBDA=lambda, MU = G, ALPHA=alpha, PHI=self.modelParameter.S, THETA=(alpha*alpha)/Kdim[dim]}-- THETA=(alpha*alpha)/Kdim[dim]} -- THETA=}

  
  print ("MassScale=".. self.elemDiscParams[1].THETA)
  
end

function DeLeeuw.PrintModelParameters(self)
  print ("Phi="..self.modelParameter.S)
  print ("K="..self.modelParameter.K)
  print ("Ku="..self.modelParameter.Ku)
  print ("Kv="..self.modelParameter.Kv)

  print ("E     ="..self.modelParameter.EYoung)
  print ("nu     ="..self.modelParameter.nu0) 
  
  print ("lambda="..self.modelParameter.lambda)
  print ("mu    ="..self.modelParameter.G)
  
  print ("Alpha =".. self.modelParameter.alpha)
  print ("Kappa ="..self.modelParameter.kappa)
  print ("p0    ="..self.modelParameter.p0)
  
  print ("MM="..self.modelParameter.MM)
  
  print ("derived:")
  local nu =self.modelParameter.nu0
  local Kdr3 = self.modelParameter.EYoung/(3.0*(1.0-2.0*nu))
  print ("K_dr(1)="..3.0*Kdr3*(1.0-nu)/(1.0+nu))
  print ("K_dr(2)="..self.modelParameter.Kv/(2.0-2.0*nu))
  print ("K_dr(3)="..Kdr3)

end

--From Verruijt: Poroelasticity (lecture notes)
function DeLeeuw.InitRoots(self, N)
  
  self.ROOT = {}
  
  local PI=math.pi;
  local MM = self.modelParameter.MM;

  local x=0;
  
  local i;
  -- local N = self.n_approx;
  for i=1,N do
  
    local a1=x+0.000001;
    local a2=x+1.2*PI;
    -- local f1=2.0*MM*a1*DeLeeuw.J0(a1)-DeLeeuw.J1(a1);
   --  local f2=2.0*MM*a2*DeLeeuw.J0(a2)-DeLeeuw.J1(a2);
    local f1=2.0*MM*a1*BesselJ0(a1)-BesselJ1(a1);
    local f2=2.0*MM*a2*BesselJ0(a2)-BesselJ1(a2);
    
    -- find root by bisection
    while (true) do
        local am=(a1+a2)/2.0;
        -- local fm=2.0*MM*am*DeLeeuw.J0(am)-DeLeeuw.J1(am);
        local fm=2.0*MM*am*BesselJ0(am)-BesselJ1(am);
        
        
        if (f1*fm<0.0) then 
          a2=am;
          f2=fm;
        elseif (f2*fm<0.0) then
          a1=am;
          f1=fm; 
        end
        
        if ((a2-a1)<0.0000001) then break end
    end

    -- proceed    
    self.ROOT[i]=(a1+a2)/2.0;
    x=a2;
    
    
     -- check:

    --local closeToZero = (1.0-eta*self.ROOT[j]*self.ROOT[j]/2.0)* math.tan(self.ROOT[j]) - self.ROOT[j]
    --if (math.abs(closeToZero)>1e-8) then
    --  print ("WARNING: rootCheck[i]="..closeToZero.." (i="..i..")")
    --end
    
    
  end
end


-- 
function DeLeeuw.GetCharTime(self)
  local  Kplus43G = (self.modelParameter.K+4.0/3.0*self.modelParameter.G)
  local beta = (self.modelParameter.alpha*self.modelParameter.alpha + self.modelParameter.S*Kplus43G)
  return (self.a*self.a*beta)/(Kplus43G * self.elemDiscParams[1].KAPPA);
end

function DeLeeuw.ComputePressure(self, r, t)

  local N = self.n_approx;
  local mc = self.modelParameter.MM;  
  
  local pp=0.0

  if (t==0) then return 1.0 end

  local jt = 0
  local i = 1
  while ((i<=N) --[[ and (jt <= 30.0)..]]) do
    local aa    = self.ROOT[i];
    local jj    = BesselJ0(aa) -- DeLeeuw.J0(aa);
    local jjx   = BesselJ0(aa*r/self.a); -- DeLeeuw.J0(aa*r/self.a);
    
    local coeff = (1.0-jjx/jj)/(1.0-mc*aa*aa-1/(4*mc));  
    jt    = aa*aa*(t/self.charTime);
    
    pp = pp + coeff*math.exp(-jt);
  
  if (false) then -- DEBUG
    print ("aa[i]= ".. aa) 
    print ("time[i]= ".. jt) 
    print ("coeff[i]= ".. coeff) 
    print ("decay[i]= ".. math.exp(-jt))
    print ("delta[i]= ".. coeff*math.exp(-jt))
    print ("pp[i]= ".. pp)
    print("==")
  end
   i = i+1
  end

  return pp;
end

local Cryer ={}
function Cryer.SetModelParameters(self, nporo, nu_poisson, dCmedium, dCfluid, dCsolid, kappa, dim)

  -- lower case: non dimensional
  local nu      = nu_poisson or 0.25
  local nporo   = nporo or 0.0
  
  local CMedium =  dCmedium or 1.0  
  local CFluid  =  dCfluid or 0.0         -- 1/Pa (default: 0.0, i.e., incompressible)
  local CSolid  =  dCsolid or 0.0         -- default: incompressible => alpha =1.0

  local cscm    = CSolid/CMedium;
  local cfcm    = CFluid/CMedium;
   
  local alpha = 1.0-cscm;                                               -- non-dimensional Biot coeff

  self.modelParameter.alpha = alpha;
  self.modelParameter.nporo = nporo; 
  
  self.modelParameter.S     = nporo * CFluid + (alpha-nporo) * CSolid;    -- Storativity S
 
  
  -- elasticity parameters
  self.modelParameter.nu0   = nu;                     -- Poisson's ratio
  self.modelParameter.K  = 1.0/CMedium;            -- compression (bulk) modulus of the medium [Pa]
  self.modelParameter.Ku = 1.0/CMedium + (alpha*alpha)/self.modelParameter.S;
  
  -- Express (lambda, mu=G) as (K,nu)
  local lambda = (3.0*self.modelParameter.K*nu)/(1.0+nu)
  local G = (3.0*self.modelParameter.K*(1.0-2.0*nu))/(2.0+2.0*nu)
  self.modelParameter.lambda = lambda
  self.modelParameter.G = G
   
  self.modelParameter.cscm = cscm;
  self.modelParameter.cfcm = cfcm;

  local KS_0  = nporo * cfcm + (alpha-nporo)*cscm;        -- K*S
  self.modelParameter.KS_0  = KS_0;
  
 -- self.modelParameter.eps0 = self.qForce/(alpha + KS_0/alpha)
  self.modelParameter.p0   = self.qForce/(alpha + KS_0/alpha)  -- initial pore pressure
  
  -- non-dimensional parameter (required for roots)
  local eta = (self.modelParameter.K+4.0/3.0*G)/(2.0*G)*(1.0+KS_0/(alpha*alpha))
  self.modelParameter.eta = eta; 
  
  
  local Kv = 2.0*G*(1.0-nu)/(1.0-2.0*nu)                                -- uni-axial drained bulk modulus 
  self.modelParameter.Kv = Kv
  
  local Kdim = {}
  Kdim[1] = Kv
  Kdim[2] = Kv/(2.0-2.0*nu)
  Kdim[3] = 1.0/CMedium

 print("Kd="..Kdim[1]..", "..Kdim[2]..", "..Kdim[3])

  self.elemDiscParams = {}
  self.elemDiscParams[1] =
  { VOLUME = "INNER",  KAPPA = kappa, LAMBDA=lambda, MU = G, ALPHA=alpha, PHI=self.modelParameter.S, THETA=(alpha*alpha)/Kdim[dim]}
  
  print ("Phi="..self.modelParameter.S)
  print ("K="..self.modelParameter.K)
  print ("Ku="..self.modelParameter.Ku)
  print ("lambda="..lambda)
  print ("mu="..self.modelParameter.G)
  print ("Alpha="..alpha)
  print ("Kappa="..kappa)
  print ("MassScale=".. self.elemDiscParams[1].THETA)
  print ("P0="..self.modelParameter.p0)
  
  
  
end



--From Verruijt: Poroelasticity (lecture notes)
function Cryer.InitRoots(self)
  self.ROOT = {}
  local N = self.n_approx;
  local eps= 1e-10 --0.0000001;

  local PI=math.pi;
  local eta = self.modelParameter.eta;

  local i,j,k;
  local a2,b2;

  -- print ("N="..N)
  for j=1,N do

    local a1 = eps+(2*j-1)*PI/2;   -- initial guess
    local da = PI
    local b1 = math.tan(a1)*(1.0-eta*a1*a1/2.0)-a1;
    
    local a2;
    for i=1,5 do 
      for k=1,N do
      
        a2 = a1+da/N;  -- next guess 
      
        local b2 = math.tan(a2)*(1.0-eta*a2*a2/2.0)-a2;
        if (b2*b1<0) then 
            k=N; da=a2-a1;
        else 
            a1=a2; b1=b2;  -- 
        end
      end
    end
    
    
    self.ROOT[j]=(a1+a2)/2.0;
    
    
    -- check:
    
    local closeToZero = (1.0-eta*self.ROOT[j]*self.ROOT[j]/2.0)* math.tan(self.ROOT[j]) - self.ROOT[j]
    if (math.abs(closeToZero)>1e-12) then
      print ("WARNING: rootCheck[j]="..closeToZero.." (j="..j..")")
    end
  end
end




function cryer2d:create_domain(numRefs, numPreRefs)

  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  
  --local refProjector = DomainRefinementProjectionHandler(dom)
  -- refProjector:set_callback("INNER", projector)
  --SaveDomain(dom, "Sphere3d.ugx");

  --SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "CryerSpheres.ugx", 1)
  return dom
end




--[[
  init all variables
  
  kperm: permeability
  nporo:
  
  cmedium, cfluid, csolid: compressibilities (1/[Pa])
  gamma_f: volumetric weight
--]]
function cryer2d:init(kperm, nporo, nu, cmedium, cfluid, csolid, gamma_f)
  gamma_f = gamma_f or 1.0 -- kg/m^3
  Cryer.SetModelParameters(self, nporo, nu, cmedium, cfluid, csolid, kperm*gamma_f, self.dim)
  Cryer.InitRoots(self)
  self.charTime = self:get_char_time(); 
  print ("charTime="..self.charTime)
end


function Cryer.ComputePressure(self, t)

local N = self.n_approx;
local eta = self.modelParameter.eta;
local pp=0.0

if (t==0) then return 1.0 end
 local jt = 0
  local i = 1
  while ((i<=N) and (jt <= 20.0)) do
  local aa  = self.ROOT[i];
  local sinaa = math.sin(aa)
  local cosaa = math.cos(aa)
  local coeff   = (sinaa-aa) / (eta*aa*cosaa/2.0 + (eta-1.0)*sinaa);
  jt  = aa*aa*t/self.charTime;
  pp = pp + eta*coeff*math.exp(-jt);
  

  print ("aa[i]= ".. aa) 
  print ("time[i]= ".. jt) 
  print ("coeff[i]= ".. coeff) 
  print ("decay[i]= ".. math.exp(-jt))
  print ("delta[i]= ".. eta*coeff*math.exp(-jt))
  print ("pp[i]= ".. pp)
  print("=====")
  i =i+1
  --if (jt>20) then i=N+1 end
end
return(pp);
end


function Cryer.GetCharTime(self)
  local  Kplus43G = (self.modelParameter.K+4.0/3.0*self.modelParameter.G)
  local beta = (self.modelParameter.alpha*self.modelParameter.alpha + self.modelParameter.S*Kplus43G)
  return (self.a*self.a*beta)/(Kplus43G * self.elemDiscParams[1].KAPPA);
end

-- pressure (at center node)
function cryer2d:compute_pressure(x,y,t)
  return Cryer.ComputePressure(self, t)
end


-- characteristic time
function cryer2d:get_char_time()
 return Cryer.GetCharTime(self)
end


-- boundary conditions
function cryer2d:add_boundary_conditions(domainDisc, bStationary)
  
  
 local doStationary = bStationary or false

 local neumannX = NeumannBoundaryFE("ux")
 neumannX:add("Sphere2dNormalX", "RIM", "INNER")
  
 local neumannY = NeumannBoundaryFE("uy")
 neumannY:add("Sphere2dNormalY", "RIM", "INNER")
 
 if doStationary then 
    neumannX:set_stationary()
    neumannY:set_stationary() 
  end
 
 domainDisc:add(neumannX)
 domainDisc:add(neumannY)

  local dirichlet = DirichletBoundary()
  dirichlet:add(0.0, "p", "RIM,EAST,WEST,NORTH,SOUTH")
  
  dirichlet:add(0.0, "ux", "CENTER, NORTH, SOUTH")
  dirichlet:add(0.0, "uy", "CENTER, EAST,WEST")
  
  --dirichlet:add(0.0, "ux", "CENTER")
 -- dirichlet:add(0.0, "uy", "CENTER")
  
  --dirichlet:add(0.0, "ux", "NORTH, SOUTH, RIM")
  --dirichlet:add(0.0, "uy", "EAST, WEST, RIM")

  domainDisc:add(dirichlet)
 
end




function cryer2d:add_elem_discs(domainDisc, bStationary)
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
end

-- initial values
function cryer2d:interpolate_start_values(u, startTime)
  Interpolate(self.modelParameter.p0, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
end


-- post processing (after each step)
function cryer2d:post_processing(u, step, time)
  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..(Integral(u, "p", "CENTER")/self.modelParameter.p0).."\t"..Cryer.ComputePressure(self, time))
end

--[[ Cf. Verruijt:

  Cylindrical soil sample (radius a) 
  - is placed between two stiff horizontal planes and,
  - subjected to radial pressure (of magnitude q)
  
  
-- ]]
-- OLD:
--[[
deleeuw2d = {

  gridName = "../grids/cryer2d-with-projectorb.ugx",
  dim = 2,
  cpu = 1,
  
  porder = 1,
  uorder = 2, -- should be 2
  mandatorySubsets ={"INNER"},
  
  a = 0.5, 
  qForce = 1.0,
  
  n_approx = 100,
  
  modelParameter = {},
  elemDiscParams = {},

  bRAP = false,
}
--]]

deleeuw2d = {

  gridName = "../grids/cryer2d-with-projectorb-tri.ugx",
  dim = 2,
  cpu = 1,
  
  porder = 1,
  uorder = 1, -- should be 2
  vStab = 0.0, --1.0/12.0,  --*h^2 /(lambda+2G)
  mandatorySubsets ={"INNER"},
  
  a = 0.5, 
  qForce = 1.0,
  
  n_approx = 4096,
  
  modelParameter = {},
  elemDiscParams = {},
  
  bRAP = false,
  bAdjustTransfers = true

}


-- Read parameters from command line.
function deleeuw2d:parse_cmd_args()

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


function deleeuw2d:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end


function deleeuw2d:init(kperm, nporo, nu, cmedium, cfluid, csolid, gamma_f)
  gamma_f = gamma_f or 1.0 -- kg/m^3
  DeLeeuw.SetModelParameters(self, nporo, nu, cmedium, cfluid, csolid, kperm*gamma_f, self.dim)
  DeLeeuw.InitRoots(self, self.n_approx)
  self.charTime = self:get_char_time(); 
  DeLeeuw.PrintModelParameters(self)
  
  -- Reset stabilization
  if (self.vStab and self.vStab ~= 0.0) then
      print ("vStab(preReset)="..self.vStab)
      local lambda = self.modelParameter.lambda
      local G = self.modelParameter.G 
      local factor = lambda+2.0*G
      self.vStab = self.vStab/factor
      print ("vStab(postReset)="..self.vStab)
  end
  
  print ("charTime="..self.charTime)
end

-- characteristic time
function deleeuw2d:get_char_time()
 return DeLeeuw.GetCharTime(self)
end


function deleeuw2d:add_elem_discs(domainDisc, bStationary)
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
  
   if (self.vStab and self.vStab~= 0.0 and self.porder==self.uorder) then     -- Add stabilization?
    print ("vStab="..self.vStab)
    CommonAddBiotStabDiscs(self, domainDisc)   
  end
end

function deleeuw2d:add_uzawa_discs(domainDisc)
  CommonAddMassMatrixDiscs(self, domainDisc)
end


function deleeuw2d:add_boundary_conditions(domainDisc, bStationary)
  
  
 local doStationary = bStationary or false
 --[[
 local unity = {}
 unity[0] = ConstUserVector()
 unity[0]:set_entry(0, -1.0)
 
 unity[1] = ConstUserVector()
 unity[1]:set_entry(1, -1.0)
 
 
 unity[0]:print()
 unity[1]:print()
 print( )

 local neumannY = NeumannBoundaryFE("uy")
 local neumannX = NeumannBoundaryFE("ux")
 
 if doStationary then 
    neumannX:set_stationary()
    neumannY:set_stationary() 
  end
  
 for i=1,#self.elemDiscParams do
  --  self.flowDisc[i], self.dispDisc[i] = CreateElemDiscs(self.elemDiscParams[i], self.dim, bStationary)
    print(self.flowDisc[i])
    local gradP ={}
    gradP[0] = ScaleAddLinkerVectorVector()
    gradP[0]:add(self.flowDisc[i]:gradient(), unity[0])
    
    gradP[1] = ScaleAddLinkerVectorVector()
    gradP[1]:add(unity[1],self.flowDisc[i]:gradient())
    gradP[1]:add(unity[1], unity[1])
    
    -- neumannY:add("Sphere2dNormalY", "RIM", "INNER")
    neumannX:add(gradP[0], "RIM", "INNER")
    --neumannX:add(gradP[0], "RIM", "INNER")

    -- neumannX:add(1, "RIM", "INNER")
   -- neumannY:add(1, "RIM", "INNER")
   -- neumannX:add(self.flowDisc[i]:gradient(), "RIM", "INNER")
 
  end
 domainDisc:add(neumannX)
 domainDisc:add(neumannY)

--]]

local neumannX = NeumannBoundaryFE("ux")
 neumannX:add("Sphere2dNormalX", "RIM", "INNER")
  
 local neumannY = NeumannBoundaryFE("uy")
 neumannY:add("Sphere2dNormalY", "RIM", "INNER")
 
 if doStationary then 
    neumannX:set_stationary()
    neumannY:set_stationary() 
  end
 
 domainDisc:add(neumannX)
 domainDisc:add(neumannY)
 
  if (self.bAdjustTransfers) then print ("Dirichlet w/ adjusting transfers") end
  
  -- local dirichlet = DirichletBoundary(false, self.bAdjustTransfers)
  local dirichlet = DirichletBoundary(false)
  dirichlet:add(0.0, "p", "RIM,EAST,WEST,NORTH,SOUTH")
  
  dirichlet:add(0.0, "ux", "CENTER, NORTH, SOUTH")
  dirichlet:add(0.0, "uy", "CENTER, EAST,WEST")
  
  --dirichlet:add(0.0, "ux", "CENTER")
  --dirichlet:add(0.0, "uy", "CENTER")
  
  domainDisc:add(dirichlet)
 
end



-- initial values
function deleeuw2d:interpolate_start_values(u, startTime)
  Interpolate(self.modelParameter.p0, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
end


-- create error estimator
function deleeuw2d:error_estimator()
  local p = self.elemDiscParams[1]
  local gamma2 = (p.LAMBDA+2*p.MU)*(p.LAMBDA+2*p.MU)
  print("deleeuw3dTet:error_estimator: gamma="..gamma2)
  
  local biotErrorEst = CompositeGridFunctionEstimator()
  biotErrorEst:add(L2ComponentSpace("p", 2))       
  biotErrorEst:add(H1SemiComponentSpace("ux", 4, gamma2, p.VOLUME))  
  biotErrorEst:add(H1SemiComponentSpace("uy", 4, gamma2, p.VOLUME))
  -- biotErrorEst:add(H1SemiComponentSpace("uz", 4, gamma2, p.VOLUME))
 
  return biotErrorEst
end

__DELEEUW_2D_GLOBAL__  = nil
function DeLeeuw2DPressure(x, y, t)
  return DeLeeuw.ComputePressure(__DELEEUW_2D_GLOBAL__, math.sqrt(x*x+y*y), t)
end

-- post processing (after each step)
function deleeuw2d:post_processing(u, step, time)
 
  local simP = Integral(u, "p", "CENTER")/self.modelParameter.p0
  local refP = DeLeeuw.ComputePressure(self, 0.0, time)
  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..simP.."\t"..refP.."\t"..simP-refP)
 
 
  __DELEEUW_2D_GLOBAL__ = self

  
  local uref = u:clone()
  Interpolate("DeLeeuw2DPressure", uref, "p", "INNER,CENTER", time)
  
  local vtk = VTKOutput()
  vtk:print("DeLeeuw2D_Sol.vtu", u, step, time)
  vtk:print("DeLeeuw2D_Ref.vtu", uref, step, time)
  
  VecScaleAdd2(uref, 1.0, uref, -1.0, u)
  vtk:print("DeLeeuw2D_Err.vtu", uref, step, time)
  
  local normP=L2Norm(u, "p", 2)
  local errP=L2Norm(uref, "p", 2)
  print ("errP:\t"..time.."\t"..time/self.charTime.."\t"..errP.."\t"..normP)
end



deleeuw3d = {

  -- gridName = "../grids/cryer3d-with-projectorb.ugx",
  gridName = "../grids/cryer3d-noprojector.ugx",
  dim = 3,
  cpu = 1,
  
  porder = 1,
  uorder = 2, -- should be 2
  vStab = 0.0, 
  
  mandatorySubsets ={"INNER"},
  
  a = 0.5, 
  qForce = 1.0,
  
  n_approx = 100,
  
  modelParameter = {},
  elemDiscParams = {},
  
  bRAP = true,

}

function deleeuw3d:parse_cmd_args()
end

function deleeuw3d:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end


function deleeuw3d:init(kperm, nporo, nu, cmedium, cfluid, csolid, gamma_f)
  gamma_f = gamma_f or 1.0 -- kg/m^3
  DeLeeuw.SetModelParameters(self, nporo, nu, cmedium, cfluid, csolid, kperm*gamma_f, self.dim)
  DeLeeuw.InitRoots(self, self.n_approx)
  self.charTime = self:get_char_time(); 
  DeLeeuw.PrintModelParameters(self)
  print ("charTime="..self.charTime)
end

-- characteristic time
function deleeuw3d:get_char_time()
 return DeLeeuw.GetCharTime(self)
end


function deleeuw3d:add_elem_discs(domainDisc, bStationary)
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
 
  if (self.vStab and self.vStab~= 0.0 and self.porder==self.uorder) then     -- Add stabilization?
   print ("vStab="..self.vStab)
    CommonAddBiotStabDiscs(self, domainDisc)   
  end
end

function deleeuw3d:add_uzawa_discs(domainDisc)
  CommonAddMassMatrixDiscs(self, domainDisc)
end


function deleeuw3d:add_boundary_conditions(domainDisc, bStationary)
  
  
local doStationary = bStationary or false

local neumannX = NeumannBoundaryFE("ux")
 neumannX:add("Cylinder3dNormalX", "RIM", "INNER")
  
 local neumannY = NeumannBoundaryFE("uy")
 neumannY:add("Cylinder3dNormalY", "RIM", "INNER")
 
 if doStationary then 
    neumannX:set_stationary()
    neumannY:set_stationary() 
  end
 
 domainDisc:add(neumannX)
 domainDisc:add(neumannY)
 
 -- local dirichlet = DirichletBoundary(false, true)
  local dirichlet = DirichletBoundary(false)
  dirichlet:add(0.0, "p", "RIM,EAST,WEST,NORTH,SOUTH")

  dirichlet:add(0.0, "ux", "CENTER")
  dirichlet:add(0.0, "uy", "CENTER")
  dirichlet:add(0.0, "uz", "CENTER")
  
  dirichlet:add(0.0, "uz", "TOP")
  dirichlet:add(0.0, "uz", "BOT")
  
  domainDisc:add(dirichlet)
 
end



-- initial values
function deleeuw3d:interpolate_start_values(u, startTime)
  Interpolate(self.modelParameter.p0, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
end


-- post processing (after each step)
function deleeuw3d:post_processing(u, step, time)
  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..(Integral(u, "p", "CENTER")/self.modelParameter.p0).."\t"..DeLeeuw.ComputePressure(self, 0.0, time))
end




deleeuw3dTet = {

  -- gridName = "../grids/cryer3d-with-projectorb.ugx",
  gridName = "../grids/deleeuw3d-cylinder-tet.ugx",
  dim = 3,
  cpu = 1,
  
  porder = 1,
  uorder = 2, -- should be 2
  vStab = 0.0, 
  mandatorySubsets ={"INNER"},
  
  a = 0.5, 
  qForce = 1.0,
  
  n_approx = 100,
  
  modelParameter = {},
  elemDiscParams = {},
  bRAP = false,

}

function deleeuw3dTet:parse_cmd_args()
end

function deleeuw3dTet:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end


function deleeuw3dTet:init(kperm, nporo, nu, cmedium, cfluid, csolid, gamma_f)
  gamma_f = gamma_f or 1.0 -- kg/m^3
  DeLeeuw.SetModelParameters(self, nporo, nu, cmedium, cfluid, csolid, kperm*gamma_f, self.dim)
  DeLeeuw.InitRoots(self, self.n_approx)
  self.charTime = self:get_char_time(); 
  DeLeeuw.PrintModelParameters(self)
  print ("charTime="..self.charTime)
end

-- characteristic time
function deleeuw3dTet:get_char_time()
 return DeLeeuw.GetCharTime(self)
end


function deleeuw3dTet:add_elem_discs(domainDisc, bStationary)
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
end

function deleeuw3dTet:add_uzawa_discs(domainDisc)
  CommonAddMassMatrixDiscs(self, domainDisc)
end


function deleeuw3dTet:add_boundary_conditions(domainDisc, bStationary)
  
  
local doStationary = bStationary or false

local neumannX = NeumannBoundaryFE("ux")
 neumannX:add("Cylinder3dNormalX", "RIM", "INNER")
  
 local neumannY = NeumannBoundaryFE("uy")
 neumannY:add("Cylinder3dNormalY", "RIM", "INNER")
 
 if doStationary then 
    neumannX:set_stationary()
    neumannY:set_stationary() 
  end
 
 domainDisc:add(neumannX)
 domainDisc:add(neumannY)
 
 -- local dirichlet = DirichletBoundary(true, false) -- dirichlet columns, adjust interpolation 
   local dirichlet = DirichletBoundary(true) -- dirichlet columns, adjust interpolation 
  dirichlet:add(0.0, "p", "RIM")

  dirichlet:add(0.0, "ux", "CENTER")
  dirichlet:add(0.0, "uy", "CENTER")
  dirichlet:add(0.0, "uz", "CENTER")
  
  dirichlet:add(0.0, "uz", "TOP")
  dirichlet:add(0.0, "uz", "BOT")
  dirichlet:add(0.0, "uz", "RIM")
  
  domainDisc:add(dirichlet)
 
end



-- initial values
function deleeuw3dTet:interpolate_start_values(u, startTime)
  Interpolate(self.modelParameter.p0, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
end


-- create error estimator
function deleeuw3dTet:error_estimator()
  local p = self.elemDiscParams[1]
  local gamma2 = (p.LAMBDA+2*p.MU)*(p.LAMBDA+2*p.MU)
  print("deleeuw3dTet:error_estimator: gamma="..gamma2)
  
  local biotErrorEst = CompositeGridFunctionEstimator()
        
  biotErrorEst:add(H1SemiComponentSpace("ux", 4, gamma2, p.VOLUME))  
  biotErrorEst:add(H1SemiComponentSpace("uy", 4, gamma2, p.VOLUME))
  biotErrorEst:add(H1SemiComponentSpace("uz", 4, gamma2, p.VOLUME))
  biotErrorEst:add(L2ComponentSpace("p", 2)) 
  
  return biotErrorEst
end


-- post processing (after each step)
function deleeuw3dTet:post_processing(u, step, time)
  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..(Integral(u, "p", "CENTER")/self.modelParameter.p0).."\t"..DeLeeuw.ComputePressure(self, 0.0, time))
end


---------------------------------------------
--
-- 3D problem
--
---------------------------------------------
cryer3d = {
 
  gridName = "../grids/eighth-sphere48-2layer.ugx",
  dim = 3,
  cpu = 1,
  
  porder = 1,
  uorder = 2,
  mandatorySubsets ={"INNER", "SPHERE1"},
  
  a = 1.0, 
  qForce = 1.0,
  
  n_approx = 110,
  
  modelParameter = {},
  elemDiscParams = {}

}


function cryer3d:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end


function cryer3d:add_elem_discs(domainDisc, bStationary)
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
end


function cryer3d:add_uzawa_discs(domainDisc)
  CommonAddMassMatrixDiscs(self, domainDisc)
end

--- Boundary conditions.
function cryer3d:add_boundary_conditions(domainDisc, bStationary)
  
  
 local doStationary = bStationary or false

 local neumannX = NeumannBoundaryFE("ux")
 local neumannY = NeumannBoundaryFE("uy")
 local neumannZ = NeumannBoundaryFE("uz")
 
 neumannX:add("Sphere3dNormalX", "SPHERE1", "INNER")
 neumannY:add("Sphere3dNormalY", "SPHERE1", "INNER")
 neumannZ:add("Sphere3dNormalZ", "SPHERE1", "INNER")
 
 if doStationary then 
    neumannX:set_stationary()
    neumannY:set_stationary() 
    neumannZ:set_stationary() 
 end
 
 domainDisc:add(neumannX)
 domainDisc:add(neumannY)
 domainDisc:add(neumannZ)

  local dirichlet = DirichletBoundary()
  dirichlet:add(0.0, "p", "SPHERE1, ARC_XY, ARC_YZ, ARC_XZ, PX, PY, PZ")
  dirichlet:add(0.0, "ux", "ARC_YZ, SURF_YZ, AXIS_Y, AXIS_Z, PY, PZ, CENTER")
  dirichlet:add(0.0, "uy", "ARC_XZ, SURF_XZ, AXIS_X, AXIS_Z, PX, PZ, CENTER")
  dirichlet:add(0.0, "uz", "ARC_XY, SURF_XY, AXIS_X, AXIS_Y, PX, PY, CENTER")

  domainDisc:add(dirichlet)
 
end

-- init all variables
function cryer3d:init(kperm, nporo, nu, cmedium, cfluid, csolid, volumetricweight)
  volumetricweight = volumetricweight or 1.0 -- kg/m^3
  Cryer.SetModelParameters(self, nporo, nu, cmedium, cfluid, csolid, kperm*volumetricweight, self.dim)
  Cryer.InitRoots(self)
  self.charTime = self:get_char_time(); 
  print ("charTime="..self.charTime)
end

function cryer3d:get_char_time()
 return Cryer.GetCharTime(self)
end

-- Initial values.
function cryer3d:interpolate_start_values(u, startTime)
  Interpolate(self.modelParameter.p0, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
  Interpolate(0.0, u, "uz", startTime)
end

-- Post processing (after each step)
function cryer3d:post_processing(u, step, time)
  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..(Integral(u, "p", "CENTER")/self.modelParameter.p0).."\t"..Cryer.ComputePressure(self, time))
end

cryer3dTet = {
 
  gridName = "../grids/cryer3d-sphere-tet.ugx",
  dim = 3,
  cpu = 1,
  
  porder = 1,
  uorder = 2,
  mandatorySubsets ={"INNER", "SPHERE"},
  
  a = 1.0, 
  qForce = 1.0,
  
  n_approx = 110,
  
  modelParameter = {},
  elemDiscParams = {}

}


function cryer3dTet:create_domain(numRefs, numPreRefs)
  local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
  return dom
end


function cryer3dTet:add_elem_discs(domainDisc, bStationary)
  CommonAddBiotElemDiscs(self, domainDisc, bStationary)
end


function cryer3dTet:add_uzawa_discs(domainDisc)
  CommonAddMassMatrixDiscs(self, domainDisc)
end

--- Boundary conditions.
function cryer3dTet:add_boundary_conditions(domainDisc, bStationary)
  
  
 local doStationary = bStationary or false

 local neumannX = NeumannBoundaryFE("ux")
 local neumannY = NeumannBoundaryFE("uy")
 local neumannZ = NeumannBoundaryFE("uz")
 
 neumannX:add("Sphere3dNormalX", "SPHERE", "INNER")
 neumannY:add("Sphere3dNormalY", "SPHERE", "INNER")
 neumannZ:add("Sphere3dNormalZ", "SPHERE", "INNER")
 
 if doStationary then 
    neumannX:set_stationary()
    neumannY:set_stationary() 
    neumannZ:set_stationary() 
 end
 
 domainDisc:add(neumannX)
 domainDisc:add(neumannY)
 domainDisc:add(neumannZ)

  local dirichlet = DirichletBoundary()
  dirichlet:add(0.0, "p", "SPHERE")
  dirichlet:add(0.0, "ux", "CENTER")
  dirichlet:add(0.0, "uy", "CENTER")
  dirichlet:add(0.0, "uz", "CENTER")

  domainDisc:add(dirichlet)
 
end

-- init all variables
function cryer3dTet:init(kperm, nporo, nu, cmedium, cfluid, csolid, volumetricweight)
  volumetricweight = volumetricweight or 1.0 -- kg/m^3
  Cryer.SetModelParameters(self, nporo, nu, cmedium, cfluid, csolid, kperm*volumetricweight, self.dim)
  Cryer.InitRoots(self)
  self.charTime = self:get_char_time(); 
  print ("charTime="..self.charTime)
end

function cryer3dTet:get_char_time()
 return Cryer.GetCharTime(self)
end

-- Initial values.
function cryer3dTet:interpolate_start_values(u, startTime)
  Interpolate(self.modelParameter.p0, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
  Interpolate(0.0, u, "uz", startTime)
end

-- Post processing (after each step)
function cryer3dTet:post_processing(u, step, time)
  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..(Integral(u, "p", "CENTER")/self.modelParameter.p0).."\t"..Cryer.ComputePressure(self, time))
end
