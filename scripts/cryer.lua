
function Cryer2dBndX(x,y,t)
  return 1.0*x/math.sqrt(x*x+y*y)
end

function Cryer2dBndY(x,y,t)
  return 1.0*y/math.sqrt(x*x+y*y)
end

cryer2d = {

 
  gridName = "../grids/cryer2d.ugx",
  dim = 2,
  cpu = 1,
  
  porder = 1,
  uorder = 2,
  mandatorySubsets ={"INNER"},
  
  a = 0.5, 
  qForce = 1.0,
  
  n_approx = 200,
  
  modelParameter = {},
  elemDiscParams = {}

}


function cryer2d:setModelParameters_(nu_poisson, nporo, cmedium, cfluid, csolid, D)

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
  local lambda = (3.0*self.modelParameter.K*nu)/(1.0+nu)
  local mu = (3.0*self.modelParameter.K*(1.0-2.0*nu))/(2.0+2.0*nu)
  self.modelParameter.mu = mu
   
  self.modelParameter.cscm = cscm;
  self.modelParameter.cfcm = cfcm;

  local KS_0  = poro * cfcm + (alpha-poro)*cscm;                        -- K*S
  self.modelParameter.KS_0  = KS_0;
  
  -- p0
  self.modelParameter.p0 = self.qForce/(alpha + KS_0/alpha)
  
  
  -- non-dimensional parameter (required for roots)
  local eta = (self.modelParameter.K+4.0/3.0*mu)/(2.0*mu)*(1.0+KS_0/(alpha*alpha))
  self.modelParameter.eta = eta; 
  
  
  
  self.elemDiscParams = {}
  self.elemDiscParams[1] =
  {VOLUME = "INNER",  KAPPA = D, LAMBDA=lambda, MU = mu, ALPHA=alpha, PHI=self.modelParameter.S}
  
  print ("K="..self.modelParameter.K)
  print ("Phi="..self.modelParameter.S)
  print ("Alpha="..alpha)
  print ("Kappa="..D)
  
  print ("charTime="..self:get_char_time())
  
end

--Verruijt: Poroelasticity (lecture notes)
function cryer2d:initRoots_()

  self.ROOT = {}
  local N = self.n_approx;
  local eps=0.000001;

  local PI=math.pi;
  local eta = self.modelParameter.eta;

  local i,j,k;
  local a2,b2;

  for j=1,N do

    local a1 = eps+(2*j-1)*PI/2;
    local a2;
    local da = PI
    local b1 = math.tan(a1)*(1.0-eta*a1*a1/2.0)-a1;
    
    
    for i=1,4 do
      for k=0,(N-1) do
      
        a2=a1+da/N;
      
        local b2 = math.tan(a2)*(1.0-eta*a2*a2/2.0)-a2;
        if (b2*b1<0) then k=N; da=a2-a1;
        else a1=a2; b1=b2;  end
      end
    end
    self.ROOT[j]=(a1+a2)/2.0;
  end

end -- function 'initRoots_'



-- init all variables
function cryer2d:init(kperm, poro, nu, cmedium, cfluid, csolid, volumetricweight)

  volumetricweight = volumetricweight or 1.0 -- kg/m^3
  cryer2d:setModelParameters_(nu, poro, cmedium, cfluid, csolid, kperm*volumetricweight)
  cryer2d:initRoots_()
  self.charTime = self:get_char_time(); 
end


-- pressure (at center node)
function cryer2d:compute_pressure(x,y,t)

local N = self.n_approx;
local eta = self.modelParameter.eta;
local pp=0.0

for i=1,N do
  local aa  = self.ROOT[i];
  local sinaa = math.sin(aa)
  local cosaa = math.cos(aa)
  local val   = (sinaa-aa) / (eta*aa*cosaa/2.0 + (eta-1.0)*sinaa);
  local jt  = aa*aa*t/self.charTime;
  pp = pp + eta*val*math.exp(-jt);
  
  if (jt>20) then i=N+1 end
end
return(self.modelParameter.p0*pp);
end

-- characteristic time
function cryer2d:get_char_time()
  local  Kplus43G = (self.modelParameter.K+4.0/3.0*self.modelParameter.mu)
  local beta = (self.modelParameter.alpha*self.modelParameter.alpha + self.modelParameter.S*Kplus43G)
  return (self.a*self.a*beta)/(Kplus43G * self.elemDiscParams[1].KAPPA);
end


-- boundary conditions
function cryer2d:add_boundary_conditions(domainDisc)
  local dirichlet = DirichletBoundary()
  dirichlet:add(0.0, "p", "RIM")
  domainDisc:add(dirichlet)
  
  
 local neumannX = NeumannBoundaryFE("ux")
 neumannX:add("Cryer2dBndX", "RIM", "INNER")
 domainDisc:add(neumannX)
  
 local neumannY = NeumannBoundaryFE("uy")
 neumannY:add("Cryer2dBndY", "RIM", "INNER")
 domainDisc:add(neumannY)
end



function cryer2d:add_elem_discs(domainDisc)
  
  self.flowDisc = {}
  self.dispDisc = {}
  
  for i=1,#self.elemDiscParams do
    self.flowDisc[i], self.dispDisc[i] = CreateElemDiscs(self.elemDiscParams[i], self.dim)
    domainDisc:add(self.flowDisc[i])
    domainDisc:add(self.dispDisc[i])
  end
end

-- initial values
function cryer2d:interpolate_start_values(u, startTime)
  Interpolate(self.modelParameter.p0, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
end


-- post processing (after each step)
function cryer2d:post_processing(u, step, time)
  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..Integral(u, "p", "CENTER").."\t"..self:compute_pressure(0.0,0.0,time))
end


