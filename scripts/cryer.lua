
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


-- Neumann boundary (2D)
function Cryer2dBndX(x,y,t)
  return 1.0*x/math.sqrt(x*x+y*y)
end

function Cryer2dBndY(x,y,t)
  return 1.0*y/math.sqrt(x*x+y*y)
end

-- Neumann boundary (3D)
function Cryer3dBndX(x,y,z,t)
  return 1.0*x/math.sqrt(x*x+y*y+z*z)
end

function Cryer3dBndY(x,y,z,t)
  return 1.0*y/math.sqrt(x*x+y*y+z*z)
end

function Cryer3dBndY(x,y,z,t)
  return 1.0*z/math.sqrt(x*x+y*y+z*z)
end

local function CryerSetModelParameters_(self, nu_poisson, nporo, cmedium, cfluid, csolid, D)

  -- lower case: non dimensional
  
  
  local nu = nu_poisson or 0.2
  local poro = nporo or 0.0
  
  local CMedium =  cmedium or 1.0  
  local CFluid  =  cfluid or 0.0           -- default: incompressible
  local CSolid  =  csolid or 0.0           -- default: incompressible => alpha =1.0

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
  print ("lambda="..lambda)
  print ("mu="..mu)
  print ("Alpha="..alpha)
  print ("Kappa="..D)
  print ("P0="..self.modelParameter.p0)
  
  
  
end


--Verruijt: Poroelasticity (lecture notes)
local function CryerInitRoots_(self)
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
end




local function createSphereProjector(dom)
     -- args: dom, center, radius, eps
    ProjectVerticesToSphere(dom, {0.0, 0.0}, 0.25, 0.05)  --
    ProjectVerticesToSphere(dom, {0.0, 0.0}, 0.5, 0.05)  --
    return SphericalFalloffProjector(dom, {0.0, 0.0}, 0.5, 0.05) -- For: grid, pos, center, inner_radius, outer_radius --
end

function cryer2d:create_domain(numRefs, numPreRefs)

local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)
--ProjectVerticesToSphere (dom, {0.0, 0.0, 0.0}, 0.5, 0.01)
--local projector = createSphereProjector(dom)

--local refProjector = DomainRefinementProjectionHandler(dom)
-- refProjector:set_callback("INNER", projector)
--SaveDomain(dom, "Sphere3d.ugx");

--SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "CryerSpheres.ugx", 1)
return dom
end




-- init all variables
function cryer2d:init(kperm, poro, nu, cmedium, cfluid, csolid, volumetricweight)
  volumetricweight = volumetricweight or 1.0 -- kg/m^3
  CryerSetModelParameters_(self, nu, poro, cmedium, cfluid, csolid, kperm*volumetricweight)
  CryerInitRoots_(self)
  self.charTime = self:get_char_time(); 
  print ("charTime="..self.charTime)
end


local function CryerComputePressure(self, x,y,t)

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


local function CryerCharTime(self)

 local  Kplus43G = (self.modelParameter.K+4.0/3.0*self.modelParameter.mu)
  local beta = (self.modelParameter.alpha*self.modelParameter.alpha + self.modelParameter.S*Kplus43G)
  return (self.a*self.a*beta)/(Kplus43G * self.elemDiscParams[1].KAPPA);
end

-- pressure (at center node)
function cryer2d:compute_pressure(x,y,t)
  return CryerComputePressure(self, x,y,t)
end


-- characteristic time
function cryer2d:get_char_time()
 return CryerCharTime(self)
end


-- boundary conditions
function cryer2d:add_boundary_conditions(domainDisc, bStationary)
  
  
 local doStationary = bStationary or false

 local neumannX = NeumannBoundaryFE("ux")
 neumannX:add("Cryer2dBndX", "RIM", "INNER")
  
 local neumannY = NeumannBoundaryFE("uy")
 neumannY:add("Cryer2dBndY", "RIM", "INNER")
 
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
  
  --dirichlet:add(0.0, "ux", "NORTH, SOUTH, RIM")
  --dirichlet:add(0.0, "uy", "EAST, WEST, RIM")

  domainDisc:add(dirichlet)
 
end

local function CryerAddElemDiscs(self, domainDisc, bStationary)
  
  self.flowDisc = {}
  self.dispDisc = {}
  
  for i=1,#self.elemDiscParams do
    self.flowDisc[i], self.dispDisc[i] = CreateElemDiscs(self.elemDiscParams[i], self.dim, bStationary)
    domainDisc:add(self.flowDisc[i])
    domainDisc:add(self.dispDisc[i])
  end
end

function cryer2d:add_elem_discs(domainDisc, bStationary)
  CryerAddElemDiscs(self, domainDisc, bStationary)
end

-- initial values
function cryer2d:interpolate_start_values(u, startTime)
  Interpolate(self.modelParameter.p0, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
end


-- post processing (after each step)
function cryer2d:post_processing(u, step, time)
  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..(Integral(u, "p", "CENTER")/self.modelParameter.p0).."\t"..CryerComputePressure(self, 0.0,0.0,time))
end


---------------------------------------------
--
-- 3D problem
--
---------------------------------------------
cryer3d = {
 
  gridName = "../grids/cryer3d-flat.ugx",
  dim = 3,
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


function cryer3d:create_domain(numRefs, numPreRefs)

-- local dom = util.CreateAndDistributeDomain(self.gridName, numRefs, numPreRefs, self.mandatorySubsets)


  -- create Instance of a Domain
  local dom = Domain()

  local gridName = self.gridName
  local mandatorySubsets = self.mandatorySubsets
  
  -- load domain
  write("Loading Domain "..gridName.." ... ") 
  LoadDomain(dom, gridName)
  write("done. ")


  -- NEW: init projections
  -- ProjectVerticesToSphere (dom, {0.0, 0.0, 0.0}, 1.0, 0.1)
  
  self.refinementProjector = DomainRefinementProjectionHandler(dom)
  local falloffProjector = SphericalFalloffProjector(dom, {0.0, 0.0, 0.0}, 1.0, 0.1)
  self.refinementProjector:set_callback("SURF", falloffProjector)  


  -- create Refiner
  ug_assert(numPreRefs <= numRefs, "numPreRefs must be smaller than numRefs. Aborting.");
  
  if numPreRefs > numRefs then
    numPreRefs = numRefs
  end
  
  -- Create a refiner instance. This is a factory method
  -- which automatically creates a parallel refiner if required.
  local refiner = nil
  if numRefs > 0 then
    refiner = GlobalDomainRefiner(dom)
     -- NEW: set callback
   -- refiner:set_refinement_callback(self.refinementProjector)
  end
  
 
  
  
  -- Performing pre-refines
  write("Pre-Refining("..numPreRefs.."): ")
  for i=1,numPreRefs do
    TerminateAbortedRun()
    write(i .. " ")
    refiner:refine()
  end
  write("done. Distributing...")
  
  
  -- Distribute the domain to all involved processes
  if util.DistributeDomain(dom, distributionMethod, verticalInterfaces, numTargetProcs, distributionLevel, wFct) == false then
    ug_error("Error while Distributing Grid. Aborting.")
  end
  
  write(" done. Post-Refining("..(numRefs-numPreRefs).."): ")
  if numRefs > 0 then
    -- Perform post-refine
    for i=numPreRefs+1,numRefs do
      TerminateAbortedRun()
      refiner:refine()
      write(i-numPreRefs .. " ")
    end
  end
  write("done.\n")
  
  -- Now we loop all subsets an search for it in the SubsetHandler of the domain
  if neededSubsets ~= nil then
    if util.CheckSubsets(dom, neededSubsets) == false then 
      ug_error("Something wrong with required subsets. Aborting.");
    end
  end
  
  
  --clean up
  if refiner ~= nil then
    delete(refiner)
  end
  
  
  SaveGridHierarchyTransformed(dom:grid(), dom:subset_handler(), "sphere3d_mg.ugx", 3)
  
  -- return the created domain
  return dom

end


function cryer3d:add_elem_discs(domainDisc, bStationary)
  CryerAddElemDiscs(self, domainDisc, bStationary)
end


-- boundary conditions
function cryer3d:add_boundary_conditions(domainDisc, bStationary)
  
  
 local doStationary = bStationary or false

 local neumannX = NeumannBoundaryFE("ux")
 local neumannY = NeumannBoundaryFE("uy")
 local neumannZ = NeumannBoundaryFE("uz")
 
 
 --neumannX:add("Cryer3dBndX", "RIM", "INNER")
 --neumannY:add("Cryer3dBndY", "RIM", "INNER")
 
 if doStationary then 
    neumannX:set_stationary()
    neumannY:set_stationary() 
    neumannZ:set_stationary() 
 end
 
 --domainDisc:add(neumannX)
 --domainDisc:add(neumannY)
 --domainDisc:add(neumannZ)

  local dirichlet = DirichletBoundary()
  dirichlet:add(0.0, "p", "SURF")
  dirichlet:add(0.0, "ux", "CENTER, XY, XZ")
  dirichlet:add(0.0, "uy", "CENTER, XY, YZ")
  dirichlet:add(0.0, "uz", "CENTER, XZ, YZ")
  
  --dirichlet:add(0.0, "ux", "NORTH, SOUTH, RIM")
  --dirichlet:add(0.0, "uy", "EAST, WEST, RIM")

  domainDisc:add(dirichlet)
 
end

-- init all variables
function cryer3d:init(kperm, poro, nu, cmedium, cfluid, csolid, volumetricweight)
  volumetricweight = volumetricweight or 1.0 -- kg/m^3
  CryerSetModelParameters_(self, nu, poro, cmedium, cfluid, csolid, kperm*volumetricweight)
  CryerInitRoots_(self)
  self.charTime = self:get_char_time(); 
  print ("charTime="..self.charTime)
end

function cryer3d:get_char_time()
 return CryerCharTime(self)
end

-- initial values
function cryer3d:interpolate_start_values(u, startTime)
  Interpolate(self.modelParameter.p0, u, "p", startTime)
  Interpolate(0.0, u, "ux", startTime)
  Interpolate(0.0, u, "uy", startTime)
  Interpolate(0.0, u, "uz", startTime)
end

-- post processing (after each step)
function cryer3d:post_processing(u, step, time)
  print ("p0:\t"..time.."\t"..time/self.charTime.."\t"..(Integral(u, "p", "CENTER")/self.modelParameter.p0).."\t"..CryerComputePressure(self, 0.0,0.0,time))
end