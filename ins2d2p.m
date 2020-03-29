function ins2d2p
% 2D 2-Phase Incompressible Navier-Stokes Solver with Level Set
%
% Numerics: * Colocated Cartesian mesh in a square domain
%           * Projection method for velocity-pressure coupling
%           * Explicit 2nd order RK for time stepping
%           * 2nd order finite volume reconstruction with Superbee 
%             limiter for advection, WENO5 for reinitialization
%             and central scheme for the rest
%           * 2nd order viscosity terms are treated implicitly
%           * Multigrid V cycle is used to solve linear systems
%
% Conventions: * In: phi > 0; Out: phi < 0; Interface: phi = 0
%
% Ref: * Q. Wang, MIT Numerical Methods for PDE lecture, 2015. - MG
%      * B. Seibold, A compact and fast Matlab code ... , 2008. - Matlab
%      * G. Tryggvason, CFD Course Lecture Slides. - CFD theory
%      * A non-oscillatory Eulerian ... by Fedkiw, 1999. - WENO5
%      * http://www.geocomputation.org/1999/082/gc_082.htm - Bilinear int.
% -----------------------------------------------------------------------
n = 32;                       % Number of cells, need to be 2^n
global Lx; Lx = 0.5;          % L in x dir.
global scale; scale = 8;      % L in y dir./L in x dir., need to be 2^n
T = 1.;                       % End time
rhoIn = 0.1694; rhoOut = 1.225;   % Densities
muIn = 0.00313; muOut = 0.00313;  % Viscosity
gy = -9.8;                    % Gravitational acc.
sigC = 0.;                    % Surface tension coefficient
cEps = 1.0;                   % Constant for interface thickness
rein = 1; itr = 3; nr = 2;    % Rein. control parameters
cDtho = 0.5;                  % Constant for rein. time step
nImp = 2; nPPE = 2;           % Number of iter. for MG
global nSmooth; nSmooth = 2;  % Number of MG smoother steps
ig = 1e2;                     % Time step for dynamic post-process
dtmax = 1e-2;                 % Max time step
% -----------------------------------------------------------------------

% Setting the grid
x = linspace(0,Lx,n+1); 
y = linspace(0,Lx*scale,scale*n+1); 
xP = 0.5*(x(:,2:end)+x(:,1:end-1)); 
yP = 0.5*(y(:,2:end)+y(:,1:end-1));
h = Lx/n;  
[X,Y] = meshgrid(xP,yP);

% Reinitialization parameters
eps = cEps*h; 
dtho = cDtho*h;

% Initial conditions
U = zeros(n,scale*n);
V = zeros(n,scale*n);
P = zeros(n,scale*n);
uf = zeros(n+1,scale*n); 
vf = zeros(n,scale*n+1);

% Initial level set function
phi = 2-Y'-0.05*cos(2*pi*X');

% Boundary conditions
uW = yP*0; uE = yP*0;
vN = xP'*0; vS = xP'*0;

% Plot the initial level set function
contour(X,Y,phi',[0 0],'r','LineWidth',1);
set(gca,'xTick',0:0.2:Lx); 
set(gca,'yTick',0:0.2:scale*Lx); 
grid('on'); axis('equal'); pause;

dum1 = 1; dum2 = 1; dum3 = 1;
t = 0; it = 0;
tic; % --------------------------------------- START OF THE TIME STEPPING
while (true)

  % Time step based on CFL
  maxU = max(max(max(abs(uf))),max(max(abs(vf)))); 
  dt = min(0.1*h/maxU,dtmax);
  
  % If new time is bigger than T, adjust it.
  if ((t+dt)>T); dt = T-t; end
  
  % Adjust time step for post-process
  if (dum1 && ((t+dt)>0.7)); dt = 0.7-t; dum1 = 0; end
  if (dum2 && ((t+dt)>0.8)); dt = 0.8-t; dum2 = 0; end
  if (dum3 && ((t+dt)>0.9)); dt = 0.9-t; dum3 = 0; end
  
  % Old solution (used for RK2)
  Un = U;  Vn = V; phin = phi;
  
  for is = 1:2 % ----------------------------------- START OF THE RK LOOP
  
    % Extend u and v to include boundary points (2 ghost cells)
    Uex = [-3*U(1,:);-U(1,:);U;-U(end,:);-3*U(end,:)]; % W-E
    Uey = [-3*U(:,1) -U(:,1) U -U(:,end) -3*U(:,end)]; % S-N
    Vex = [V(1,:);V(1,:);V;V(end,:);V(end,:)]; % W-E
    Vey = [-3*V(:,1) -V(:,1) V -V(:,end) -3*V(:,end)]; % S-N  
    
    % Extend phi to include boundary points (2 ghost cells)
    phiex = [repmat(phi(1,:),2,1);phi;repmat(phi(end,:),2,1)]; % W-E
    phiey = [repmat(phi(:,1),1,2) phi repmat(phi(:,end),1,2)]; % S-N    
    
    % Smoothed Heaviside function
    H = 0.5*(1+tanh(0.5*phi/eps));
    
    % grad(H) - needed for surface tension term
    Hex = [H(1,:);H;H(end,:)];
    Hey = [H(:,1) H H(:,end)];
    HxC = 0.5*(Hex(3:end,:)-Hex(1:end-2,:))/h;
    HyC = 0.5*(Hey(:,3:end)-Hey(:,1:end-2))/h;

    % Calculate density and viscosity
    rho = H.*rhoIn+(1-H).*rhoOut;
    mu = H.*muIn+(1-H).*muOut;
       
    % Calculate interface normal and curvature: \kappa = -div(\vec(n)) 
    absGrad = sqrt(HxC.^2+HyC.^2);
    nx = HxC./(absGrad+1e-10);
    ny = HyC./(absGrad+1e-10);
    nxe = [nx(1,:);nx;nx(end,:)];
    nye = [ny(:,1) ny ny(:,end)];
    kappa = -(0.5/h)*(nxe(3:end,:)-nxe(1:end-2,:)+...
                      nye(:,3:end)-nye(:,1:end-2));
    
    % Face velocities 
    uf = [uW;0.5*(U(2:end,:)+U(1:end-1,:));uE];
    vf = [vS 0.5*(V(:,2:end)+V(:,1:end-1)) vN];  
  
    % U velocity gradients
    UX = diff(Uex); UY = diff(Uey,1,2);
    
    % Treat convective terms - U velocity 
    resU = fv2(Uex,Uey,UX,UY,uf,vf);
   
    % V velocity gradients
    VX = diff(Vex); VY = diff(Vey,1,2);
    
    % Treat convective terms - V velocity 
    resV = fv2(Vex,Vey,VX,VY,uf,vf);

    % Treat body forces
    fx = (sigC*kappa.*HxC)./rho;
    fy = gy+(sigC*kappa.*HyC)./rho;
    
    % Treat explicit part of viscous term
    muxe = [mu(1,:);mu;mu(end,:)];
    muye = [mu(:,1) mu mu(:,end)];
    mux = muxe(3:end,:)-muxe(1:end-2,:);
    muy = muye(:,3:end)-muye(:,1:end-2);
    UxC = Uex(4:end-1,:)-Uex(2:end-3,:);
    UyC = Uey(:,4:end-1)-Uey(:,2:end-3);
    VxC = Vex(4:end-1,:)-Vex(2:end-3,:);
    VyC = Vey(:,4:end-1)-Vey(:,2:end-3);
    
    visx = ((0.5./rho).*mux.*UxC+(0.25./rho).*muy.*(UyC+VxC))/(h^2);
    visy = ((0.5./rho).*muy.*VyC+(0.25./rho).*mux.*(UyC+VxC))/(h^2);
    
    % Update velocities with convective, viscous and body terms
    U = U+(dt/h)*resU+dt*fx+dt*visx;
    V = V+(dt/h)*resV+dt*fy+dt*visy;  
    
    % Treat viscosity terms (2nd order part) implicitly with MG V-cycle
    Uo = U; Vo = V;
    for k = 1:nImp 
      U = mgVcycleU(U,Uo,dt*mu./rho);
      V = mgVcycleV(V,Vo,dt*mu./rho);
    end
    
    % Face velocities with updated velocities
    uf = [uW;0.5*(U(2:end,:)+U(1:end-1,:));uE];
    vf = [vS 0.5*(V(:,2:end)+V(:,1:end-1)) vN];
    
    % Divergence of velocity (RHS of PPE)
    div = (diff(uf,1,1)+diff(vf,1,2))/h/dt; 
  
    % Pressure Poisson Equation (PPE) with MG V-cycle
    for k = 1:nPPE
      P = mgVcycle(P,div,1./rho); 
    end
    
    % Pressure Correction
    Pex = [P(1,:);P;P(end,:)];
    Pey = [P(:,1) P P(:,end)];
  
    U = U-(0.5*dt/h)*(Pex(3:end,:)-Pex(1:end-2,:))./rho;
    V = V-(0.5*dt/h)*(Pey(:,3:end)-Pey(:,1:end-2))./rho; 
    
    % Face velocities with updated velocities
    uf = [uW;0.5*(U(2:end,:)+U(1:end-1,:));uE];
    vf = [vS 0.5*(V(:,2:end)+V(:,1:end-1)) vN];   
    
    % Phi gradients
    phiX = diff(phiex); phiY = diff(phiey,1,2);
    
    % Level set advection
    resPhi = fv2(phiex,phiey,phiX,phiY,uf,vf);
    
    % Update phi
    phi = phi+(dt/h)*resPhi;
    
  end % ---------------------------------------------- END OF THE RK LOOP    

  % Average for RK2
  U = 0.5*(U+Un); V = 0.5*(V+Vn); phi = 0.5*(phi+phin);
  
  if (rein && floor(it/itr) == it/itr) % ---------- REINITIALIZATION LOOP

    % Smoothed sign function
    S = phi./sqrt(phi.^2+h^2);
    s = sign(S);

    for ir = 1:nr % ------------------------- # of REINITIALIZATION STEPS

      % Old solution for RK2
      phin = phi;

      for isr = 1:2 % ------------------------------------------ RK2 LOOP

        % Ghost cell values (No change in phi)
        phiex = [repmat(phi(1,:),3,1);phi;repmat(phi(end,:),3,1)];
        phiey = [repmat(phi(:,1),1,3) phi repmat(phi(:,end),1,3)];

        % Phi gradients
        phiX = diff(phiex)/h; phiY = diff(phiey,1,2)/h;

        % Right and left states from WENO5
        [dxM,dxP,dyM,dyP] = weno5(phiX,phiY);

        % RHS of rein. equation
        dx2 = max(max(s.*dxM,-s.*dxP),0).^2;
        dy2 = max(max(s.*dyM,-s.*dyP),0).^2;
        absGrad = sqrt(dx2+dy2);
        rRhs = S.*(1-absGrad);

        % Update phi
        phi = phi+dtho*rRhs;
        
      end % --------------------------------------------- END OF RK2 LOOP

      % Average for RK2
      phi = 0.5*(phi+phin);

    end  % --------------------------- END OF # of REINITIALIZATION STEPS 

  end % ------------------------------------ END OF REINITIALIZATION LOOP  
  
  t = t+dt;
  it = it+1; disp(it)
  
  % Dynamic post-process
  if it == 1 || floor(it/ig) == it/ig || abs(t-T)<1e-10
    clf;
    contour(X,Y,phi',[0 0],'r','LineWidth',1);
    set(gca,'xTick',0:0.2:Lx); 
    set(gca,'yTick',0:0.2:scale*Lx);
    axis('equal'); grid on; drawnow;
  end 
  
  % Calculate area for in phase
  area(it) = sum(sum(H*h*h)); 
  tit(it) = t;

  if (abs(t-T)<1e-10); break; end
 
end; ctime = toc; % ---------------------------  END OF THE TIME STEPPING  
disp(ctime)

% Save the final bubble shape
eval(['print -dpng ','n',num2str(n)]);

% Normalize area variation with initial area
area = area/area(1);

% Extra figures
figure(2)
contourf(X,Y,phi',20);
set(gca,'xTick',0:0.1:Lx); 
set(gca,'yTick',0:0.1:scale*Lx);
axis('equal'); grid on;

figure(3)
plot(tit,area,'LineWidth',3);

% =======================================================================
function [rc] = f2c(rf) 
% Fine to coarse residual interpolation (restriction) - Simple averaging
rc = 0.25*(rf(1:2:end-1,1:2:end-1)+rf(2:2:end,1:2:end-1)+...
           rf(1:2:end-1,2:2:end)+rf(2:2:end,2:2:end));

% =======================================================================
function [ef] = c2f(ec)
% Bilinear interpolation
global scale;

% Ghost cell values (No change in e on the boundary)
ece = [ec(1,:);ec;ec(end,:)];
ece = [ece(:,1) ece ece(:,end)];

n = size(ec,1);
ef = zeros(2*n,scale*2*n);

% Coarse to fine error interpolation (prolongation)
ef(1:2:end-1,1:2:end-1) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(1:end-2,2:end-1)+3*ece(2:end-1,1:end-2)+ece(1:end-2,1:end-2));
ef(2:2:end,1:2:end-1) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(3:end,2:end-1)+3*ece(2:end-1,1:end-2)+ece(3:end,1:end-2));
ef(1:2:end-1,2:2:end) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(1:end-2,2:end-1)+3*ece(2:end-1,3:end)+ece(1:end-2,3:end));
ef(2:2:end,2:2:end) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(3:end,2:end-1)+3*ece(2:end-1,3:end)+ece(3:end,3:end));
  
% =======================================================================
function [ef] = c2fU(ec)
% Bilinear interpolation
global scale;

% Ghost cell values (e=0 on the boundary)
ece = [-ec(1,:);ec;-ec(end,:)];
ece = [-ece(:,1) ece -ece(:,end)];

n = size(ec,1);
ef = zeros(2*n,scale*2*n);

% Coarse to fine error interpolation (prolongation)
ef(1:2:end-1,1:2:end-1) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(1:end-2,2:end-1)+3*ece(2:end-1,1:end-2)+ece(1:end-2,1:end-2));
ef(2:2:end,1:2:end-1) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(3:end,2:end-1)+3*ece(2:end-1,1:end-2)+ece(3:end,1:end-2));
ef(1:2:end-1,2:2:end) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(1:end-2,2:end-1)+3*ece(2:end-1,3:end)+ece(1:end-2,3:end));
ef(2:2:end,2:2:end) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(3:end,2:end-1)+3*ece(2:end-1,3:end)+ece(3:end,3:end));
  
% =======================================================================
function [ef] = c2fV(ec)
% Bilinear interpolation
global scale;

% Ghost cell values (e=0 on the boundary)
ece = [ec(1,:);ec;ec(end,:)];
ece = [-ece(:,1) ece -ece(:,end)];

n = size(ec,1);
ef = zeros(2*n,scale*2*n);

% Coarse to fine error interpolation (prolongation)
ef(1:2:end-1,1:2:end-1) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(1:end-2,2:end-1)+3*ece(2:end-1,1:end-2)+ece(1:end-2,1:end-2));
ef(2:2:end,1:2:end-1) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(3:end,2:end-1)+3*ece(2:end-1,1:end-2)+ece(3:end,1:end-2));
ef(1:2:end-1,2:2:end) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(1:end-2,2:end-1)+3*ece(2:end-1,3:end)+ece(1:end-2,3:end));
ef(2:2:end,2:2:end) = 0.0625*(9*ece(2:end-1,2:end-1)+...
    3*ece(3:end,2:end-1)+3*ece(2:end-1,3:end)+ece(3:end,3:end));
  
% =======================================================================
function [u] = Jacobi(u,b,k)
% Jacobi iteration for div(k grad(u)) = b
global Lx;

n = size(u,1);
h = Lx/n;

kfx = 0.5*(k(2:end,:)+k(1:end-1,:));
kfy = 0.5*(k(:,2:end)+k(:,1:end-1));

kfxe = [kfx(1,:);kfx;kfx(end,:)];
kfye = [kfy(:,1) kfy kfy(:,end)];

% Ghost cell values (No change in u)
ue = [u(1,:);u;u(end,:)];
ue = [ue(:,1) ue ue(:,end)];

% Jacobi iter.
u = (kfxe(2:end,:).*ue(3:end,2:end-1)+...
     kfye(:,2:end).*ue(2:end-1,3:end)+...
     kfxe(1:end-1,:).*ue(1:end-2,2:end-1)+...
     kfye(:,1:end-1).*ue(2:end-1,1:end-2)-b*h*h)./...
    (kfxe(2:end,:)+kfxe(1:end-1,:)+kfye(:,2:end)+kfye(:,1:end-1));
               
% =======================================================================
function [r] = residual(u,b,k) 
% Residual calc. for div(k grad(u)) = b
global Lx;

n = size(u,1);
h = Lx/n;

kfx = 0.5*(k(2:end,:)+k(1:end-1,:));
kfy = 0.5*(k(:,2:end)+k(:,1:end-1));

kfxe = [kfx(1,:);kfx;kfx(end,:)];
kfye = [kfy(:,1) kfy kfy(:,end)];

% Ghost cell values (No change in u)
ue = [u(1,:);u;u(end,:)];
ue = [ue(:,1) ue ue(:,end)];

r = b-(kfxe(2:end,:).*ue(3:end,2:end-1)+...
       kfye(:,2:end).*ue(2:end-1,3:end)+...
       kfxe(1:end-1,:).*ue(1:end-2,2:end-1)+...
       kfye(:,1:end-1).*ue(2:end-1,1:end-2)-...
      (kfxe(2:end,:)+kfxe(1:end-1,:)+kfye(:,2:end)+kfye(:,1:end-1))...
      .*ue(2:end-1,2:end-1))/(h^2);

% =======================================================================
function [u] = mgVcycle(u,b,kappa) 
% Multigrid V cycle for PPE
global nSmooth;
global scale;

n = size(u,1);

% Pre-smoothing
for i=1:nSmooth
    u = 0.8*Jacobi(u,b,kappa)+0.2*u;
end

if n>2
    r = residual(u,b,kappa);
    rc = f2c(r);
    kappac = f2c(kappa);
    ec = mgVcycle(zeros(n/2,scale*n/2),rc,kappac);
    e = c2f(ec);
    u = u+e;
end

% Post-smoothing
for i=1:nSmooth
    u = 0.8*Jacobi(u,b,kappa)+0.2*u;
end

% =======================================================================
function [u] = JacobiU(u,b,a)
% Jacobi iteration for u-a*Laplacian(u) = b
global Lx;

n = size(u,1);
h = Lx/n;

% Ghost cell values (u=0 on the boundary)
ue = [-u(1,:);u;-u(end,:)];
ue = [-ue(:,1) ue -ue(:,end)];    

% Jacobi iter.
u = (a.*(ue(3:end,2:end-1)+ue(2:end-1,3:end)+...
    ue(1:end-2,2:end-1)+ue(2:end-1,1:end-2))+b*h*h)./(h^2+4*a);

% =======================================================================
function [r] = residualU(u,b,a) 
% Residual calculation for u-a*Laplacian(u) = b
global Lx;

n = size(u,1);
h = Lx/n;

% Ghost cell values (u=0 on the boundary)
ue = [-u(1,:);u;-u(end,:)];
ue = [-ue(:,1) ue -ue(:,end)];    

r = b-u+a.*(ue(3:end,2:end-1)+ue(2:end-1,3:end)+ue(1:end-2,2:end-1)+ ...
    ue(2:end-1,1:end-2)-4*ue(2:end-1,2:end-1))/(h^2);

% =======================================================================
function [u] = mgVcycleU(u,b,a) 
% Multigrid V cycle for U velocity
n = size(u,1);
global nSmooth;
global scale;

% Pre-smoothing
for i=1:nSmooth
  u = 0.8*JacobiU(u,b,a)+0.2*u;
end

if n>2
  r = residualU(u,b,a);
  rc = f2c(r);
  ac = f2c(a);
  ec = mgVcycleU(zeros(n/2,scale*n/2),rc,ac);
  e = c2fU(ec);
  u = u+e;
end

% Post-smoothing
for i=1:nSmooth
  u = 0.8*JacobiU(u,b,a)+0.2*u;
end

% =======================================================================
function [u] = JacobiV(u,b,a)
% Jacobi iteration for u-a*Laplacian(u) = b
global Lx;

n = size(u,1);
h = Lx/n;

% Ghost cell values (u=0 on the boundary)
ue = [u(1,:);u;u(end,:)];
ue = [-ue(:,1) ue -ue(:,end)];
    
% Jacobi iter.
u = (a.*(ue(3:end,2:end-1)+ue(2:end-1,3:end)+...
    ue(1:end-2,2:end-1)+ue(2:end-1,1:end-2))+b*h*h)./(h^2+4*a);

% =======================================================================
function [r] = residualV(u,b,a) 
% Residual calculation for u-a*Laplacian(u) = b
global Lx;

n = size(u,1);
h = Lx/n;

% Ghost cell values (u=0 on the boundary)
ue = [u(1,:);u;u(end,:)];
ue = [-ue(:,1) ue -ue(:,end)];

r = b-u+a.*(ue(3:end,2:end-1)+ue(2:end-1,3:end)+ue(1:end-2,2:end-1)+ ...
    ue(2:end-1,1:end-2)-4*ue(2:end-1,2:end-1))/(h^2);
  
% =======================================================================
function [u] = mgVcycleV(u,b,a) 
% Multigrid V cycle for V velocity
n = size(u,1);
global nSmooth;
global scale;

% Pre-smoothing
for i=1:nSmooth
  u = 0.8*JacobiV(u,b,a)+0.2*u;
end

if n>2
  r = residualV(u,b,a);
  rc = f2c(r);
  ac = f2c(a);
  ec = mgVcycleV(zeros(n/2,scale*n/2),rc,ac);
  e = c2fV(ec);
  u = u+e;
end

% Post-smoothing
for i=1:nSmooth
  u = 0.8*JacobiV(u,b,a)+0.2*u;
end

% =======================================================================
function [resPhi] = fv2(phiex,phiey,phiX,phiY,uf,vf)
% 2nd order FV reconstruction w/ Superbee limiter

% Gradient ratios needed for limiters
rLx = phiX(2:end-1,:)./phiX(1:end-2,:);
rRx = phiX(3:end,:)./phiX(2:end-1,:);
rLy = phiY(:,2:end-1)./phiY(:,1:end-2);
rRy = phiY(:,3:end)./phiY(:,2:end-1);   

% Left and right reconstructed values
recLx = phiex(2:end-2,:)+...
        0.5*max(max(0,min(2*rLx,1)),min(rLx,2)).*phiX(1:end-2,:);
recRx = phiex(3:end-1,:)-...
        0.5*max(max(0,min(2*rRx,1)),min(rRx,2)).*phiX(2:end-1,:);
recLy = phiey(:,2:end-2)+...
        0.5*max(max(0,min(2*rLy,1)),min(rLy,2)).*phiY(:,1:end-2);
recRy = phiey(:,3:end-1)-...
        0.5*max(max(0,min(2*rRy,1)),min(rRy,2)).*phiY(:,2:end-1);

% Upwind flux
xFlux = 0.5.*uf.*(recLx+recRx)+0.5.*abs(uf).*(recLx-recRx);
yFlux = 0.5.*vf.*(recLy+recRy)+0.5.*abs(vf).*(recLy-recRy);

% Residual for phi
resPhi = -diff(xFlux)-diff(yFlux,1,2);

% =======================================================================
function [dxM,dxP,dyM,dyP] = weno5(vx,vy)
% WENO5 scheme

v1x = vx(1:end-5,:); v2x = vx(2:end-4,:);
v3x = vx(3:end-3,:); v4x = vx(4:end-2,:);
v5x = vx(5:end-1,:); v6x = vx(6:end,:);

v1y = vy(:,1:end-5); v2y = vy(:,2:end-4);
v3y = vy(:,3:end-3); v4y = vy(:,4:end-2);
v5y = vy(:,5:end-1); v6y = vy(:,6:end);

S1mx = (13/12)*((v1x-2*v2x+v3x).^2)+0.25*((v1x-4*v2x+3*v3x).^2);
S2mx = (13/12)*((v2x-2*v3x+v4x).^2)+0.25*((v2x-v4x).^2);
S3mx = (13/12)*((v3x-2*v4x+v5x).^2)+0.25*((3*v3x-4*v4x+v5x).^2);
S1px = (13/12)*((v6x-2*v5x+v4x).^2)+0.25*((v6x-4*v5x+3*v4x).^2);
S2px = (13/12)*((v5x-2*v4x+v3x).^2)+0.25*((v5x-v3x).^2);
S3px = (13/12)*((v4x-2*v3x+v2x).^2)+0.25*((3*v4x-4*v3x+v2x).^2);
S1my = (13/12)*((v1y-2*v2y+v3y).^2)+0.25*((v1y-4*v2y+3*v3y).^2);
S2my = (13/12)*((v2y-2*v3y+v4y).^2)+0.25*((v2y-v4y).^2);
S3my = (13/12)*((v3y-2*v4y+v5y).^2)+0.25*((3*v3y-4*v4y+v5y).^2);
S1py = (13/12)*((v6y-2*v5y+v4y).^2)+0.25*((v6y-4*v5y+3*v4y).^2);
S2py = (13/12)*((v5y-2*v4y+v3y).^2)+0.25*((v5y-v3y).^2);
S3py = (13/12)*((v4y-2*v3y+v2y).^2)+0.25*((3*v4y-4*v3y+v2y).^2);

a1mx = (1/10)./(((1e-6)+S1mx).^2);
a2mx = (6/10)./(((1e-6)+S2mx).^2);
a3mx = (3/10)./(((1e-6)+S3mx).^2);
a1px = (1/10)./(((1e-6)+S1px).^2);
a2px = (6/10)./(((1e-6)+S2px).^2);
a3px = (3/10)./(((1e-6)+S3px).^2);
a1my = (1/10)./(((1e-6)+S1my).^2);
a2my = (6/10)./(((1e-6)+S2my).^2);
a3my = (3/10)./(((1e-6)+S3my).^2);
a1py = (1/10)./(((1e-6)+S1py).^2);
a2py = (6/10)./(((1e-6)+S2py).^2);
a3py = (3/10)./(((1e-6)+S3py).^2);

w1mx = a1mx./(a1mx+a2mx+a3mx);
w2mx = a2mx./(a1mx+a2mx+a3mx);
w3mx = a3mx./(a1mx+a2mx+a3mx);
w1px = a1px./(a1px+a2px+a3px);
w2px = a2px./(a1px+a2px+a3px);
w3px = a3px./(a1px+a2px+a3px);
w1my = a1my./(a1my+a2my+a3my);
w2my = a2my./(a1my+a2my+a3my);
w3my = a3my./(a1my+a2my+a3my);
w1py = a1py./(a1py+a2py+a3py);
w2py = a2py./(a1py+a2py+a3py);
w3py = a3py./(a1py+a2py+a3py);

dxM = w1mx.*((1/3)*v1x-(7/6)*v2x+(11/6)*v3x)+...
    w2mx.*((-1/6)*v2x+(5/6)*v3x+(1/3)*v4x)+...
    w3mx.*((1/3)*v3x+(5/6)*v4x-(1/6)*v5x);
dxP = w1px.*((1/3)*v6x-(7/6)*v5x+(11/6)*v4x)+...
    w2px.*((-1/6)*v5x+(5/6)*v4x+(1/3)*v3x)+...
    w3px.*((1/3)*v4x+(5/6)*v3x-(1/6)*v2x);
dyM = w1my.*((1/3)*v1y-(7/6)*v2y+(11/6)*v3y)+...
    w2my.*((-1/6)*v2y+(5/6)*v3y+(1/3)*v4y)+...
    w3my.*((1/3)*v3y+(5/6)*v4y-(1/6)*v5y);
dyP = w1py.*((1/3)*v6y-(7/6)*v5y+(11/6)*v4y)+...
    w2py.*((-1/6)*v5y+(5/6)*v4y+(1/3)*v3y)+...
    w3py.*((1/3)*v4y+(5/6)*v3y-(1/6)*v2y);



