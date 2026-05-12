% \infty+m parameters, mu and psi only depend on m=2
lamy = @(x,y) ones(size(x));
mu1 = @(x) 2*ones(size(x));
mu2 = @(x) 1*ones(size(x));
sigy = @(x,y,h) x.^3.*(x+1).*(y-0.5).*h.*(h-1);
Wy1 = @(x,y) 3*(y-0.5).*ones(size(x));
Wy2 = @(x,y) 2*(y-0.5).*ones(size(x));
thy1 = @(x,y) -3*y.*(y-1).*ones(size(x));
thy2 = @(x,y) -2*y.*(y-1).*ones(size(x));
qy1 = @(y) 8*(y-0.5);
qy2 = @(y) -8*(y-2);
ry1 = @(y) cos(2*pi*y);
ry2 = @(y) 2*y.*(y+5);
psi12 = @(x) 0*ones(size(x));
psi21 = @(x) 0*ones(size(x));

% discretization in x
xg = 128; % grid size for x (minus one end point)
tmp = linspace(0,1,xg+1); % auxiliary grid
x0 = tmp(1:xg); x1 = tmp(2:xg+1); % grids for u and v
x0h = x0(1:xg/2); xh0 = x0(xg/2+1:xg); % split grid for x0
x1h = x1(1:xg/2); xh1 = x1(xg/2+1:xg); % split grid for x1

% finite-difference approximation for observer from continuum parameters
% emulate continuum observer as 60+2 system
nc = 60; yc = linspace(1/nc,1,nc); % y-grid based on nc
[Ac,Bc,Cc,BQc,CRc] = n2FD(lamy,mu1,mu2,sigy,Wy1,Wy2,thy1,thy2,...
  qy1,qy2,ry1,ry2,psi12,psi21,x0,x1,xg,nc);

% continuum observer kernels (zeros ignored)
N22 = @(x,xi) exp(x-xi);
N21p2 = @(x,xi) exp(x-0.5*xi-1);
M1p1 = @(x,xi,y) Wy1(x,y)/3; 
M1p2 = @(x,xi,y) exp(x-0.5*xi-1).*Wy1(x,y)/3; 
M2 = @(x,xi,y) exp(x-xi).*Wy2(x,y)/2;
% vectorize M, sample to nc components
M1p1nc = zeros(xg*nc, 1); M1p2nc = M1p1nc; M2nc = M1p1nc;
for k = 1:nc
  M1p1nc((k-1)*xg+1:(k-1)*xg+xg/2) = M1p1(x1h, 0, yc(k))';
  M1p2nc((k-1)*xg+xg/2+1:k*xg) = M1p2(xh1, 0, yc(k))';
  M2nc((k-1)*xg+1:k*xg) = M2(x1, 0, yc(k))';
end
% output injection gains based on continuum kernels
Pc = [(M1p1nc+M1p2nc)*mu1(0), M2nc*mu2(0); zeros(xg,2); ...
  [[zeros(xg/2,1); mu1(0)*N21p2(xh0, 0)'], mu2(0)*N22(x0, 0)']];

% continuum control kernels (zeros ignored)
L22 = @(x,xi) -2*exp(2*(x-xi));
L12p2 = @(x,xi) -2*exp(x-2*xi); 
K1p1 = @(x,xi,y) -thy1(x,y)/3;
K1p2 = @(x,xi,y) -exp(x-2*xi).*thy1(x,y)/3;
K2 = @(x,xi,y) -exp(2*(x-xi)).*thy2(x,y)/2;
% vectorize K, discretize to nc components
K1p1nc = zeros(xg*nc, 1); K1p2nc = K1p1nc; K2nc = K1p1nc;
for k = 1:nc
  K1p1nc((k-1)*xg+xg/2+1:k*xg) = K1p1(1, xh1, yc(k))';
  K1p2nc((k-1)*xg+1:(k-1)*xg+xg/2) = K1p2(1, x1h, yc(k))';
  K2nc((k-1)*xg+1:k*xg) = K2(1, x1, yc(k))';
end

% initialize for-loop
nv = 53:2:59; % considered values of n
nnm = numel(nv);
usolc = cell(1,nnm);
Usolfc = cell(1,nnm);
Yerrc = cell(1,nnm);

for nn = 1:nnm
% finite-difference approximation for n+2 system from continuum parameters
n = nv(nn);
[A,B,C,BQ,CR] = n2FD(lamy,mu1,mu2,sigy,Wy1,Wy2,thy1,thy2,...
  qy1,qy2,ry1,ry2,psi12,psi21,x0,x1,xg,n);

% closed-loop system of plant and observer, and control law
Ae = [A+B*CR+BQ*C, -B*CRc; Pc*C+BQc*C, Ac-Pc*Cc];
Be = [B; Bc];
% backstepping controller, approximate integrals with trapezoidal rule
Ue1 = @(z) (trapz(K1p1nc.*z((n+2)*xg+1:(n+2)*xg+nc*xg))/nc + ...
  trapz(K1p2nc.*z((n+2)*xg+1:(n+2)*xg+nc*xg))/nc + ...
  trapz(L12p2(1,x0h').*z((n+nc+3)*xg+1:(n+nc+3)*xg+xg/2)))/xg;
Ue2 = @(z) (trapz(K2nc.*z((n+2)*xg+1:(n+2)*xg+nc*xg))/nc + ...
  trapz(L22(1,x0').*z((n+nc+3)*xg+1:(n+nc+4)*xg)))/xg;

% simulate for t \in [0,15]
opts = odeset('Jacobian', Ae); % pass Jacobian to ODE solver
T = 15; % simulation end time
z0 = [repmat(.5*sin(2*pi*x1)',n,1); repmat(.5*sin(2*pi*x0)',2,1); ...
  zeros((nc+2)*xg,1)]; % initial conditions for plant and observer
% solve closed-loop ODE
sol = ode45(@(t,z) Ae*z + Be*[Ue1(z); Ue2(z)], [0, T], z0, opts);
% create time grid and evaluate solution at the grid points
tg = 513; TT = linspace(0,T,tg); 
usol = deval(sol, TT); 
% evaluate backstepping controller at the time grid points
Usol1 = zeros(1,tg);
Usol2 = Usol1;
for k=1:tg
  Usol1(k) = Ue1(usol(:,k));
  Usol2(k) = Ue2(usol(:,k));
end
% full control input including the -Ru(1) term
Usolf = [Usol1; Usol2] + [zeros(2,(n+2)*xg), -CRc]*usol;
% output estimation error
Yerr = [-C Cc]*usol;
% store results
usolc{nn} = usol;
Usolfc{nn} = Usolf;
Yerrc{nn} = Yerr;
end

% plot observer error norms
noe = zeros(nnm,tg);
uosic = cell(1,4);
for nn = 1:nnm
  usol = usolc{nn};
  n = nv(nn);
  idata = zeros(xg,n);
  uosi = zeros(n*xg,tg);
  y = linspace(1/n,1,n);
  % sample (interpolate in y) observer u state to a compatible u^n state
  for kk = 1:tg
    data = reshape(usol((n+2)*xg+1:(n+2+nc)*xg,kk),xg,nc);
    for ll = 1:xg
      idata(ll,:) = interp1(yc, data(ll,:), y);
    end
    uosi(:,kk) = idata(:);
  end
  uosic{nn} = uosi;
  % compute observer error norm
  noe(nn,:) = sqrt(sum((usol(1:n*xg,:)-uosi).^2/(n*xg),1)) + ...
     sqrt(sum((usol(n*xg+1:(n+2)*xg,:) - ...
     usol((n+nc+2)*xg+1:(n+nc+4)*xg,:)).^2/(2*xg),1));
end
figure(1)
plot(TT, noe,'linewidth',2)
set(gca,'tickdir', 'out', 'fontsize',11,'xtick',0:5:15,'ytick',0:2:8)
set(gca,'position',get(gca,'position')+[0 0 0.04 0])
xlabel('$t$', 'interpreter', 'latex', 'fontsize',12)
ylabel('$\|(\tilde{\mathbf{u}}(t),\tilde{\mathbf{v}}(t))\|_E$',...
  'interpreter','latex','fontsize',12)
legend({'$n=53$','$n=55$','$n=57$','$n=59$'}, 'interpreter', 'latex',...
  'fontsize',12,'location','northeast','numcolumns',2)

% plot inputs 
figure(2)
subplot(211)
for nn = 1:nnm
Usolf2 = Usolfc{nn};
hold on
plot(TT, Usolf2(1,:), 'linewidth',2,'DisplayName', ...
  ['$U^1_{',num2str(nv(nn)),'}(t)$'])
end
set(gca,'tickdir', 'out', 'fontsize',11,'xtick',0:5:15,'xticklabel',[])
set(gca,'position',get(gca,'position')+[0.01 0.0 0.07 0.05])
ylabel('$U^1(t)$','interpreter','latex','fontsize',12, ...
  'rotation',0)
xlim([0 15])
legend('interpreter', 'latex','fontsize',12,'location','northeast', ...
  'numcolumns',2)
subplot(212)
for nn = 1:nnm
Usolf2 = Usolfc{nn};
hold on
plot(TT, Usolf2(2,:), 'linewidth',2,'DisplayName', ...
  ['$U^2_{',num2str(nv(nn)),'}(t)$']);
end
set(gca,'tickdir', 'out', 'fontsize',11,'xtick',0:5:15)
xlabel('$t$', 'interpreter', 'latex', 'fontsize',12)
set(gca,'position',get(gca,'position')+[0.01 0 0.07 0.05])
ylabel('$U^2(t)$','interpreter','latex','fontsize',12, ...
   'rotation',0)
legend('interpreter', 'latex','fontsize',12,'location','northeast', ...
  'numcolumns',2)
xlim([0 15])

% auxiliary function for finite-difference approximation of n+m systems
function [A,B,C,BQ,CR] = n2FD(lamy,mu1,mu2,sigy,Wy1,Wy2,thy1,thy2,...
  qy1,qy2,ry1,ry2,psi12,psi21,x0,x1,xg,n)

% sample n+2 parameters from continuum
y = linspace(1/n,1,n); % y-grid based on n
lam = @(i,x) lamy(x,y(i));
sig = @(i,j,x) sigy(x,y(i),y(j));
W1 = @(i,x) Wy1(x,y(i));
W2 = @(i,x) Wy2(x,y(i));
th1 = @(j,x) thy1(x,y(j));
th2 = @(j,x) thy2(x,y(j));
q1 = qy1(y);
q2 = qy2(y);
r1 = ry1(y);
r2 = ry2(y);

% finite difference approximation \dot x = Ax+BU, Y=Cx, for n+2 system
D = eye(xg) - diag(ones(1,xg-1), -1); % backward difference in x
A = zeros((n+2)*xg);
B = zeros((n+2)*xg, 2);
C = zeros(2,(n+2)*xg);
B((n+1)*xg,1) = mu1(1)*xg; % control operator, v1(1) = U1
B((n+2)*xg,2) = mu2(1)*xg; % control operator, v2(1) = U2
C(1,n*xg+1) = 1; C(2,(n+1)*xg+1) = 1; % output operator, Y = v(0)
% auxiliary input and output operators for boundary couplings
BQ = zeros((n+2)*xg,2); CR = zeros(2,(n+2)*xg);
BQ((0:n-1)*xg+1,1) = xg*lam(1:n,0).*q1;
BQ((0:n-1)*xg+1,2) = xg*lam(1:n,0).*q2;
CR(1,(1:n)*xg) = r1/n;
CR(2,(1:n)*xg) = r2/n;

% fill A, ignore the boundary conditions as those are accounted for by the
% coupling terms BQ+C and B*CR corresponding to u(0)=Qv(0) and v(1)=Ru(1),
% respectively
for k = 1:n
  Ik = (k-1)*xg+1:k*xg; % index set
  A(Ik,Ik) = -xg*lam(k,x1).*D; % tranport term
  for l = 1:n
    Il = (l-1)*xg+1:l*xg; % index set
    A(Ik,Il) = A(Ik,Il) + diag(sig(k,l,x1))/n; % sigma terms (scale 1/n)
  end
  Il = n*xg+1:(n+1)*xg; % index set
  A(Ik,Il) = A(Ik,Il) + diag(W1(k,x0)); % W1 term
  Il = (n+1)*xg+1:(n+2)*xg; % index set
  A(Ik,Il) = A(Ik,Il) + diag(W2(k,x0)); % W2 term
end
k=n+1;
Ik = (k-1)*xg+1:k*xg; % index set
A(Ik,Ik) = -xg*mu1(x0).*D'; % transport term
for l = 1:n
  Il = (l-1)*xg+1:l*xg; % index set
  A(Ik,Il) = A(Ik,Il) + diag(th1(l,x1))/n; % theta1 terms (scale 1/n)
end
k=n+2;
Ik = (k-1)*xg+1:k*xg; % index set
A(Ik,Ik) = -xg*mu2(x0).*D'; % transport term
for l = 1:n
  Il = (l-1)*xg+1:l*xg; % index set
  A(Ik,Il) = A(Ik,Il) + diag(th2(l,x1))/n; % theta2 terms (scale 1/n)
end
Il = xg*n+1:xg*(n+1); % index set
A(Ik,Il) = A(Ik,Il) + diag(psi21(x0)); % psi21 term
A(Il,Ik) = A(Il,Ik) + diag(psi12(x0)); % psi12 term
end