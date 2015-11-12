% This problem is taken from 'Numerical Mathematics and Computing' 6th
% edition by Ward Cheney and David Kincaid, Thomson Brooks/Cole 2008

% This code solves the heat diffusion equation for a thin 1-dimensional
% rod of length 1. The temperature, u(x,t), is solved for and plotted in
% time and space. The boundary conditions are u(x=0,t) = u(x=1,t) = 0.0
% and u(x,t=0) = sin(pi*x).
% The analytical solution is u(x,t) = exp(-pi^2 t) sin(pi*x). The rod
% exponentially decays to zero temperature.

clear;
h = 0.05; % step size in space variable
k = 0.005; % step size in time variable
xgrid = h:h:1-h;
tgrid = 0:k:0.5;
n = length(xgrid);
m = length(tgrid);
[X, T] = meshgrid(tgrid,xgrid);
u = zeros(n,m); % numerical solution
v = zeros(n,m); % analytical solution
b = zeros(n); % RHS of matrix equation
s = h^2/k; 
r = 2.0 + s;
% set up the matrix A to be inverted
% A is tridiagonal and diagonally dominant
A = zeros(n);
for i=2:n-1
    A(i,i) = r;
    A(i,i+1) = -1;
    A(i,i-1) = -1;
end
A(1,1) = r;
A(1,2) = -1;
A(n,n) = r;
A(n,n-1) = -1;

% initial conditions
u(:,1) = sin(pi*xgrid);

% loop over time. At each time step the system Ax=b is solved. The method
% is based on Crank-Nicolson.
for k = 1:m-1
    b = s*u(:,k);
    u(:,k+1) = A \ b;
end

% analytic solution
for j=1:m
   for i=1:n
       v(i,j) = exp(-pi^2*tgrid(j))*sin(pi*xgrid(i));
   end
end

figure
subplot(1,2,1)
surf(X,T,u)
title('Numerical solution','fontsize',24)
xlabel('Time','fontsize',24)
ylabel('Length','fontsize',24)
zlabel('Temperature','fontsize',24)
shading interp;
subplot(1,2,2)
surf(X,T,v)
title('Analytical solution','fontsize',24)
xlabel('Time','fontsize',24)
ylabel('Length','fontsize',24)
zlabel('Temperature','fontsize',24)
shading interp;

% the percent difference at each grid point can be plotted
% figure
% surf(X,T,abs(v-u)./v*100)
% title('difference')
% shading interp;
