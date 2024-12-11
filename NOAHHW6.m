clc;
clear;

% parameters
% step size
% grid points dimension
dt = 0.001;  
Da = 10;         
Dh = 200;     
N = 100;         
% Generate matrix for Laplacian operator
L = laplacian_2d(N,N);
% initial conditions
a = 0.1 + 0.9*rand(N,N);  
h = 0.1*ones(N,N);  
%time, counter plotS
t = 0;  
k = 0;  
%store
a_prev = a;
h_prev = h;


% Simulation
while t < 100
% Updating variables
a = a + dt*(xa(a_prev,h_prev) + Da*reshape((L*reshape(a,N*N,1)),N,N));
h = h + dt*(xb(a_prev,h_prev) + Dh*reshape((L*reshape(h,N*N,1)),N,N));
% add time and counter
t = t + dt;   
k = k + 1;   
    
if mod(k,100) == 0
% Plot
figure(1);
pcolor(a);
colormap('gray');
shading flat;
end
    
% update previous values
a_prev = a;
h_prev = h;
end
% c1=c2=1
function xa = xa(a, h)
%mue = 1;
c1 = 1;
K1 = 0.03;
xa = c1*a.^2./(h.*(1+K1*a.^2)) - a;
end

function xb = xb(a, h)
c2 = 1;
v = 1.2;
xb = c2*a.^2 - v*h;
end

%attached document listed on canvas
% Laplacian operator matrix
function L = laplacian_2d(nx, ny)
% 1D Laplacian matrices
Lx = laplacian_1d(nx);
Ly = laplacian_1d(ny);
    
% Neumann boundary conditions
Lx(1,1) = -1;
Lx(1,2) = 1;
Lx(nx,nx-1) = 1;
Lx(nx,nx) = -1;
Ly(1,1) = -1;
Ly(1,2) = 1;
Ly(ny,ny-1) = 1;
Ly(ny,ny) = -1;
    
% 2D Laplacian matrix 
Ix = speye(nx);
Iy = speye(ny);
L = kron(Iy, Lx) + kron(Ly, Ix);
end

% 1D Laplacian matrix
function L = laplacian_1d(n)
e = ones(n,1);
L = spdiags([e -2*e e], [-1 0 1], n, n);
end
