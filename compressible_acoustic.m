clc; close all; clear all;

% Script to solve 1-D compressible acoustic waves
% with a Hamiltonian finite volume scheme.

% We assume that the equations have been scaled such that 
% R_0 = c^2_0 = 1 
%         u_t = -r_x
%         p_t = -u_x

% Boundaries at x = a & x = L are taken to be solid walls,
% so no normal flow through the boundary. A similar condition is taken 
% for the test function, to preserve skew symmetry, which leads to a
% condition for the density at the boundary.


% Currently attempts to save animation as movies fail.
% Uncomment pause code in the timestep loop to view wave evolution.

% Need to add clause for Nx = 2 , 3 in matrix loops

%% Mesh Constants

a = 0; % Mesh starting point
L = 1; % Mesh end point
Nx = 16; % Number of elements
dx = (L-a)/Nx; % Element size
dt = dx; % Timestep discretisation
theta = 0.5; % Flux constant, 0 < theta < 1
periods = 100; % End time for simulation, as period = 1 non dim time
t = 0; % Starting time
%% Storage
% Vector for solution
% U = [U,R]^T
U = zeros(2*Nx,1);

% Discrete Energy
H = zeros(ceil(periods/dt),1);

% Storage for flux matrices
DIV = zeros(2*Nx, 2*Nx);

%% Initial Condition
X = linspace(a+dx, L-dx, Nx);

% For each element we use the cell centre of the initial condition
for i = 1:Nx
    U(i) = -sin(2*pi*X(i))*sin(2*pi*0.125);
    U(Nx+i) = -cos(2*pi*X(i))*cos(2*pi*0.125);
end 

% Initial Energy
H0 = sum(0.5*dx*(U(1:Nx).*U(1:Nx) + U(Nx+1:2*Nx).*U(Nx+1:2*Nx)));

%% Build flux matrices
% First we construct the system d/dt(U) = DIV(U),
% DIV(U) will be the discrete divergence .
% In this example DIV is the sum of the numerical line fluxes for
% each element.

% Internal cell contributions - this will form a tridiagonal system
for i = 2:Nx-1
    % Flux for velocity U, which depends on R
    DIV(i, Nx + i -1) =  theta ;
    DIV(i, Nx + i ) =  -theta + (1- theta) ;
    DIV(i, Nx + i +1) =  -(1 - theta) ;
    % Flux for density R, which depends on U
    DIV(Nx + i,  i -1) = 1 - theta ;
    DIV(Nx + i,  i ) =  theta - (1-theta);
    DIV(Nx + i,  i +1) = - theta;
end    

% Solid wall boundaries
% Left wall
DIV(1,Nx+1) =  (1-theta);
DIV(1,Nx+2) = - (1-theta);
DIV(Nx+1,1) = -(1-theta);
DIV(Nx+1,2) = -theta;
% Right wall
DIV(Nx,2*Nx-1) = theta;
DIV(Nx,2*Nx) = -theta;
DIV(2*Nx,Nx-1) = (1-theta);
DIV(2*Nx,Nx) = theta;

% Create timestepper matrices
% Using implicit mid-point timestepper creates matrices acting on 
% U^(n+1) and U^n.
% The system we will solve will be 
% P U^(n+1) = Q U^n
% P = I - 0.5*dt*DIV
% Q = I + 0.5*dt*DIV

P = eye(2*Nx) - dt*0.5*DIV;
Q = eye(2*Nx) + dt*0.5*DIV;

% Inverse of P\Q is not time dependant so we will build it here
Inverse = P\Q;

%F(ceil(periods/dt)) = struct('cdata',[],'colormap',[]);
%% Timestep Loop
count_energy = 0;
figure
while t < periods
    t =t+dt;
    % Calculate  U^(n+1) from  U^n and advance in time
    U = Inverse*U;
    % Calculate Energy
    count_energy = count_energy+1;
    H(count_energy,1) = sum(0.5*dx*(U(1:Nx).*U(1:Nx) + U(Nx+1:2*Nx).*U(Nx+1:2*Nx)));
    
    
    % Plot velocity 
    plot(X,U(1:Nx))
    title(['Solution after ', num2str(dt*count_energy), ' periods - Velocity u'])
    axis([  a L -1 1])
    xlabel('x')
    ylabel('u')
    %pause(0.01)
    hold off
    
    
    
end
%F(count_energy)=getframe(gcf);
%movie2avi(F, 'velocity.avi')

% Plot Error in Energy
figure
hold all
plot(linspace(0,periods,periods/dt), (H-H0)/H0)
title('Relative Error in Energy')
xlabel('Number of Periods')
hold off

