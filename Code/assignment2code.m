% Ahmad Osman
% 101070948
%% QUESTION 1
% Solution to simple 2D Laplace case

W = 2;
L = 3;
V0 = 1;

dx = 0.2; % x mesh spacing
dy = 0.2; % y mesh spacing
nx = L/dx; % Number of points along the x axis
ny = W/dy; % Number of points along the y axis

%% 
% The finite difference can be implemented using a matrix $GF=F$.
% V = voltage at discrete points
% F = force
% G = relation between voltages at different points

c1 = -2*(1/dx^2 + 1/dy^2);
c2 = 1/(dx^2);
c3 = 1/(dy^2);

%% 
% The mapCoordinate function takes an x,y coordinate and converts it to an index in the V array.

G = zeros(nx*ny,nx*ny);

for x=2:(nx-1)
    for y=2:(ny-1)
        i = coordinate(x,y,nx);
        G(i,i) = c1;
        G(i,coordinate(x-1,y,nx)) = c2;
        G(i,coordinate(x+1,y,nx)) = c2;
        G(i,coordinate(x,y-1,nx)) = c3;
        G(i,coordinate(x,y+1,nx)) = c3;
    end
end

%% 
% The F matrix is generated with boundary set to V0, where x = 0 and x = L.

F = zeros(nx*ny,1);

for y=1:ny
    i = coordinate(1,y,nx);
    G(i,i) = 1;
    
    F(i) = V0;
    
    i = coordinate(nx,y,nx);
    G(i,i) = 1;
end

%% 
% Set up the boundary conditions for analytical solution @ V=0 located at the
% corners.

for x=2:(nx-1)
    i = coordinate(x,1,nx);
    G(i,i) = 1;
    G(i,coordinate(x,2,nx)) = -1;
    
    i = coordinate(x,ny,nx);
    G(i,i) = 1;
    G(i,coordinate(x,ny-1,nx)) = -1;
end

%% 
% Solution from matrices

sol = G\F;
sol = reshape(sol,[],ny)';

figure(1);
surf(linspace(0,L,nx),linspace(0,W,ny),sol);
xlabel('x');
ylabel('y');
title(sprintf('Finite Differences Solution with Grid Spacing of %.2f', dx));
set(gca, 'View', [45 45])

set(gca,'Color', [0 0 0]);
a_x = gca; 
a_x.GridAlpha = 0.5;
a_x.GridColor = [1, 1, 1];


%% 
% Analytical solution for comparison with Finite Difference solution.

analyticalSol = zeros(ny, nx);
x1 = repmat(linspace(-L/2,L/2,nx),ny,1);
y1 = repmat(linspace(0,W,ny),nx,1)';
iter = 100;
avgError = zeros(iter,1);


for i=1:iter
    n = 2*i - 1;
    analyticalSol = analyticalSol + 1./n.*cosh(n.*pi.*x1./W) ...
        ./cosh(n.*pi.*(L./2)./W).*sin(n.*pi.*y1./W);

    avgError(i) = mean(mean(abs(analyticalSol.*4.*V0./pi - sol)));
end

analyticalSol = analyticalSol.*4.*V0./pi;

figure(2);
surf(linspace(0,L,nx),linspace(0,W,ny),analyticalSol);
xlabel('x');
ylabel('y');
title(sprintf('Analytical Solution with %d iterations', iter));

set(gca,'Color', [0 0 0]);
a_x = gca; 
a_x.GridAlpha = 0.5;
a_x.GridColor = [1, 1, 1];

%% QUESTION 2
% Solving current flow in a rectangular region using Finite Diffrenece
% Method.
%
% Part A
nx = 75;
ny = 50;
Lb = 20;
Wb = 10;
V1 = 1; 
figure(4);
hold on;

% Generating the map of conductivity of the area
sigma_conduct = 1;
sigma_insulate = 10e-2;
cMap = sigma_conduct*ones(nx, ny);
cMap(1:Wb,(1:Lb)+ny/2-Lb/2) = sigma_insulate;
cMap((1:Wb)+nx-Wb,(1:Lb)+ny/2-Lb/2) = sigma_insulate;
surf(linspace(0,1.5,ny), linspace(0,1,nx), cMap,'EdgeColor','none','LineStyle','none');
title("AO_101070948");
xlabel('x');
ylabel('y');
zlabel('Conduction (Mho)');
view([120 25])


% Numeric solution
V = numericSolution(nx, ny, cMap, Inf, Inf, 0, V1);
figure(5);
hold on;
surf(linspace(0,1.5,ny), linspace(0,1,nx), V,'EdgeColor','none','LineStyle','none');
title("AO_101070948");
xlabel('x');
ylabel('y');
zlabel('Voltage (V)');
view([120 25])
colorbar

% Electric field
[Ex, Ey] = gradient(V);
Ex = -Ex;
Ey = -Ey;
figure(6);
quiver(linspace(0,1.5,ny), linspace(0,1,nx), Ex, Ey);
ylim([0 1]);
xlim([0 1.5]);
xlabel('x');
ylabel('y');

% Current density
Jx = cMap.*Ex;
Jy = cMap.*Ey;
J = sqrt(Jx.^2 + Jy.^2);
figure(7);
hold on;
contourf(linspace(0,1.5,ny), linspace(0,1,nx), J,'EdgeColor','none','LineStyle','none');
quiver(linspace(0,1.5,ny), linspace(0,1,nx), Jx, Jy);
title("AO_101070948");
xlabel('x');
ylabel('y');
colorbar

set(gca,'Color', [0 0 0]);
a_x = gca; 
a_x.GridAlpha = 0.5;
a_x.GridColor = [1, 1, 1];

%%
% Part B
figure(8);
hold on;
range = 20:5:100;
I = [];
for x = range
    I = [I totalI(x, ny, V1, sigma_conduct, sigma_insulate, Wb, Lb)];
end
plot(range, I);
title("AO_101070948");
ylabel('Total Current (A)');
xlabel('Width mesh size');

set(gca,'Color', [0 0 0]);
a_x = gca; 
a_x.GridAlpha = 0.5;
a_x.GridColor = [1, 1, 1];

%%
% Part C
figure(9);
range = 0:1:50;
I = [];
for W = range
    I = [I totalI(nx, ny, V1, sigma_conduct, sigma_insulate, W, Lb)];
end
plot(range, I);
title("AO_101070948");
ylabel('Total Current (A)');
xlabel('Box width');

set(gca,'Color', [0 0 0]);
a_x = gca; 
a_x.GridAlpha = 0.5;
a_x.GridColor = [1, 1, 1];

%%
% Part D
figure(10);
hold on;
range = logspace(-5,0, 50);
I = [];
for sigma = range
    I = [I totalI(nx, ny, V1, sigma_conduct, sigma, Wb, Lb)];
end
plot(range, I);
title("AO_101070948");
ylabel('Total Current (A)');
xlabel('Box Conduction (Mho)');

set(gca,'Color', [0 0 0]);
a_x = gca; 
a_x.GridAlpha = 0.5;
a_x.GridColor = [1, 1, 1];

%% Functions
% the following functions set are used to operate the code for questions 1
% and 2


% The following functions is used to map the x and y coordinates
function i = coordinate(x,y,nx)
i = (y-1).*nx + x;

end
%%
% The following function numericSolution is used to calculate numerical
% solutions using the finite difference method
function V = numericSolution(nx, ny, crecCoordinates, bc_left, bc_right, bc_top, bc_bottom)

    global C;
    G = sparse(nx*ny, ny*nx);
    B = zeros(1, nx*ny);
    for i=1:nx
        for j=1:ny
            n = recCoordinates(i,j, nx, ny);
            nxm = recCoordinates(i-1,j, nx, ny);
            nxp = recCoordinates(i+1,j, nx, ny);
            nym = recCoordinates(i,j-1, nx, ny);
            nyp = recCoordinates(i,j+1, nx, ny);

            if (i == 1 && j == 1)
                if (bc_left == Inf)
                    rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                    ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;

                    G(n,n)   = -(rxp+ryp);
                    G(n,nxp) =  rxp;
                    G(n,nyp) =  ryp;
                else
                    G(n,n) = 1;
                    B(n) = bc_left;
                end
            elseif (i == 1 && j == ny)
                if (bc_left == Inf)
                    rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                    rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;

                    G(n,n)   = -(rxp+rym);
                    G(n,nxp) =  rxp;
                    G(n,nym) =  rym;
                else
                    G(n,n) = 1;
                    B(n) = bc_left;
                end
            elseif i == nx && j == 1 % Right side
                if (bc_right == Inf)
                    rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                    ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;
                    G(n,n)   = -(rxm+ryp);
                    G(n,nxm) =  rxm;
                    G(n,nyp) =  ryp;
                else
                    G(n,n) = 1;
                    B(n) = bc_right;
                end
            elseif i == nx && j == ny % Right side
                if (bc_right == Inf)
                    rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                    rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;
                    G(n,n)   = -(rxm+rym);
                    G(n,nxm) =  rxm;
                    G(n,nym) =  rym;
                else
                    G(n,n) = 1;
                    B(n) = bc_right;
                end
            elseif (i == 1) % Left Side
                if (bc_left == Inf)
                    rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                    rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;
                    ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;

                    G(n,n)   = -(rxp+rym+ryp);
                    G(n,nxp) =  rxp;
                    G(n,nym) =  rym;
                    G(n,nyp) =  ryp;
                else
                    G(n,n) = 1;
                    B(n) = bc_left;
                end
            elseif i == nx % Right side
                if (bc_right == Inf)
                    rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                    rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;
                    ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;
                    G(n,n)   = -(rxm+rym+ryp);
                    G(n,nxm) =  rxm;
                    G(n,nym) =  rym;
                    G(n,nyp) =  ryp;
                else
                    G(n,n) = 1;
                    B(n) = bc_right;
                end
            elseif j == 1 % Top Side
                if (bc_top == Inf)
                    rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                    rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                    ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;
                    G(n,n) = -(rxm+rxp+ryp);
                    G(n,nxm) =  rxm;
                    G(n,nxp) =  rxp;
                    G(n,nyp) =  ryp;
                else
                    G(n,n) = 1;
                    B(n) = bc_top;
                end
            elseif j == ny % Bottom Side
                if (bc_bottom == Inf)
                    rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                    rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                    rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;
                    G(n,n) = -(rxm+rxp+rym);
                    G(n,nxm) =  rxm;
                    G(n,nxp) =  rxp;
                    G(n,nym) =  rym;
                else
                    G(n,n) = 1;
                    B(n) = bc_bottom;
                end
            else % Bulk Area
                rxm = (crecCoordinates(i,j) + crecCoordinates(i-1,j))/2;
                rxp = (crecCoordinates(i,j) + crecCoordinates(i+1,j))/2;
                rym = (crecCoordinates(i,j) + crecCoordinates(i,j-1))/2;
                ryp = (crecCoordinates(i,j) + crecCoordinates(i,j+1))/2;
                
                G(n,n)   = -(rxm+rxp+rym+ryp);
                G(n,nxm) =  rxm;
                G(n,nxp) =  rxp;
                G(n,nym) =  rym;
                G(n,nyp) =  ryp;
            end
        end
    end
    
    V_temp = G\B';
    
    V = zeros(nx,ny,1);
    for i=1:nx
        for j=1:ny
            V(i,j) = V_temp(recCoordinates(i,j, nx, ny));
        end
    end
end
%%
% The following function maps the coordinates between the linear and
% rectangular plot
function [n] = recCoordinates(i,j, nx, ny)

    global C;
    n = j + (i - 1)*ny;
end
%%
% The following function calculates the total current for the needed
% parameters
function I = totalI(nx, ny, V1, sigma_conduct, sigma_insulate, Wb, Lb)


    cMap = sigma_conduct*ones(nx, ny);
    cMap(1:Wb,(1:Lb)+ny/2-Lb/2) = sigma_insulate;
    cMap((1:Wb)+nx-Wb,(1:Lb)+ny/2-Lb/2) = sigma_insulate;
    V = numericSolution(nx, ny, cMap, Inf, Inf, 0, V1);
    [Ex, Ey] = gradient(V);
    Ex = -Ex;
    Ey = -Ey;
    Jx = cMap.*Ex;
    I = (abs(sum(Jx(1,:))) + abs(sum(Jx(nx,:))))/2;
end
