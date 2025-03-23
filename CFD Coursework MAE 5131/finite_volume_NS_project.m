%finite volume NS implementation
clear all; close all;

%number of spatial points
Nx = 100;
Ny = Nx;

%spatial grid spacing
h = 1/(Nx-2);
% h = 1/(Nx - 1);

%initial and final times
t0 = 0;
tEnd = 2;

%time step size
dt = 0.001/8;

%number of time steps
Nt = (tEnd - t0)/dt;

%Reynold's number
Re = 200;

%spatial coordinates
x = [0 h/2:h:1-h/2 1];
% x = linspace(0,1,Nx);
y = x;

%meshgrid plots for x and y for plotting purposes
[Xplot, Yplot] = meshgrid(x,y);

%initial conditions for velocity field
u0 = zeros(Nx,Ny);
u0(:,end) = 1;

v0 = zeros(Nx,Ny);


%Laplacian operator for solving pressure poisson equation
onesVec = ones(Nx,1);

%second derivative matrix
D2 = 1/h^2*spdiags([onesVec -2*onesVec onesVec], [-1 0 1], Nx, Ny);

I = speye(Nx,Nx);

%5 point laplacian matrix
Lop = kron(D2,I) + kron(I,D2);

%nonlinear convection terms
N_x = @(ue,uw,un,us,ve,vw,vn,vs) 1/h*( ue^2 - uw^2 + un*vn - us*vs );
N_y = @(ue,uw,un,us,ve,vw,vn,vs) 1/h*( ue*ve - uw*vw + vn^2 - vs^2 );

%diffusion terms
%use ghost cells for boundary conditions
L_x = @(uP,uE,uW,uN,uS) 1/(Re*h^2)*( uN + uW + uS + uE - 4*uP );
L_y = @(vP,vE,vW,vN,vS) 1/(Re*h^2)*( vN + vW + vS + vE - 4*vP );

%RHS of pressure equation
R = @(ue,uw,un,us,ve,vw,vn,vs) 1/h*( ue - uw + vn - vs );


%time stepping
for k = 1:Nt
    %initialize the starred variablesx
    unew = zeros(Nx,Ny);
    vnew = zeros(Nx,Ny);

    ustar = zeros(Nx,Ny);
    vstar = zeros(Nx,Ny);

    pHold = zeros(Nx,Ny);
    
    for i = 1:Nx
        for j = 1:Ny
            
            %interior nodes
            if (i > 1 && i < Nx && j > 1 && j < Ny)

                %interpolate nodes 
                ue = (u0(i+1,j) + u0(i,j))/2;
                uw = (u0(i-1,j) + u0(i,j))/2;
                un = (u0(i,j+1) + u0(i,j))/2;
                us = (u0(i,j-1) + u0(i,j))/2;

                ve = (v0(i+1,j) + v0(i,j))/2;
                vw = (v0(i-1,j) + v0(i,j))/2;
                vn = (v0(i,j+1) + v0(i,j))/2;
                vs = (v0(i,j-1) + v0(i,j))/2;

                
                N1 = N_x(ue, uw, un, us, ve, vw, vn, vs);
                N2 = N_y(ue, uw, un, us, ve, vw, vn, vs);

                L1 = L_x(u0(i,j), u0(i+1,j), u0(i-1,j), u0(i,j+1), u0(i,j-1));
                L2 = L_y(v0(i,j), v0(i+1,j), v0(i-1,j), v0(i,j+1), v0(i,j-1));

                R1 = R(ue, uw, un, us, ve, vw, vn, vs);

                ustar(i,j) = u0(i,j) - dt*N1 + dt*L1;
                vstar(i,j) = v0(i,j) - dt*N2 + dt*L2;
                
            %Bottom boundary NO CORNERS
            elseif (i > 1 && i < Nx && j == 1)
                ue = (u0(i+1,j) + u0(i,j))/2;
                uw = (u0(i-1,j) + u0(i,j))/2;
                un = (u0(i,j+1) + u0(i,j))/2;
                us = 0;

                ve = (v0(i+1,j) + v0(i,j))/2;
                vw = (v0(i-1,j) + v0(i,j))/2;
                vn = (v0(i,j+1) + v0(i,j))/2;
                vs = 0;

                
                N1 = N_x(ue, uw, un, us, ve, vw, vn, vs);
                N2 = N_y(ue, uw, un, us, ve, vw, vn, vs);

                L1 = L_x(u0(i,j), u0(i+1,j), u0(i-1,j), u0(i,j+1), 0);
                L2 = L_y(v0(i,j), v0(i+1,j), v0(i-1,j), v0(i,j+1), 0);

                R1 = R(ue, uw, un, us, ve, vw, vn, vs);

                ustar(i,j) = u0(i,j) - dt*N1 + dt*L1;
                vstar(i,j) = v0(i,j) - dt*N2 + dt*L2;

            %top boundary NO CORNERS
            elseif (i > 1 && i < Nx && j == Ny)
                ue = (u0(i+1,j) + u0(i,j))/2;
                uw = (u0(i-1,j) + u0(i,j))/2;
                un = 1;
                us = (u0(i,j-1) + u0(i,j))/2;

                ve = (v0(i+1,j) + v0(i,j))/2;
                vw = (v0(i-1,j) + v0(i,j))/2;
                vn = 0;
                vs = (v0(i,j-1) + v0(i,j))/2;

                
                N1 = N_x(ue, uw, un, us, ve, vw, vn, vs);
                N2 = N_y(ue, uw, un, us, ve, vw, vn, vs);

                L1 = L_x(u0(i,j), u0(i+1,j), u0(i-1,j), 1, u0(i,j-1));
                L2 = L_y(v0(i,j), v0(i+1,j), v0(i-1,j), 0, v0(i,j-1));

                R1 = R(ue, uw, un, us, ve, vw, vn, vs);

                ustar(i,j) = u0(i,j) - dt*N1 + dt*L1;
                vstar(i,j) = v0(i,j) - dt*N2 + dt*L2;
            
            %left boundary NO CORNERS
            elseif (j > 1 && j < Ny && i == 1)
                ue = (u0(i+1,j) + u0(i,j))/2;
                uw = 0;
                un = (u0(i,j+1) + u0(i,j))/2;
                us = (u0(i,j-1) + u0(i,j))/2;

                ve = (v0(i+1,j) + v0(i,j))/2;
                vw = 0;
                vn = (v0(i,j+1) + v0(i,j))/2;
                vs = (v0(i,j-1) + v0(i,j))/2;

                
                N1 = N_x(ue, uw, un, us, ve, vw, vn, vs);
                N2 = N_y(ue, uw, un, us, ve, vw, vn, vs);

                L1 = L_x(u0(i,j), u0(i+1,j), 0, u0(i,j+1), u0(i,j-1));
                L2 = L_y(v0(i,j), v0(i+1,j), 0, v0(i,j+1), v0(i,j-1));

                R1 = R(ue, uw, un, us, ve, vw, vn, vs);

                ustar(i,j) = u0(i,j) - dt*N1 + dt*L1;
                vstar(i,j) = v0(i,j) - dt*N2 + dt*L2;


            %right boundary NO CORNERS
            elseif (j > 1 && j < Ny && i == Nx)
                ue = 0;
                uw = (u0(i-1,j) + u0(i,j))/2;
                un = (u0(i,j+1) + u0(i,j))/2;
                us = (u0(i,j-1) + u0(i,j))/2;

                ve = 0;
                vw = (v0(i-1,j) + v0(i,j))/2;
                vn = (v0(i,j+1) + v0(i,j))/2;
                vs = (v0(i,j-1) + v0(i,j))/2;

                
                N1 = N_x(ue, uw, un, us, ve, vw, vn, vs);
                N2 = N_y(ue, uw, un, us, ve, vw, vn, vs);

                L1 = L_x(u0(i,j), 0, u0(i-1,j), u0(i,j+1), u0(i,j-1));
                L2 = L_y(v0(i,j), 0, v0(i-1,j), v0(i,j+1), v0(i,j-1));

                R1 = R(ue, uw, un, us, ve, vw, vn, vs);

                ustar(i,j) = u0(i,j) - dt*N1 + dt*L1;
                vstar(i,j) = v0(i,j) - dt*N2 + dt*L2;

            %bottom left corner
            elseif (i == 1 && j == 1)
                ue = (u0(i+1,j) + u0(i,j))/2;
                uw = 0;
                un = (u0(i,j+1) + u0(i,j))/2;
                us = 0;

                ve = (v0(i+1,j) + v0(i,j))/2;
                vw = 0;
                vn = (v0(i,j+1) + v0(i,j))/2;
                vs = 0;

                
                N1 = N_x(ue, uw, un, us, ve, vw, vn, vs);
                N2 = N_y(ue, uw, un, us, ve, vw, vn, vs);

                L1 = L_x(u0(i,j), u0(i+1,j), 0, u0(i,j+1), 0);
                L2 = L_y(v0(i,j), v0(i+1,j), 0, v0(i,j+1), 0);

                R1 = R(ue, uw, un, us, ve, vw, vn, vs);

                ustar(i,j) = u0(i,j) - dt*N1 + dt*L1;
                vstar(i,j) = v0(i,j) - dt*N2 + dt*L2;

            %top left corner
            elseif (i == 1 && j == Ny)
                ue = (u0(i+1,j) + u0(i,j))/2;
                uw = 0;
                un = 1;
                us = (u0(i,j-1) + u0(i,j))/2;

                ve = (v0(i+1,j) + v0(i,j))/2;
                vw = 0;
                vn = 0;
                vs = (v0(i,j-1) + v0(i,j))/2;

                
                N1 = N_x(ue, uw, un, us, ve, vw, vn, vs);
                N2 = N_y(ue, uw, un, us, ve, vw, vn, vs);

                L1 = L_x(u0(i,j), u0(i+1,j), 0, 1, u0(i,j-1));
                L2 = L_y(v0(i,j), v0(i+1,j), 0, 0, v0(i,j-1));

                R1 = R(ue, uw, un, us, ve, vw, vn, vs);

                ustar(i,j) = u0(i,j) - dt*N1 + dt*L1;
                vstar(i,j) = v0(i,j) - dt*N2 + dt*L2;

            %bottom right corner
            elseif (i == Nx && j == 1)
                ue = 0;
                uw = (u0(i-1,j) + u0(i,j))/2;
                un = (u0(i,j+1) + u0(i,j))/2;
                us = 0;

                ve = 0;
                vw = (v0(i-1,j) + v0(i,j))/2;
                vn = (v0(i,j+1) + v0(i,j))/2;
                vs = 0;

                
                N1 = N_x(ue, uw, un, us, ve, vw, vn, vs);
                N2 = N_y(ue, uw, un, us, ve, vw, vn, vs);

                L1 = L_x(u0(i,j), 0, u0(i-1,j), u0(i,j+1), 0);
                L2 = L_y(v0(i,j), 0, v0(i-1,j), v0(i,j+1), 0);

                R1 = R(ue, uw, un, us, ve, vw, vn, vs);

                ustar(i,j) = u0(i,j) - dt*N1 + dt*L1;
                vstar(i,j) = v0(i,j) - dt*N2 + dt*L2;

            %top right corner
            elseif (i == Nx && j == Ny)
                ue = 0;
                uw = (u0(i-1,j) + u0(i,j))/2;
                un = 1;
                us = (u0(i,j-1) + u0(i,j))/2;

                ve = 0;
                vw = (v0(i-1,j) + v0(i,j))/2;
                vn = 0;
                vs = (v0(i,j-1) + v0(i,j))/2;

                
                N1 = N_x(ue, uw, un, us, ve, vw, vn, vs);
                N2 = N_y(ue, uw, un, us, ve, vw, vn, vs);

                L1 = L_x(u0(i,j), 0, u0(i-1,j), 1, u0(i,j-1));
                L2 = L_y(v0(i,j), 0, v0(i-1,j), 0, v0(i,j-1));

                R1 = R(ue, uw, un, us, ve, vw, vn, vs);

                ustar(i,j) = u0(i,j) - dt*N1 + dt*L1;
                vstar(i,j) = v0(i,j) - dt*N2 + dt*L2;
            end



             %solve for the pressure
            pHold(i,j) = 1/dt*R(ue,uw,un,us,ve,vw,vn,vs);
           
        end %end y indexing for loop
    end %end x indexing for loop


    pHold = reshape(pHold, Nx*Ny, 1);

    pHold = Lop\pHold;

    pHold = reshape(pHold, Nx,Ny);

    for i = 2:Nx-1
        for j = 2:Ny-1
            Px = 1/(2*h)*( pHold(i+1,j) - pHold(i-1,j) );
            Py = 1/(2*h)*( pHold(i,j+1) - pHold(i,j-1) );

            %update velocity field
            unew(i,j) = ustar(i,j) - dt*Px;
            vnew(i,j) = vstar(i,j) - dt*Py;

        end
    end

    %update velocity values
    u0 = unew;
    u0(:,end) = 1;

    v0 = vnew;

    %display progress of number of time steps
    prog = k/Nt

end %end time stepping loop


%fine grid
NFine = Nx^2/4;
xFine = linspace(0,1,NFine);
yFine = xFine;

XFine = Xplot;
YFine = Yplot;
pFine = pHold;
% [XFine, YFine] = meshgrid(xFine, yFine);
% pFine = interp2(Xplot,Yplot, pHold, XFine, YFine, 'spline');


contourf(XFine, YFine, pFine', -1:0.01:1, 'LineStyle', 'none')
hold on
qspace = 10;

quiver(Xplot(1:qspace:end,1:qspace:end), Yplot(1:qspace:end,1:qspace:end), u0(1:qspace:end,1:qspace:end)',v0(1:qspace:end,1:qspace:end)', 1, 'k-', 'linewidth', 0.5)
cb = colorbar();
colormap('Turbo')
% hold off

axis([0 1 0 1])
xlabel('$x$', 'fontsize', 45, 'interpreter', 'latex')
ylabel('$y$', 'fontsize', 45, 'interpreter', 'latex')
ylabel(cb, 'Pressure $P$', 'fontsize', 45, 'interpreter', 'latex')

strmln1 = stream2(Xplot,Yplot, u0',v0', 0.7755,0.7755);
strmln2 = stream2(Xplot,Yplot, u0', v0', 0.995, 0.005);
strmln3 = stream2(Xplot,Yplot, u0', v0', 0.7602, 0.413);
strmln4 = stream2(Xplot,Yplot, u0', v0', 0.6786, 0.7092);

strmPts = 2000;

strmDbl1 = strmln1{1};
strmDbl2 = strmln2{1};
strmDbl3 = strmln3{1};
strmDbl4 = strmln4{1};

plot(strmDbl1(1:strmPts,1), strmDbl1(1:strmPts,2), 'k-', 'linewidth', 0.8)
plot(strmDbl2(:,1), strmDbl2(:,2), 'k-', 'linewidth', 0.8)
% plot(strmDbl3(:,1), strmDbl3(:,2), 'k-', 'linewidth', 0.8)
% plot(strmDbl4(:,1), strmDbl4(:,2), 'k-', 'linewidth', 0.8)

hold off
title("Velocity field for lid driven cavity ($N = $ " + Nx + ", $\Delta t = $ " + dt + ")", 'fontsize', 25, 'interpreter', 'latex')