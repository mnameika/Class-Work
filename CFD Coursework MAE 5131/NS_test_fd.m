%finite difference approach to NS
clear all; close all;

%number of point in x and y, respectively
Nx = 50;
Ny = Nx;

%x and y vectors
x = linspace(0,1,Nx);
y = linspace(0,1,Nx);

[yPlot, xPlot] = meshgrid(y,x);

%reynold's number
Re = 200;

%spatial step size
h = x(2) - x(1);

%final time
tEnd = 20;

%time step size
dt = 0.001;

%number of time steps
Nt = tEnd/dt;


%initial condition of u component of velocity
u0 = zeros(Nx,Ny);
u0(:,end) = 1;

%initial condition of v component of velocity
v0 = zeros(Nx,Ny);


onesVec = ones(Nx,1);

%second order finite difference matrix for second derivative
D2 = 1/h^2*spdiags([onesVec -2*onesVec onesVec], [-1 0 1], Nx, Ny);

%identity matrix
I = speye(Nx,Ny);

%5 point laplacian
Lop = kron(I,D2) + kron(D2,I);


animFile = 'C:\Users\Michael\Desktop\Computational Fluid Dynamics\NS_animation.gif';

[yFine, xFine] = meshgrid(0:0.01:1, 0:0.01:1);

%time stepping
for k = 1:Nt
    ustar = zeros(Nx,Ny);
    vstar = zeros(Nx,Ny);

    unew = zeros(Nx,Ny);
    vnew = zeros(Nx,Ny);

    pHold = zeros(Nx,Ny);

    for i = 1:Nx
        for j = 1:Ny
            
            %interior points
            if (i > 1 && i < Nx && j > 1 && j < Ny)
                ustar(i,j) = u0(i,j) - dt*( 1/(2*h)*( u0(i+1,j)^2 - u0(i-1,j)^2  + u0(i,j+1)*v0(i,j+1) - u0(i,j-1)*v0(i,j-1)) - 1/(Re*h^2)*( u0(i+1,j) + u0(i-1,j) + u0(i,j+1) + u0(i,j-1) - 4*u0(i,j) ) );
    
                vstar(i,j) = v0(i,j) - dt*( 1/(2*h)*( v0(i,j+1)^2 - v0(i,j-1)^2 + u0(i+1,j)*v0(i+1,j) - u0(i-1,j)*v0(i-1,j)) - 1/(Re*h^2)*( v0(i+1,j) + v0(i-1,j) + v0(i,j-1) + v0(i,j+1) - 4*v0(i,j) ) );

            %left boundary, no corner. u(i-1) = 0, v(i-1) = 0
            elseif (i == 1 && j > 1 && j < Ny)
                ustar(i,j) = u0(i,j) - dt*( 1/(2*h)*( u0(i+1,j)^2 - 0^2  + u0(i,j+1)*v0(i,j+1) - u0(i,j-1)*v0(i,j-1)) - 1/(Re*h^2)*( u0(i+1,j) + 0 + u0(i,j+1) + u0(i,j-1) - 4*u0(i,j) ) );
    
                vstar(i,j) = v0(i,j) - dt*( 1/(2*h)*( v0(i,j+1)^2 - v0(i,j-1)^2 + u0(i+1,j)*v0(i+1,j) - 0*0) - 1/(Re*h^2)*( v0(i+1,j) + 0 + v0(i,j-1) + v0(i,j+1) - 4*v0(i,j) ) );

            %right boundary, no corner. u(i+1) = 0, v(i+1) = 0
            elseif (i == Nx && j > 1 && j < Ny)
                ustar(i,j) = u0(i,j) - dt*( 1/(2*h)*( 0^2 - u0(i-1,j)^2  + u0(i,j+1)*v0(i,j+1) - u0(i,j-1)*v0(i,j-1)) - 1/(Re*h^2)*( 0 + u0(i-1,j) + u0(i,j+1) + u0(i,j-1) - 4*u0(i,j) ) );
    
                vstar(i,j) = v0(i,j) - dt*( 1/(2*h)*( v0(i,j+1)^2 - v0(i,j-1)^2 + 0*0 - u0(i-1,j)*v0(i-1,j)) - 1/(Re*h^2)*( 0 + v0(i-1,j) + v0(i,j-1) + v0(i,j+1) - 4*v0(i,j) ) );

            %bottom boundary, no corner. u(j-1) = 0, v(j-1) = 0
            elseif (i > 1 && i < Nx && j == 1)
                ustar(i,j) = u0(i,j) - dt*( 1/(2*h)*( u0(i+1,j)^2 - u0(i-1,j)^2  + u0(i,j+1)*v0(i,j+1) - 0*0) - 1/(Re*h^2)*( u0(i+1,j) + u0(i-1,j) + u0(i,j+1) + 0 - 4*u0(i,j) ) );
    
                vstar(i,j) = v0(i,j) - dt*( 1/(2*h)*( v0(i,j+1)^2 - 0^2 + u0(i+1,j)*v0(i+1,j) - u0(i-1,j)*v0(i-1,j)) - 1/(Re*h^2)*( v0(i+1,j) + v0(i-1,j) + 0 + v0(i,j+1) - 4*v0(i,j) ) );

            %top boundary, no corner. u(j+1) = 1, v(j+1) = 0
            elseif (i > 1 && i < Nx && j == Ny)
                ustar(i,j) = u0(i,j) - dt*( 1/(2*h)*( u0(i+1,j)^2 - u0(i-1,j)^2  + 1*0 - u0(i,j-1)*v0(i,j-1)) - 1/(Re*h^2)*( u0(i+1,j) + u0(i-1,j) + 1 + u0(i,j-1) - 4*u0(i,j) ) );
    
                vstar(i,j) = v0(i,j) - dt*( 1/(2*h)*( 0^2 - v0(i,j-1)^2 + u0(i+1,j)*v0(i+1,j) - u0(i-1,j)*v0(i-1,j)) - 1/(Re*h^2)*( v0(i+1,j) + v0(i-1,j) + v0(i,j-1) + 0 - 4*v0(i,j) ) );

            %bottom left corner. u(j-1) = u(i-1) = 0, v(j-1) = v(i-1) = 0
            elseif (i == 1 && j == 1)
                ustar(i,j) = u0(i,j) - dt*( 1/(2*h)*( u0(i+1,j)^2 - 0^2 + u0(i,j+1)*v0(i,j+1) - 0*0) - 1/(Re*h^2)*( u0(i+1,j) + 0 + u0(i,j+1) + 0 - 4*u0(i,j) ) );
    
                vstar(i,j) = v0(i,j) - dt*( 1/(2*h)*( v0(i,j+1)^2 - 0^2 + u0(i+1,j)*v0(i+1,j) - 0*0) - 1/(Re*h^2)*( v0(i+1,j) + 0 + 0 + v0(i,j+1) - 4*v0(i,j) ) );

            %top left corner. u(j+1) = 1, u(i-1) = 0, v(j+1) = v(i-1) = 0
            elseif (i == 1 && j == Ny)
                ustar(i,j) = u0(i,j) - dt*( 1/(2*h)*( u0(i+1,j)^2 - 0^2  + 1*0 - u0(i,j-1)*v0(i,j-1)) - 1/(Re*h^2)*( u0(i+1,j) + 0 + 1 + u0(i,j-1) - 4*u0(i,j) ) );
    
                vstar(i,j) = v0(i,j) - dt*( 1/(2*h)*( 0^2 - v0(i,j-1)^2 + u0(i+1,j)*v0(i+1,j) - 0*0) - 1/(Re*h^2)*( v0(i+1,j) + 0 + v0(i,j-1) + 0 - 4*v0(i,j) ) );

            %top right corner. u(j+1) = 1, u(i+1) = 0, v(j+1) = v(i+1) = 0
            elseif (i == Nx && j == Ny)
                ustar(i,j) = u0(i,j) - dt*( 1/(2*h)*( 0^2 - u0(i-1,j)^2  + 1*0 - u0(i,j-1)*v0(i,j-1)) - 1/(Re*h^2)*( 0 + u0(i-1,j) + 1 + u0(i,j-1) - 4*u0(i,j) ) );
    
                vstar(i,j) = v0(i,j) - dt*( 1/(2*h)*( 0^2 - v0(i,j-1)^2 + 0*0 - u0(i-1,j)*v0(i-1,j)) - 1/(Re*h^2)*( 0 + v0(i-1,j) + v0(i,j-1) + 0 - 4*v0(i,j) ) );

            %bottom right corner. u(j-1) = u(i+1) = 0, v(j-1) = v(i+1) = 0
            elseif (i == Nx && j == 1)
                ustar(i,j) = u0(i,j) - dt*( 1/(2*h)*( 0^2 - u0(i-1,j)^2  + u0(i,j+1)*v0(i,j+1) - 0*0) - 1/(Re*h^2)*( 0 + u0(i-1,j) + u0(i,j+1) + 0 - 4*u0(i,j) ) );
    
                vstar(i,j) = v0(i,j) - dt*( 1/(2*h)*( v0(i,j+1)^2 - 0^2 + 0*0 - u0(i-1,j)*v0(i-1,j)) - 1/(Re*h^2)*( 0 + v0(i-1,j) + 0 + v0(i,j+1) - 4*v0(i,j) ) );
            end %end if

        end %end for loop for j index
    end %end for loop for i index
    

    %loop for solving pressure equation
    for i = 2:Nx-1

        for j = 2:Ny-1
            pHold(i,j) = 1/(2*dt*h)*( ustar(i+1,j) - ustar(i-1,j) + vstar(i,j+1) - vstar(i,j-1) );
            
        end %end for loop for j index
    end %end for loop for i index

    pHold = reshape(pHold, (Nx)*(Ny),1);

    pHold = Lop\pHold;

    pHold = reshape(pHold, Nx, Ny);

    
    for i=2:Nx-1
        for j = 2:Ny-1
            unew(i,j) = ustar(i,j) - dt/(2*h)*( pHold(i+1,j) - pHold(i-1,j) );
            vnew(i,j) = vstar(i,j) - dt/(2*h)*( pHold(i,j+1) - pHold(i,j-1) );
        end
    end

    u0 = unew;
    u0(:,end) = 1;
    v0 = vnew;

    % if (rem(k-2, 1000) == 0)
    %     pFine = interp2(yPlot,xPlot, pHold, yFine, xFine, 'spline');
    % 
    %     contourf(xFine, yFine, pFine, -15:0.01:15, 'LineStyle', 'none')
    %     hold on
    %     quiver(xPlot,yPlot,u0,v0,'-k', 'linewidth', 0.9)
    %     axis([0 1 0 1])
    %     colorbar
    %     colormap('Turbo')
    %     % hold off
    %     % drawnow
    % 
    %     exportgraphics(gca, animFile, 'Append', true)
    % end
    prog = k/Nt
end %end for loop for time stepping

% [yFine, xFine] = meshgrid(0:1/32:1, 0:1/32:1);

pFine = interp2(yPlot,xPlot, pHold, yFine, xFine, 'cubic');



contourf(xFine, yFine, pFine, -15:0.01:15, 'LineStyle', 'none')
hold on
quiver(xPlot,yPlot,u0,v0,'-k', 'linewidth', 0.9)
axis([0 1 0 1])
colorbar
colormap('Turbo')


hold off