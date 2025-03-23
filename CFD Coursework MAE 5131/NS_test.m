%attempt to solve Navier-Stokes on [0,1] x [0,1] with dirichlet boundary
%conditions

clear all; close all;

%number of spatial points
Nx = 50;
Ny = Nx;
%spatial discretization
x = linspace(0,1,Nx);
y = linspace(0,1,Ny);

%Reynold's number
Re = 200;

%grid spacing
dx = x(2) - x(1);
dy = y(2) - y(1);

%grid spacing
h = dx;

%meshes for plotting
[yPlot, xPlot] = meshgrid(y,x);

%time step size
dt = 0.000001;

%initial time
t0 = 0;

%final time
tend = 1;

onesVec = ones(Nx,1);

%second derivative matrix
D2 = 1/h^2*spdiags([onesVec -2*onesVec onesVec], [-1 0 1], Nx, Ny);

%identity matrix
I = eye(Nx,Ny);

%5 point laplacian operator
Lop = speye(Nx*Ny, Nx*Ny) + dt/(2*Re)*(kron(D2,I) + kron(I,D2));

onesVecP = ones(Nx-1,1);

D2_P = 1/h^2*spdiags([onesVecP -2*onesVecP onesVecP], [-1 0 1], Nx-2, Ny-2);
IP = eye(Nx-2, Ny-2);

Lop_p = speye((Nx-2)*(Ny-2), (Nx-2)*(Ny-2)) + dt/(2*Re)*(kron(D2_P,IP) + kron(IP,D2_P));





%number of time steps
Nt = tend/dt;

%initial condition x direction
u0 = zeros(Nx,Ny);
u0(:,end) = 1;

%initial condition y direction
v0 = zeros(Nx,Ny);

%second initial condition: for now just copy the first initial condition
%for sanity's sake
u1 = u0;
u1(:,end-1) = 0.1;
u1(end,end-1) = 0;
v1 = v0;
v1(:,end-1) = -0.01;

%time stepping loop
for k = 1:Nt

    unew = zeros(Nx, Ny);
    unew(:,end) = 1;

    vnew = zeros(Nx,Ny);

    ustar = zeros(Nx,Ny);
    vstar = zeros(Nx,Ny);

    pHold = zeros(Nx,Ny);
    
    %i component
    for i = 1:Nx
        %j component
        for j = 1:Ny

            %interior nodes
            if ((i > 1 && i < Nx) && (j > 1 && j < Ny))
                %solve for ustar (before it's hit with the laplacian)
                uDiffuse = 1/(2*Re*dx^2)*(u1(i+1,j) + u1(i-1,j) + u1(i,j+1) + u1(i,j-1) - 4*u1(i,j));
                uCN1 = 3/2*( 1/(4*dx)*( ( u1(i+1,j) + u1(i,j) )^2 - ( u1(i,j) + u1(i-1,j) )^2 ) + 1/(4*dy)*( (u1(i,j+1) + u1(i,j))*(v1(i,j+1) + v1(i,j)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i,j-1)) ) );
                uCN0 = -1/2*( 1/(4*dx)*( ( u0(i+1,j) + u0(i,j) )^2 - ( u0(i,j) + u0(i-1,j) )^2 ) + 1/(4*dy)*( (u0(i,j+1) + u0(i,j))*(v0(i,j+1) + v0(i,j)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i,j-1)) ) );
    
                ustar(i,j) = u1(i,j) - dt*uDiffuse - dt*(uCN1 + uCN0);
    
                %solve for vstar (before it's hit with the laplacian)
                vDiffuse = 1/(2*Re*dx^2)*(v1(i+1,j) + v1(i-1,j) + v1(i,j+1) + v1(i,j-1) - 4*v1(i,j));
                vCN1 = 3/2*( 1/(4*dx)*( (u1(i+1,j) + u1(i,j))*(v1(i+1,j) + v1(i,j)) - (u1(i,j) + u1(i-1,j))*(v1(i,j) + v1(i-1,j)) ) + 1/(4*dy)*( (v1(i,j+1) + v1(i,j))^2 - (v1(i,j) + v1(i,j-1))^2 ) );
                vCN0 = -1/2*( 1/(4*dx)*( (u0(i+1,j) + u0(i,j))*(v0(i+1,j) + v0(i,j)) - (u0(i,j) + u0(i-1,j))*(v0(i,j) + v0(i-1,j)) ) + 1/(4*dy)*( (v0(i,j+1) + v0(i,j))^2 - (v0(i,j) + v0(i,j-1))^2 ) );
                     
                vstar(i,j) = v1(i,j) - dt*vDiffuse - dt*(vCN1 + vCN0);

            %case we're at the top boundary but NOT at the corners
            %u(j+1)=1, v(j+1) = 0
            elseif (j == Ny && i > 1 && i < Nx)
                %solve for ustar (before it's hit with the laplacian)
                uDiffuse = 1/(2*Re*dx^2)*(u1(i+1,j) + u1(i-1,j) + 1 + u1(i,j-1) - 4*u1(i,j));
                uCN1 = 3/2*( 1/(4*dx)*( ( u1(i+1,j) + u1(i,j) )^2 - ( u1(i,j) + u1(i-1,j) )^2 ) + 1/(4*dy)*( (1 + u1(i,j))*(0 + v1(i,j)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i,j-1)) ) );
                uCN0 = -1/2*( 1/(4*dx)*( ( u0(i+1,j) + u0(i,j) )^2 - ( u0(i,j) + u0(i-1,j) )^2 ) + 1/(4*dy)*( (1 + u0(i,j))*(0 + v0(i,j)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i,j-1)) ) );
    
                ustar(i,j) = u1(i,j) - dt*uDiffuse - dt*(uCN1 + uCN0);
    
                %solve for vstar (before it's hit with the laplacian)
                vDiffuse = 1/(2*Re*dx^2)*(v1(i+1,j) + v1(i-1,j) + 0 + v1(i,j-1) - 4*v1(i,j));
                vCN1 = 3/2*( 1/(4*dx)*( (u1(i+1,j) + u1(i,j))*(v1(i+1,j) + v1(i,j)) - (u1(i,j) + u1(i-1,j))*(v1(i,j) + v1(i-1,j)) ) + 1/(4*dy)*( (0 + v1(i,j))^2 - (v1(i,j) + v1(i,j-1))^2 ) );
                vCN0 = -1/2*( 1/(4*dx)*( (u0(i+1,j) + u0(i,j))*(v0(i+1,j) + v0(i,j)) - (u0(i,j) + u0(i-1,j))*(v0(i,j) + v0(i-1,j)) ) + 1/(4*dy)*( (0 + v0(i,j))^2 - (v0(i,j) + v0(i,j-1))^2 ) );
                     
                vstar(i,j) = v1(i,j) - dt*vDiffuse - dt*(vCN1 + vCN0);

            %case we're at the bottom boundary but NOT at the corners
            %(j-1)=0
            elseif (j == 1 && i > 1 && i < Nx)
                %solve for ustar (before it's hit with the laplacian)
                uDiffuse = 1/(2*Re*dx^2)*(u1(i+1,j) + u1(i-1,j) + u1(i,j+1) + 0 - 4*u1(i,j));
                uCN1 = 3/2*( 1/(4*dx)*( ( u1(i+1,j) + u1(i,j) )^2 - ( u1(i,j) + u1(i-1,j) )^2 ) + 1/(4*dy)*( (u1(i,j+1) + u1(i,j))*(v1(i,j+1) + v1(i,j)) - (u1(i,j) + 0)*(v1(i,j) + 0) ) );
                uCN0 = -1/2*( 1/(4*dx)*( ( u0(i+1,j) + u0(i,j) )^2 - ( u0(i,j) + u0(i-1,j) )^2 ) + 1/(4*dy)*( (u0(i,j+1) + u0(i,j))*(v0(i,j+1) + v0(i,j)) - (u0(i,j) + 0)*(v0(i,j) + 0) ) );
    
                ustar(i,j) = u1(i,j) - dt*uDiffuse - dt*(uCN1 + uCN0);
    
                %solve for vstar (before it's hit with the laplacian)
                vDiffuse = 1/(2*Re*dx^2)*(v1(i+1,j) + v1(i-1,j) + v1(i,j+1) + 0 - 4*v1(i,j));
                vCN1 = 3/2*( 1/(4*dx)*( (u1(i+1,j) + u1(i,j))*(v1(i+1,j) + v1(i,j)) - (u1(i,j) + u1(i-1,j))*(v1(i,j) + v1(i-1,j)) ) + 1/(4*dy)*( (v1(i,j+1) + v1(i,j))^2 - (v1(i,j) + 0)^2 ) );
                vCN0 = -1/2*( 1/(4*dx)*( (u0(i+1,j) + u0(i,j))*(v0(i+1,j) + v0(i,j)) - (u0(i,j) + u0(i-1,j))*(v0(i,j) + v0(i-1,j)) ) + 1/(4*dy)*( (v0(i,j+1) + v0(i,j))^2 - (v0(i,j) + 0)^2 ) );
                     
                vstar(i,j) = v1(i,j) - dt*vDiffuse - dt*(vCN1 + vCN0);

            %Case we're at the left boundary but NOT at the corners (i-1)=0
            elseif (i == 1 && j > 1 && j < Ny)
                %solve for ustar (before it's hit with the laplacian)
                uDiffuse = 1/(2*Re*dx^2)*(u1(i+1,j) + 0 + u1(i,j+1) + u1(i,j-1) - 4*u1(i,j));
                uCN1 = 3/2*( 1/(4*dx)*( ( u1(i+1,j) + u1(i,j) )^2 - ( u1(i,j) + 0 )^2 ) + 1/(4*dy)*( (u1(i,j+1) + u1(i,j))*(v1(i,j+1) + v1(i,j)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i,j-1)) ) );
                uCN0 = -1/2*( 1/(4*dx)*( ( u0(i+1,j) + u0(i,j) )^2 - ( u0(i,j) + 0 )^2 ) + 1/(4*dy)*( (u0(i,j+1) + u0(i,j))*(v0(i,j+1) + v0(i,j)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i,j-1)) ) );
    
                ustar(i,j) = u1(i,j) - dt*uDiffuse - dt*(uCN1 + uCN0);
    
                %solve for vstar (before it's hit with the laplacian)
                vDiffuse = 1/(2*Re*dx^2)*(v1(i+1,j) + 0 + v1(i,j+1) + v1(i,j-1) - 4*v1(i,j));
                vCN1 = 3/2*( 1/(4*dx)*( (u1(i+1,j) + u1(i,j))*(v1(i+1,j) + v1(i,j)) - (u1(i,j) + 0)*(v1(i,j) + 0) ) + 1/(4*dy)*( (v1(i,j+1) + v1(i,j))^2 - (v1(i,j) + v1(i,j-1))^2 ) );
                vCN0 = -1/2*( 1/(4*dx)*( (u0(i+1,j) + u0(i,j))*(v0(i+1,j) + v0(i,j)) - (u0(i,j) + 0)*(v0(i,j) + 0) ) + 1/(4*dy)*( (v0(i,j+1) + v0(i,j))^2 - (v0(i,j) + v0(i,j-1))^2 ) );
                     
                vstar(i,j) = v1(i,j) - dt*vDiffuse - dt*(vCN1 + vCN0);

            %Case we're at the right boundary but NOT at the corners
            %(i+1)=0
            elseif (i == Nx && j > 1 && j < Ny)
                %solve for ustar (before it's hit with the laplacian)
                uDiffuse = 1/(2*Re*dx^2)*(0 + u1(i-1,j) + u1(i,j+1) + u1(i,j-1) - 4*u1(i,j));
                uCN1 = 3/2*( 1/(4*dx)*( ( 0 + u1(i,j) )^2 - ( u1(i,j) + u1(i-1,j) )^2 ) + 1/(4*dy)*( (u1(i,j+1) + u1(i,j))*(v1(i,j+1) + v1(i,j)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i,j-1)) ) );
                uCN0 = -1/2*( 1/(4*dx)*( ( 0 + u0(i,j) )^2 - ( u0(i,j) + u0(i-1,j) )^2 ) + 1/(4*dy)*( (u0(i,j+1) + u0(i,j))*(v0(i,j+1) + v0(i,j)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i,j-1)) ) );
    
                ustar(i,j) = u1(i,j) - dt*uDiffuse - dt*(uCN1 + uCN0);
    
                %solve for vstar (before it's hit with the laplacian)
                vDiffuse = 1/(2*Re*dx^2)*(0 + v1(i-1,j) + v1(i,j+1) + v1(i,j-1) - 4*v1(i,j));
                vCN1 = 3/2*( 1/(4*dx)*( (0 + u1(i,j))*(0 + v1(i,j)) - (u1(i,j) + u1(i-1,j))*(v1(i,j) + v1(i-1,j)) ) + 1/(4*dy)*( (v1(i,j+1) + v1(i,j))^2 - (v1(i,j) + v1(i,j-1))^2 ) );
                vCN0 = -1/2*( 1/(4*dx)*( (0 + u0(i,j))*(0 + v0(i,j)) - (u0(i,j) + u0(i-1,j))*(v0(i,j) + v0(i-1,j)) ) + 1/(4*dy)*( (v0(i,j+1) + v0(i,j))^2 - (v0(i,j) + v0(i,j-1))^2 ) );
                     
                vstar(i,j) = v1(i,j) - dt*vDiffuse - dt*(vCN1 + vCN0);

            %case bottom left corner (i-1) and (j-1) = 0
            elseif (i == 1 && j == 1)
               %solve for ustar (before it's hit with the laplacian)
                uDiffuse = 1/(2*Re*dx^2)*(u1(i+1,j) + 0 + u1(i,j+1) + 0 - 4*u1(i,j));
                uCN1 = 3/2*( 1/(4*dx)*( ( u1(i+1,j) + u1(i,j) )^2 - ( u1(i,j) + 0 )^2 ) + 1/(4*dy)*( (u1(i,j+1) + u1(i,j))*(v1(i,j+1) + v1(i,j)) - (u1(i,j) + 0)*(v1(i,j) + 0) ) );
                uCN0 = -1/2*( 1/(4*dx)*( ( u0(i+1,j) + u0(i,j) )^2 - ( u0(i,j) + 0 )^2 ) + 1/(4*dy)*( (u0(i,j+1) + u0(i,j))*(v0(i,j+1) + v0(i,j)) - (u0(i,j) + 0)*(v0(i,j) + 0) ) );
    
                ustar(i,j) = u1(i,j) - dt*uDiffuse - dt*(uCN1 + uCN0);
    
                %solve for vstar (before it's hit with the laplacian)
                vDiffuse = 1/(2*Re*dx^2)*(v1(i+1,j) + 0 + v1(i,j+1) + 0 - 4*v1(i,j));
                vCN1 = 3/2*( 1/(4*dx)*( (u1(i+1,j) + u1(i,j))*(v1(i+1,j) + v1(i,j)) - (u1(i,j) + 0)*(v1(i,j) + 0) ) + 1/(4*dy)*( (v1(i,j+1) + v1(i,j))^2 - (v1(i,j) + 0)^2 ) );
                vCN0 = -1/2*( 1/(4*dx)*( (u0(i+1,j) + u0(i,j))*(v0(i+1,j) + v0(i,j)) - (u0(i,j) + 0)*(v0(i,j) + 0) ) + 1/(4*dy)*( (v0(i,j+1) + v0(i,j))^2 - (v0(i,j) + 0)^2 ) );
                     
                vstar(i,j) = v1(i,j) - dt*vDiffuse - dt*(vCN1 + vCN0);

            %case top left corner u(j+1) = 1, v(j+1) = 0, (i-1) = 0
            elseif (i == 1 && j == Ny)
                %solve for ustar (before it's hit with the laplacian)
                uDiffuse = 1/(2*Re*dx^2)*(u1(i+1,j) +0 + 1 + u1(i,j-1) - 4*u1(i,j));
                uCN1 = 3/2*( 1/(4*dx)*( ( u1(i+1,j) + u1(i,j) )^2 - ( u1(i,j) + 0 )^2 ) + 1/(4*dy)*( (1 + u1(i,j))*(0 + v1(i,j)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i,j-1)) ) );
                uCN0 = -1/2*( 1/(4*dx)*( ( u0(i+1,j) + u0(i,j) )^2 - ( u0(i,j) + 0 )^2 ) + 1/(4*dy)*( (1 + u0(i,j))*(0 + v0(i,j)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i,j-1)) ) );
    
                ustar(i,j) = u1(i,j) - dt*uDiffuse - dt*(uCN1 + uCN0);
    
                %solve for vstar (before it's hit with the laplacian)
                vDiffuse = 1/(2*Re*dx^2)*(v1(i+1,j) + 0 + 0 + v1(i,j-1) - 4*v1(i,j));
                vCN1 = 3/2*( 1/(4*dx)*( (u1(i+1,j) + u1(i,j))*(v1(i+1,j) + v1(i,j)) - (u1(i,j) + 0)*(v1(i,j) + 0) ) + 1/(4*dy)*( (0 + v1(i,j))^2 - (v1(i,j) + v1(i,j-1))^2 ) );
                vCN0 = -1/2*( 1/(4*dx)*( (u0(i+1,j) + u0(i,j))*(v0(i+1,j) + v0(i,j)) - (u0(i,j) + 0)*(v0(i,j) + 0) ) + 1/(4*dy)*( (0 + v0(i,j))^2 - (v0(i,j) + v0(i,j-1))^2 ) );
                     
                vstar(i,j) = v1(i,j) - dt*vDiffuse - dt*(vCN1 + vCN0);

            %case top right corner u(j+1) = 1 v(j+1) = 0 (i+1) = 0
            elseif (i == Nx && j == Ny)
                %solve for ustar (before it's hit with the laplacian)
                uDiffuse = 1/(2*Re*dx^2)*(0 + u1(i-1,j) + 1 + u1(i,j-1) - 4*u1(i,j));
                uCN1 = 3/2*( 1/(4*dx)*( ( 0 + u1(i,j) )^2 - ( u1(i,j) + u1(i-1,j) )^2 ) + 1/(4*dy)*( (1 + u1(i,j))*(0 + v1(i,j)) - (u1(i,j) + u1(i,j-1))*(v1(i,j) + v1(i,j-1)) ) );
                uCN0 = -1/2*( 1/(4*dx)*( ( 0 + u0(i,j) )^2 - ( u0(i,j) + u0(i-1,j) )^2 ) + 1/(4*dy)*( (1 + u0(i,j))*(0 + v0(i,j)) - (u0(i,j) + u0(i,j-1))*(v0(i,j) + v0(i,j-1)) ) );
    
                ustar(i,j) = u1(i,j) - dt*uDiffuse - dt*(uCN1 + uCN0);
    
                %solve for vstar (before it's hit with the laplacian)
                vDiffuse = 1/(2*Re*dx^2)*(0 + v1(i-1,j) + 0 + v1(i,j-1) - 4*v1(i,j));
                vCN1 = 3/2*( 1/(4*dx)*( (0 + u1(i,j))*(0 + v1(i,j)) - (u1(i,j) + u1(i-1,j))*(v1(i,j) + v1(i-1,j)) ) + 1/(4*dy)*( (0 + v1(i,j))^2 - (v1(i,j) + v1(i,j-1))^2 ) );
                vCN0 = -1/2*( 1/(4*dx)*( (0 + u0(i,j))*(0 + v0(i,j)) - (u0(i,j) + u0(i-1,j))*(v0(i,j) + v0(i-1,j)) ) + 1/(4*dy)*( (0 + v0(i,j))^2 - (v0(i,j) + v0(i,j-1))^2 ) );
                     
                vstar(i,j) = v1(i,j) - dt*vDiffuse - dt*(vCN1 + vCN0);

            %case bottom right corner (j-1) = (i+1) = 0
            elseif ( i == Nx && j == 1)
                %solve for ustar (before it's hit with the laplacian)
                uDiffuse = 1/(2*Re*dx^2)*(0 + u1(i-1,j) + u1(i,j+1) + 0 - 4*u1(i,j));
                uCN1 = 3/2*( 1/(4*dx)*( ( 0 + u1(i,j) )^2 - ( u1(i,j) + u1(i-1,j) )^2 ) + 1/(4*dy)*( (u1(i,j+1) + u1(i,j))*(v1(i,j+1) + v1(i,j)) - (u1(i,j) + 0)*(v1(i,j) + 0) ) );
                uCN0 = -1/2*( 1/(4*dx)*( ( 0 + u0(i,j) )^2 - ( u0(i,j) + u0(i-1,j) )^2 ) + 1/(4*dy)*( (u0(i,j+1) + u0(i,j))*(v0(i,j+1) + v0(i,j)) - (u0(i,j) + 0)*(v0(i,j) + 0) ) );
    
                ustar(i,j) = u1(i,j) - dt*uDiffuse - dt*(uCN1 + uCN0);
    
                %solve for vstar (before it's hit with the laplacian)
                vDiffuse = 1/(2*Re*dx^2)*(0 + v1(i-1,j) + v1(i,j+1) + 0 - 4*v1(i,j));
                vCN1 = 3/2*( 1/(4*dx)*( (0 + u1(i,j))*(0 + v1(i,j)) - (u1(i,j) + u1(i-1,j))*(v1(i,j) + v1(i-1,j)) ) + 1/(4*dy)*( (v1(i,j+1) + v1(i,j))^2 - (v1(i,j) + 0)^2 ) );
                vCN0 = -1/2*( 1/(4*dx)*( (0 + u0(i,j))*(0 + v0(i,j)) - (u0(i,j) + u0(i-1,j))*(v0(i,j) + v0(i-1,j)) ) + 1/(4*dy)*( (v0(i,j+1) + v0(i,j))^2 - (v0(i,j) + 0)^2 ) );
                     
                vstar(i,j) = v1(i,j) - dt*vDiffuse - dt*(vCN1 + vCN0);
            end


           




        end %end for loop j component
    end %end for loop i component

    ustar = reshape(ustar, Nx*Ny,1);
    vstar = reshape(vstar, Nx*Ny,1);

    ustar = Lop\ustar;
    vstar = Lop\vstar;

    ustar = reshape(ustar, Nx,Ny);
    vstar = reshape(vstar, Nx,Ny);


    for i = 2:Nx-1
        for j = 2:Ny-1
            pHold(i,j) = 1/(dt)*( 1/dx*(ustar(i+1,j) - ustar(i-1,j)) + 1/dy*(vstar(i,j+1) - vstar(i,j-1)) );
        end
    end

    pHoldL = reshape(pHold(2:end-1,2:end-1), (Nx-2)*(Ny-2),1);
    pHoldL = Lop_p\pHoldL;
    pHold(2:end-1,2:end-1) = reshape(pHoldL, Nx-2, Ny-2);

    for i=2:Nx-1
        for j=2:Ny-1
            unew(i,j) = ustar(i,j) - dt/(2*dx)*( pHold(i+1,j) - pHold(i-1,j) );
            vnew(i,j) = vstar(i,j) - dt/(2*dy)*( pHold(i,j+1) - pHold(i,j-1) );
        end
    end

    u0 = u1;
    v0 = v1;

    u1 = unew;
    v1 = vnew;
    
 
    quiver(xPlot,yPlot,u1,v1,5)
    drawnow
    pause(0.1)
    
end

quiver(xPlot, yPlot, u1,v1)