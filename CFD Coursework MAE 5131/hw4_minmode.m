%FEM upwind
clear all; close all;

%spatial grid spacing
dx = 0.01;

%spatial grid endpoints
x0 = 0;
xend = 5;

%final time
tend = 2;

%x discretization
x = x0+dx/2:dx:xend-dx/2;

%number of spatial points
Nx = length(x);


%initially initialize initial condition 
u = zeros(length(x),1);

%temporal grid spacing
dt = dx;


method =1;

%Limiter functions
switch method
    %Upwind scheme
    case 1
        psi = @(r) 0*r;
        methodStr = "Upwind";

    %Minmode scheme
    case 2
        psi = @(r) max(0,min(1,r));
        methodStr = "Minmode";

    %Superbee scheme
    case 3
        psi = @(r) max([0, min(2*r,1), min(r,2)]);
        methodStr = "Superbee";

    %Van Leer scheme
    case 4
        psi = @(r) (r + abs(r))./(1 + r);
        methodStr = "Van Leer";

    %QUICK scheme
    case 5
        psi = @(r) 1/4*(3 + r);
        methodStr = "QUICK";

    %MUSCL scheme
    case 6
        psi = @(r) max([0, min([2*r, (r+1)/2,2])]);
        methodStr = "MUSCL";
end


%wave speed
c = 1;

%exact solution (for Gaussian - need to add square wave)
% uEx = @(t) exp(-100*(x - 3/10 - c*t).^2);
uEx = zeros(length(x),1);

%build initial condition
for i = 1:length(x)
    if (0 <= x(i) && x(i) <= 0.6)
       u(i) = exp(-100*(x(i) - 0.3).^2);

    elseif (0.6 < x(i) && x(i) <= 0.8)
        u(i) = 1;
    end
end


%courant number or something or other
nu = c*dt/dx;

%enforce left BC
u(1) = 0;

%t0 for keeping track of time for exact solution
t0 = 0;

%number of time steps
Nt = tend/dt;


%include in denominator to avoid division by 0
eps = 1e-10;

%initialize placeholder for updating the approximate solution
unew = zeros(size(u));

%time stepping
for j=1:Nt
    t0 = t0 + dt;
    %minmode limiter function
    for i = 3:Nx-1
        
        %r at the east edge
        re = (u(i) - u(i-1))/(u(i+1) - u(i) + eps);
        psi_e = psi(re);

        %r at the west edge
        rw = (u(i-1) - u(i-2))/(u(i) - u(i-1) + eps);
        psi_w = psi(rw);


        ue = u(i) + 1/2*psi_e*(u(i+1) - u(i));

        uw = u(i-1) + 1/2*psi_w*(u(i) - u(i-1));

        %update step
        unew(i) = u(i) - nu*(ue - uw);
        % unew(i) = u(i) - nu*(u(i) + psi_e*(u(i+1) - u(i)) - u(i-1) - psi_w*(u(i) - u(i-1)));

    end

    %enforce boundary conditions
    unew(1) = 0;
    unew(2) = 0;
    % unew(end) = 0;
    
    %update u
    u = unew;

    uEx = zeros(size(x));

    %update exact solution
    for k = 1:length(x)
        if (c*t0 <= x(k) && x(k) <= 0.6 + c*t0)
           uEx(k) = exp(-100*(x(k) - 0.3 - c*t0).^2);

        elseif (0.6 + c*t0 < x(k) && x(k) <= 0.8 + c*t0)
            uEx(k) = 1;
        end
    end

    %plot the approximate solution at each time step
    plot(x,u, 'k--', x, uEx, 'b-', 'linewidth', 1.4)
    xlabel('$x$', 'fontsize', 25, 'interpreter', 'latex')
    ylabel('$u(x,2)$', 'fontsize', 25, 'interpreter', 'latex')
    title(methodStr + " scheme for $\frac{\partial u}{\partial t} + c\frac{\partial u}{\partial x} = 0$, $t = 2$" + " $\Delta t = $" + num2str(dt), 'fontsize', 20, 'interpreter', 'latex')
    legend(methodStr + " approx.", "Exact soln.", 'fontsize', 15, 'interpreter', 'latex')
    grid on
    drawnow
  
end

for k = 1:length(x)
        if (c*t0 <= x(k) && x(k) <= 0.6 + c*t0)
           uEx(k) = exp(-100*(x(k) - 0.3 - c*t0).^2);
    
        elseif (0.6 + c*t0 < x(k) && x(k) <= 0.8 + c*t0)
            uEx(k) = 1;
        end
end

plot(x, u, 'k--', x, uEx, 'b-', 'linewidth',1.4)
xlabel('$x$', 'fontsize', 25, 'interpreter', 'latex')
ylabel('$u(x,2)$', 'fontsize', 25, 'interpreter', 'latex')
title(methodStr + " scheme for $\frac{\partial u}{\partial t} + c\frac{\partial u}{\partial x} = 0$, $t = 2$" + " $\Delta t = $" + num2str(dt), 'fontsize', 20, 'interpreter', 'latex')
legend(methodStr + " approx.", "Exact soln.", 'fontsize', 15, 'interpreter', 'latex')
grid on
axis tight

err = dx*norm(u - uEx.',2)