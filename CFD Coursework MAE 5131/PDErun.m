clear all; close all;

dt = 0.001;
dx = 0.001;

x0 = 0;
xend = 1;

t0 = 0;
tend = 0.2;

tVec = 0:dt:tend;

N = 100;

x = x0:dx:xend;
x(end) = [];
x = x.';

u0 = sin(2*pi*x);

Nx = floor((xend - x0)/dx);

Nt = floor((tend - t0)/dt);
u = zeros(Nx, 1);

plotMat = zeros(Nt+1,Nx);
plotMat(1,:) = u0;

intNum = 0;
t = t0;
for i=1:Nt
    t = t + dt;

    u = zeros(Nx, 1);
    for j = 1:N
        %first numerically approximate the integral
        intNum = dx*sum(sin(2*pi*x).*sin(j*pi*x));

        u = u + 2*(intNum*exp(-j^2*pi^2*t)*sin(j*pi*x));

    end
    plotMat(i+1,:) = u;
    
    
    
end

[X,T] = meshgrid(x,tVec);

mesh(X,T,plotMat)
xlabel('$x$', 'fontsize', 20, 'interpreter', 'latex')
ylabel('$t$', 'fontsize', 20, 'interpreter', 'latex')
zlabel('$u(x,t)$', 'fontsize', 20, 'interpreter', 'latex')
title('Numerically approximated solution ($\delta t = 0.001 = \delta x$)', 'fontsize', 16, 'interpreter', 'latex')
colorbar