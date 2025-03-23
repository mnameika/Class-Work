%solve nonlinear Duffing equation with (small) diffusion term:
% y'' + y - eps(y^3 + y') = 0.
clear all; close all;
y0 = 1;
g0 = 0;

eps = 0.01;
Nt = 10000;
t = linspace(0,200,Nt+1);


F = @(y) [y(2); -y(1) + eps*y(1).^3 + eps*y(2)];

h = t(2) - t(1);

% Nt = t(end)/dt;

ysoln = zeros(2,Nt+1);
ysoln(1,1) = y0;
ysoln(2,1) = g0;



for i=1:Nt
    tn = t(i+1);

    k1 = F(ysoln(:,i));
    k2 = F(ysoln(:,i) + h*k1/2);
    k3 = F(ysoln(:,i) + h*k2/2);
    k4 = F(ysoln(:,i) + h*k3);

    ysoln(:,i+1) = ysoln(:,i) + h/6*(k1 + 2*k2 + 2*k3 + k4);
    % gsoln(i+1) = solnUpd(2);
end

plot(t,ysoln(1,:), 'k-', 'linewidth', 1.1)

yasmp = y0*exp(eps*t/2).*cos(t + y0^2/8*(1 - exp(eps*t)));
hold on
plot(t,yasmp, 'b-', 'linewidth', 1.1)
hold off
grid on
xlabel('$t$', 'fontsize', 25, 'interpreter', 'latex')
ylabel('$y(t)$', 'fontsize', 25, 'interpreter', 'latex')
legend('Numerical Solution', 'Asymptotic Approximation', 'fontsize',20, 'interpreter', 'latex', 'location', 'northwest')
title("Numerical and Approximate solutions to $\frac{d^2y}{dt^2} + y - \varepsilon\left(y^3 +" + ...
    " \frac{dy}{dt}\right) = 0$ for $\varepsilon = $ " + num2str(eps), 'fontsize', 25, 'interpreter', 'latex')
