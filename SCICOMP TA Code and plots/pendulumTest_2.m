clear all; close all;
a = 0;
b = 200*pi;
% b = 4*pi;

N = 8*64;

h = (b - a)/(N + 1);

T = linspace(a,b,N+2);
T = T.';

alpha = 0.7;
beta = 0.7;

u0 = alpha*cos(T/1000);
% u0 = alpha + 0*T;


tol = 1e-6;
maxCount = 1000;
err = 10;
count = 0;

onesVec = ones(N-1,1);
twosVec = -2*ones(N,1);

% A = 1/h^2*(diag(onesVec, -1) + diag(twosVec, 0) + diag(onesVec, 1));

%want to use Newton's method to hit the interior nodes.
file = 'newton.gif';
for j = 1:5
plot(T, u0)
xlabel('$t$', 'fontsize', 25, 'interpreter', 'latex')
    ylabel('$\Theta_n(t)$', 'fontsize', 25, 'interpreter', 'latex')
    title("Initial Guess $\Theta_0(t/20) = \alpha\cos(t)$ " + "(Iteration \# " + num2str(count) + ")", 'fontsize', 15, 'interpreter', 'latex')
    axis tight
    
drawnow
exportgraphics(gca,'newton.gif', 'Append', true)
end

while (count < maxCount && err > tol)
    J = 1/h^2*(diag(onesVec,-1) + diag(twosVec) + h^2*diag(cos(u0(2:end-1))) + diag(onesVec,1));
    G = 1/h^2*(u0(1:end-2) -2*u0(2:end-1) + u0(3:end)) + sin(u0(2:end-1));
    
    delta = -J\G;

    u0(2:end-1) = u0(2:end-1) + delta;

    err = norm(delta,inf)

    count = count + 1;
    for j=1:5
        plot(T,u0)
        xlabel('$t$', 'fontsize', 25, 'interpreter', 'latex')
        ylabel('$\Theta_n(t)$', 'fontsize', 25, 'interpreter', 'latex')
        title("Initial Guess $\Theta_0(t/20) = \alpha\cos(t)$ " + "(Iteration \# " + num2str(count) + ")", 'fontsize', 15, 'interpreter', 'latex')
        axis tight
        drawnow
        exportgraphics(gca,'newton.gif', 'Append', true)
    
    end
    
end

for i=1:30
    plot(T,u0)
    xlabel('$t$', 'fontsize', 25, 'interpreter', 'latex')
    ylabel('$\Theta_n(t)$', 'fontsize', 25, 'interpreter', 'latex')
    title("Initial Guess $\Theta_0(t/20) = \alpha\cos(t)$ " + "(Iteration \# " + num2str(count) + ")", 'fontsize', 15, 'interpreter', 'latex')
   
    axis tight
    drawnow
    exportgraphics(gca,'newton.gif', 'Append', true)

end

