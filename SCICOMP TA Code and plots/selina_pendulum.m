clear all; close all;
T = 2*pi; % Total time interval (one period of the pendulum)
m = 20; % Number of intervals
h = T / (m + 1); 
alpha = 0.7; 
beta=0.7;
% theta = alpha*cos(linspace(0, T, m+2))'; % Initial guess for theta
theta = alpha + 0*linspace(0,T,m+2)';
b = zeros(m+2, 1);
 b(1) = alpha; % Boundary condition 
b(end) = beta; 

% Construct matrix A similar to the previous example
onevec = ones(m-1, 1); 
twovec = -2 * ones(m, 1); 
A= 1 / h^2 * (diag(onevec, -1) + diag(twovec, 0) + diag(onevec, 1));

tol = 1e-10; 
max_iter = 100; 
plot(theta)
drawnow
pause(0.5)
for k = 1:max_iter
 
    G=1/h^2*(theta(1:end-2)-2*theta(2:end-1)+theta(3:end))+sin(theta(2:end-1));

    % Construct the Jacobian matrix J
    J = A + diag(cos(theta(2:end-1)));
    delta_theta = -J \ G;
    theta(2:end-1) = theta(2:end-1) + delta_theta;%solution update)

    % Check for convergence
    if norm(delta_theta, inf) < tol
        disp(['Converged in ', num2str(k), ' iterations']);
        break;
    end
    plot(theta)
    drawnow
    pause(0.5)
end

% Plot the solution
x = linspace(0, T, m+2); % Time points
plot(x, theta, '-o');
xlabel('Time t');
ylabel('\theta(t)');
title('Solution to Nonlinear Pendulum BVP');