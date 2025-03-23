%Solve the differential equation u'' = -4pi^2sin(x), u(0) = 0, u(1) = 0 on
%[0,1] via Gauss-Seidel, Jacobi, or Successive Over-Relaxation

clear all; close all;
format long
%Select which method to use: 1 = Jacobi, 2 = Gauss-Seidel, 3 = SOR
method = 1;

%end points of interval
a = 0;
b = 1;

%number of spatial points
Nx = 200;

%step size
h = 1/(Nx - 1);

%spatial vector
x = linspace(a,b,Nx);
x = x.';

% f = -4*pi^2*sin(2*pi*x);
f = 9*exp(3*x);

%initial guess
u0 = 1 + 3*x + 9/2*x.^2;

%max number of iterations
nMax = 4000;

xSmooth = linspace(a,b, 1001);
% uEx = sin(2*pi*xSmooth);
uEx = exp(3*xSmooth);


% uExact = sin(2*pi*x);
uExact = exp(3*x);
errSave = zeros(nMax, 1);
for i = 1:nMax
    switch method

        case 1
            %Jacobi method
            methodType = "Jacobi";
            onesVec = ones(Nx-1,1);
            J = diag(onesVec, -1) + diag(onesVec, 1);
            unew = 1/2*J*u0 - h^2/2*f;
            % unew(1) = 1;
            % unew(end) = exp(3);
            
        case 2
            %Gauss-Seidel
            methodType = "Gauss-Seidel";
            
            onesVec = ones(Nx-1,1);
            twosVec = -2*ones(Nx,1);
            % A = diag(onesVec, 1);
            % B = A';
            % 
            % I = eye(Nx,Nx);

            % RHS = 1/2*B*u0 - h^2/2*f;
            % LHS = I - 1/2*A;
            % unew = LHS\RHS;



            M = 1/h^2*(diag(twosVec) + diag(onesVec,-1));
            N = -1/h^2*diag(onesVec,1);

            unew = inv(M)*N*u0 + inv(M)*f;
            % unew(1) = 0;
            % unew(end) = 0;
            % 
        case 3
            %SOR
            methodType = "SOR";
            omega = 2/(1 + sin(pi*h));
            % omega = 1.5;
            onesVec = ones(Nx-1,1);
            twosVec = -2*ones(Nx,1);

            D = 1/h^2*diag(twosVec);
            L = -1/h^2*diag(onesVec, -1);
            U = L.';

            M = (1/omega)*(D - omega*L);
            N = (1/omega)*((1 - omega)*D + omega*U);
            RHS = N*u0 + f;

            unew = inv(M)*RHS;
            % unew(1) = 0;
            % unew(end) = 0;

    end

    plot(x, u0, 'bo--',xSmooth, uEx, 'k-', 'linewidth', 1.1)
    xlabel('$x$', 'fontsize', 20, 'interpreter', 'latex')
    ylabel('$u(x)$', 'fontsize', 20, 'interpreter', 'latex')
    title("Approx. \& Ex. Soln. to $u'' = -4\pi^2\sin(x)$ via " + methodType, 'fontsize', 16, 'interpreter', 'latex')
    legend('$u_{approx}$', '$u_{exact}$', 'fontsize', 14, 'interpreter', 'latex')
    grid on

    % errSave(i) = norm(unew - uExact, inf);
    % figure(2)
    % semilogy(errSave, 'k--o')
    drawnow

    u0 = unew;
end