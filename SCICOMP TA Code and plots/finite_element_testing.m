%quadratic finite element test to solve -u_xx = f(x), f(x) = e^{-x}, u(0) =
%u(1) = 0.
clear all; close all;

%number of elements
N = 3;

basis = 2;

xFine = linspace(0,1,1e5+1);

uEx = -exp(xFine) + (exp(1) - 1)*xFine + 1;

switch basis

    %linear bases
    case 1
        Nx = N + 1;
        x = linspace(0,1,Nx);
        h = x(2) - x(1);

        b = zeros(Nx,1);

        %linear differentiation matrix
        A = 1/h*(2*diag(ones(Nx,1),0) - diag(ones(Nx-1,1),1) - diag(ones(Nx-1,1),-1));

        % linear RHS
        for i = 1:Nx
            b(i) = 1/h*exp(x(i))*(exp(-h) - 2 + exp(h));
        end
        
        A(1,:) = 0;
        A(1,1) = 1;
        A(end,:) = 0;
        A(end,end) = 1;

        
        
        b(1) = 0;
        b(end) = 0;
        
        coeffs = A\b;

        plot(x,coeffs, 'ko--', 'markersize',6)

    %quadratic bases
    case 2
        Nx = 3 + 2*(N-1)x;

        x = linspace(0,1,Nx);
        h = x(2) - x(1);

        b = zeros(Nx,1);

        mainDiag = zeros(Nx,1);
        offDiag1 = -4/3*ones(Nx-1,1);
        offDiag2 = 1/6*ones(Nx-2,1);


        for i=1:Nx
            if (rem(i,2) == 0)
                mainDiag(i) = 8/3;
            else
                mainDiag(i) = 7/3;
            end
        
            if (rem(i,2) == 0 && i < Nx-1)
                offDiag2(i) = 0;
            end
        end

        A = 1/h*(diag(mainDiag,0) + diag(offDiag1,1) + diag(offDiag1,-1) + diag(offDiag2,-2) + diag(offDiag2,2));


        for i=1:Nx
            if (rem(i,2) == 1)
                b(i) = 1/h^2*exp(-x(i))*(-h*(3+cosh(2*h))+2*sinh(2*h));
            else
                b(i) = 4/h^2*exp(-x(i))*(h*cosh(h)-sinh(h));
            end
        end

        A(:,1) = 0;
        A(:,end) = 0;
        A(1,:) = 0;
        A(1,1) = 1;
        A(end,:) = 0;
        A(end,end) = 1;

        
        
        
        
        b(1) = 0;
        b(end) = 0;
        
        coeffs = A\b;
        
        phiSave = quadPhi(xFine,x);
        phiSave = phiSave';

        approxSmooth = phiSave*coeffs;

        plot(x,coeffs, 'ko', xFine, approxSmooth, 'k--', 'linewidth', 1.2, 'markersize', 6)
        

    %cubic bases
    case 3
        Nx = 4 + 3*(N-1);
        x = linspace(0,1,Nx);
        h = x(2) - x(1);

        b = zeros(Nx,1);
        

end





hold on
plot(xFine, uEx, 'b-', 'linewidth', 0.9)
xlabel('$x$', 'fontsize', 25, 'interpreter', 'latex')
ylabel('$u(x)$', 'fontsize', 25, 'interpreter', 'latex')
legend('Approx. Soln.', 'Exact Soln.', 'fontsize', 16, 'interpreter', 'latex')
grid on


%% SUB FUNCTIONS

function phiQSave = quadPhi(xFine,xGrid)

numGrid = length(xGrid);
numFine = length(xFine);

h = xGrid(2) - xGrid(1);

phiQSave = zeros(numGrid, numFine);


%quadratic basis functions
for i = 1:numGrid
    
    %quadratic tent
    if (rem(i,2) == 1)

        for j = 1:numFine

            if (i < numGrid - 1 && xFine(j) <= xGrid(i + 2) && xFine(j) >= xGrid(i))

                phiQSave(i,j) = 1/(2*h^2)*(xFine(j) - xGrid(i) - h)*(xFine(j) - xGrid(i) - 2*h);

            elseif (i > 2 && xFine(j) >= xGrid(i - 2) && xFine(j) <= xGrid(i))

                phiQSave(i,j) = 1/(2*h^2)*(xFine(j) - xGrid(i) + h)*(xFine(j) - xGrid(i) + 2*h);

            end %end if

        end %end for

    %"normal" quadratic
    else
        for j = 1:numFine

            if (xFine(j) >= xGrid(i-1) && xFine(j) <= xGrid(i + 1))
    
                phiQSave(i,j) = -1/h^2*(xFine(j) - xGrid(i)+ h)*(xFine(j) - xGrid(i) - h);

            end %end if
    
        end %end for
        
    end %end if/else


end %end main for loop
    
    
end %end function quadPhi()

function cubicPoly = cubicPhi(xFine, xGrid)
    numFine = length(xFine);
    numGrid = length(xGrid);

    h = xGrid(2) - xGrid(1);

    cubicPoly = zeros(numGrid,numFine);
    phiHold = zeros(1,numFine);




    for j = 1:3:numGrid
        phiHold = zeros(1,numFine);
        for i = 1:numFine
            if (j+3 <= numGrid && xFine(i) >= xGrid(j) && xFine(i) <= xGrid(j+3))
                    phiHold(i) = -1/(6*h^3)*(xFine(i) - xGrid(j+1)).*(xFine(i) - xGrid(j+2)).*(xFine(i) - xGrid(j+3));
                    
            end

            if (j >= 4 && xFine(i) >= xGrid(j-3) && xFine(i) <= xGrid(j))
                    phiHold(i) = 1/(6*h^3)*(xFine(i) - xGrid(j-3)).*(xFine(i) - xGrid(j-2)).*(xFine(i) - xGrid(j-1));
                  
            end
        end
        cubicPoly(j,:) = phiHold;
    end


    %these 'middle' basis functions have a maximum value of not quite 1
    %check this
    for j = 2:3:numGrid
        phiHold = zeros(1,numFine);
        for i=1:numFine
            if (xFine(i) >= xGrid(j-1) && xFine(i) <= xGrid(j+2))
                phiHold(i) = 1/(2*h^3)*(xFine(i) - xGrid(j-1)).*(xFine(i) - xGrid(j+1)).*(xFine(i) - xGrid(j+2));
            end

        end
        cubicPoly(j,:) = phiHold;
    end

    for j = 3:3:numGrid
        phiHold = zeros(1,numFine);
        for i = 1:numFine
            if (xFine(i) >= xGrid(j-2) && xFine(i) <= xGrid(j+1))
                phiHold(i) = -1/(2*h^3)*(xFine(i) - xGrid(j-2)).*(xFine(i) - xGrid(j-1)).*(xFine(i) - xGrid(j+1));
            end

        end
        cubicPoly(j,:) = phiHold;
    end
        
end