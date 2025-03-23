%quadratic basis function testing
clear all; close all;
numGrid = 11;
xGrid = linspace(0,1,numGrid);
h = xGrid(2) - xGrid(1);

numFine = 1e3 + 1;
xFine = linspace(0,1,numFine);

phiSave = zeros(numGrid, numFine);

% j = 1;

%quadratic tent function
for i = 1:numGrid
    
    %quadratic tent
    if (rem(i,2) == 1)

        for j = 1:numFine
            %right side of tent
            if (i < numGrid - 1 && xFine(j) <= xGrid(i + 2) && xFine(j) >= xGrid(i))

                phiSave(i,j) = 1/(2*h^2)*(xFine(j) - xGrid(i) - h)*(xFine(j) - xGrid(i) - 2*h);

            %left side of tent
            elseif (i > 2 && xFine(j) >= xGrid(i - 2) && xFine(j) <= xGrid(i))

                phiSave(i,j) = 1/(2*h^2)*(xFine(j) - xGrid(i) + h)*(xFine(j) - xGrid(i) + 2*h);

            end %end if

        end %end for

    %"normal" quadratic
    else
        for j = 1:numFine

            if (xFine(j) >= xGrid(i-1) && xFine(j) <= xGrid(i + 1))
    
                phiSave(i,j) = -1/h^2*(xFine(j) - xGrid(i)+ h)*(xFine(j) - xGrid(i) - h);

            end %end if
    
        end %end for
        
    end %end if/else


end %end main for loop

for i = 1:numGrid
    plot(xFine,phiSave(i,:))
    hold on

end
hold off