clear all; close all;
perhaps_img = imread('C:\Users\Michael\Desktop\perhaps.png');

perhaps_img = double(perhaps_img(:,1:end-1,3))/255;

% imshow(perhaps_img)

dims = size(perhaps_img);

uEx = zeros(dims + 2);

uEx(2:end-1, 2:end-1) = perhaps_img;


luDiag = ones(dims(1)-1, 1);

mDiag = ones(dims(1),1);

h = 1/(dims(1)-1);

D2 = 1/h^2*spdiags([mDiag -2*mDiag mDiag], [-1 0 1], dims(1), dims(2));

I = speye(dims(1),dims(2));

Lop = kron(D2,I) + kron(I,D2);


imVec = reshape(perhaps_img, dims(1)^2,1);

lapIm = Lop*imVec;
% lapIm = reshape(lapIm, dims(1), dims(2));
% 
% imshow(lapIm)

% approxIm = Lop\lapIm;
% approxIm = reshape(approxIm, dims(1),dims(2));
omega = 1.8;

D = diag(diag(Lop));

% D1 = inv(D);



L = tril(Lop) - D;
U = L';


LHS = D + omega*L;

tol = 1e-4;

maxCount = 5000;

x0 = ones(dims(1)^2,1);

for i=1:maxCount
    RHS = -(omega*U + (omega - 1)*D)*x0 + omega*lapIm;
    xnew = LHS\RHS;
    
    imshow(reshape(xnew,dims(1),dims(2)));

    x0 = xnew;
    % i
end

imshow(reshape(x0,dims(1),dims(2)))