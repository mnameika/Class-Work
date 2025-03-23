%streamline test

clear all; close all;

Nx = 50;
Ny = 50;

x = linspace(0,1, Nx);
y = linspace(0,1, Ny);

[Y, X] = meshgrid(y,x);
[Xstr, Ystr] = meshgrid(x,y);

U = X.^2;
V = 2*X.*Y;

quiver(X,Y,U,V,4)
hold on
strm_ln = stream2(X',Y',U,V,0.35,0.244)';
streamline(strm_ln)

strm_ln1 = strm_ln{1};
% strm_ln2 = strm_ln{2};
% plot(strm_ln1(:,1), strm_ln1(:,2))

norm_vec = [-1/sqrt(2);1/sqrt(2)];

reflect_mat = eye(2)  - 2*(norm_vec*norm_vec');

corrected_strm_ln = zeros(length(strm_ln1(:,1)),2);

for i = 1:length(strm_ln1(:,1))
    corrected_strm_ln(i,:) = reflect_mat*(strm_ln1(i,:)');
    
end
plot(corrected_strm_ln(:,1), corrected_strm_ln(:,2))
