function [unew,r] =NavierStokes(u,n,Re)
% spectral collocation for stream function 
% formulation of 2D steady Navier Stokes equations of the regularized driven cavity flow. 
% by  Sarah Nataj(sarah.nataj@gmail.com).
% Feb 2023
syms x 
Psiy = (1-x^2)^2;
Psiy= matlabFunction(Psiy);
N = (n-1)*(n-1);
[D, node]=PSDirv(n,1);
[X,Y] = meshgrid(node(2:n),node(2:n));
h = double(Psiy(X)); h = h(1,:); h=transpose(h); % nonhomo normal deriv on top

%in x direction, Y(x)=(1-x^2)Z(x)
D2 = D^2; D3 = D^3; D4 = D^4;
M = diag(1-node(2:n).^2);
W = diag(node(2:n));
B2x = (M*D2(2:n,2:n) - 4*W*D(2:n,2:n) - 2*eye(n-1))/M;
B3x = (M*D3(2:n,2:n) - 6*W*D2(2:n,2:n) - 6*D(2:n,2:n))/M;
B4x = (M*D4(2:n,2:n) - 8*W*D3(2:n,2:n) - 12*D2(2:n,2:n))/M;

%in y direction, Y(y)=(1+y)Z(y)
M = diag(1+node);
invM = diag([1./(1+node(1:n)); 0]);
B2y = (M*D2 + 2*D)*invM;
B3y = (M*D3 + 3*D2)*invM;
B4y = (M*D4 + 4*D3)*invM;

% linear part
A4 = kron(B4x,eye(n-1)) + kron(eye(n-1),B4y(2:n,2:n)) + 2*kron(B2x,B2y(2:n,2:n));
% nonlinear part
Ayyy = kron(eye(n-1),B3y(2:n,2:n));
Axxx = kron(B3x,eye(n-1));
Ay = kron(eye(n-1),D(2:n,2:n));
Ax = kron(D(2:n,2:n),eye(n-1));
Ayxx = kron(B2x,D(2:n,2:n));
Axyy = kron(D(2:n,2:n),B2y(2:n,2:n));
Dh = kron(eye(n-1),D(1,2:n));

A = (1/Re)*A4 - diag(Ay*u)*Axxx + diag(Ax*u)*Ayyy + ...
    diag(Ax*u)*Ayxx - diag(Ay*u)*Axyy;
fh = zeros(N,1);
fh = [fh; h];
Ah = [A; Dh];
unew = Ah\fh;
r= unew-u;


