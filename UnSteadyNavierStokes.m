function [unew,r] =UnSteadyNavierStokes(u0,n,Re)
% space-time spectral collocation for stream function 
% formulation of 2D unsteady Navier Stokes equations of the regularized driven cavity flow. 
% by  Sarah Nataj(sarah.nataj@gmail.com).
% Feb 2023
syms x y t
U0(t,x,y)=-0.25*(1-x^2)^2*(1+y)*(1-y^2);
Psiy(t,x) = (1-x^2)^2;
[D, node]=PSDirv(n,1);

[X,Y,T] = meshgrid(node(2:n),node(2:n),node(1:n));
psi0 = double(U0(-1,X,Y));psi0 = psi0(:,:,1); psi0 = reshape(psi0,(n-1)*(n-1),1);
psi1 = double(U0(-1,X,1)); psi1 = psi1(1,:,:); psi1=reshape(psi1,n*(n-1),1);
h = double(Psiy(T,X));h = h(1,:,:);h = reshape(h,n*(n-1),1);

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

dh = D(1:n,n+1);
f = - kron(dh,-kron(B2x,eye(n-1))- kron(eye(n-1),B2y(2:n,2:n)))*psi0...
    - kron(D(1:n,1:n),kron(eye(n-1),B2y(2:n,1)))*psi1;

A=-kron(D(1:n,1:n),kron(B2x,eye(n-1))+kron(eye(n-1),B2y(2:n,2:n)))...
+(1/Re)*kron(eye(n),kron(B4x,eye(n-1)) + kron(eye(n-1),B4y(2:n,2:n))...
    + 2*kron(B2x,B2y(2:n,2:n)));
A = A-diag(kron(eye(n),kron(eye(n-1),D(2:n,2:n)))*u0)*kron(eye(n),kron(B3x,eye(n-1)))...
    + diag(kron(eye(n),kron(D(2:n,2:n),eye(n-1)))*u0)*kron(eye(n),kron(eye(n-1),B3y(2:n,2:n)))...
    + diag(kron(eye(n),kron(D(2:n,2:n),eye(n-1)))*u0)*kron(eye(n),kron(B2x,D(2:n,2:n)))...
    - diag(kron(eye(n),kron(eye(n-1),D(2:n,2:n)))*u0)*kron(eye(n),kron(D(2:n,2:n),B2y(2:n,2:n)));
f = [f; h];
A = [A; kron(eye(n),kron(eye(n-1),D(1,2:n)))];
unew = A\f;
r= unew-u0;


