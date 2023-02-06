function [unew,r] =NonlinSch(x0,n,epsil)
% space-time spectral collocation for nonlinear Schrodinger equation
% i u_t=-u_xx+ epsil |u|^2 u+f, u(pm 1,t)=0
% by  Sarah Nataj(sarah.nataj@gmail.com).
% Feb 2023
n1=n+1; N=n*(n-1);
[D, node]=PSDirv(n,1);
A=D*D; i=sqrt(-1);
A=i*kron(eye(n-1),D(1:n,1:n))+kron(A(2:n,2:n),eye(n));
[X,T]=meshgrid(node(2:n),node(1:n));
f=exp(X+T).*(sin(pi*X)*(1+i-pi^2)+2*pi*cos(pi*X))-epsil*exp(3*X+3*T).*(sin(pi*X).^3);
u0=exp(node(2:n)-1).*sin(pi*node(2:n)); % IC
%uex=exp(node(2:n)+1).*sin(pi*node(2:n)); % exact soln at t=1
f=reshape(f,N,1)-i*kron(u0,D(1:n,n1));
%f=-i*kron(u0,D(1:n,n1));
u=x0;
unew=(A-epsil*diag(abs(u).^2))\f;
r= unew-u;
end
