function [unew,r] =NonlinHeat(x0,n,epsil,lambda)
% space-time spectral collocation for nonlinear reaction diffusion
% u_t=epsilon u_xx+ lambda exp(u)+f, u(0,t)=u(1,t)=0 u(x,0)=u0.
% by  Sarah Nataj(sarah.nataj@gmail.com).
% Feb 2023
n1=n+1;N=n*(n-1);
[D, node]=PSDirv(n,2);

u0=exp(node(2:n)-1).*cos(pi/2*node(2:n)); % IC
f=-kron(u0,D(1:n,n1));

% [X,T]=meshgrid(node(2:n),node(1:n));
% %UEX=exp(X+T).*cos(pi/2*X);
% f=(1-epsil+epsil*(pi^2)/4)*exp(X+T).*cos(pi/2*X)+epsil*pi*exp(X+T).*sin(pi/2*X)...
%     -lambda*exp(exp(X+T).*cos(pi/2*X));
% u0=exp(node(2:n)-1).*cos(pi/2*node(2:n)); % IC
% %uex=exp(node(2:n)+1).*cos(pi/2*node(2:n)); % exact soln at t=1
% f=reshape(f,N,1)-kron(u0,D(1:n,n1));

A=D*D;
A=kron(eye(n-1),D(1:n,1:n))-epsil*kron(A(2:n,2:n),eye(n));
u=x0;
unew=A\(f+lambda.*exp(u));
r= unew-u;
end

