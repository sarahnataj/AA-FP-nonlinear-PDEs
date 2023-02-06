% Description: this code is used to plot error of fixed point (FP) iteration
% and fixed point iteration with Anderson Acceleration of window size w(AA-FP) for
% nonlinear system arising from spectral doscretization of nonlinear PDEs including
% nonlinear reaction diffusion, nonlinear Schrodinger and  incompressible Navier Stokes equations
%
% This code calls AA_nonlinear.m which forms the AA-FP and FP methods
%                 [ITN,Err,sol]=AA_nonlinear(n,ecase,With_AA,var,var2,w)
%                 INPUT n: the number of nodes, 
%                 ecase: test problem
%                 With_AA: with or without AA, WithAA=1 for AA--FP; WithAA=0 for FP
%                 var and var2: problem physical parameters
%                 w: AA window size
%                 OUTPUT  ITN: iteration number; Err: residual error; sol: numerical solution
%
% This consider AA(w), where w is the windowsize and w is often very small, for example, w=1, 2, 3, 4, 5...
%
% The results have been proposed in the manuscript:
% Anderson acceleration for nonlinear PDEs discretized by space--time spectral methods
% by  Sarah Nataj(sarah.nataj@gmail.com) and Yunhui He.
% Feb 2023

clear
close all
ecase=1;   %choose between 1,2,3,4

if ecase==1  %Nonlinear Heat
    n=60;
    epsilon=0.0001;lambda=0.32;
    [ITN0,Err0,sol]=AA_nonlinear(n,ecase,0,epsilon,lambda,0);
    [ITN1,Err1,solAA]=AA_nonlinear(n,ecase,1,epsilon,lambda,5);

    if (IsNAN(sol)==1) && (IsNAN(solAA)==1)
        [D, node]=PSDirv(n,1);
        u=reshape(sol,n,n-1)';
        uAA=reshape(solAA,n,n-1)';
        figure(2), plot( node,[0;uAA(:,1);0],'k-');
        title('Numerical solution, AA-FP','interpreter','latex','fontsize', 12, 'fontweight','bold')
        plotformat(1.5,6)
    end
    plot2(ITN0,Err0, ITN1,Err1)
elseif ecase==2   %steady Navier Stokes
    n=31;
    Re=1000;
    [ITN0,Err0,sol]=AA_nonlinear(n,ecase,0,Re,[],0);
    [ITN1,Err1,solAA]=AA_nonlinear(n,ecase,1,Re,[],5);
    if (IsNAN(sol)==1) && (IsNAN(solAA)==1)
        [D, node]=PSDirv(n,1);
        [Xh,Yh] = meshgrid(node(2:n),node(2:n));
        UN = reshape(solAA,n-1,n-1);
        figure(2); surf(Xh,Yh,UN); title('Numerical solution, AA-FP','interpreter','latex','fontsize', 12, 'fontweight','bold');
        plotformat(1.5,6)
        UN = [zeros(1,n+1); zeros(n-1,1) UN zeros(n-1,1); zeros(1,n+1)];
        [Xh,Yh] = meshgrid(node,node);
        [Xq,Yq] = meshgrid(-1:0.01:1);
        Uq = interp2(Xh,Yh,UN,Xq,Yq,'spline');
        figure(3); contour(Xq,Yq,Uq, min(min(Uq)):0.003:max(max(Uq)));
        title('Stream function','interpreter','latex','fontsize', 12, 'fontweight','bold');colormap(jet);
        plotformat(1.5,6)
    end
    plot2(ITN0,Err0, ITN1,Err1)
elseif ecase==3   %unsteady Navier Stokes
    n=15;
    Re=5000;
    %[ITN0,Err0,sol]=AA_nonlinear(n,ecase,0,Re,[],0);
    [ITN1,Err1,solAA]=AA_nonlinear(n,ecase,1,Re,[],3);
    %if (IsNAN(sol)==1) && (IsNAN(solAA)==1)
    if (IsNAN(solAA)==1)
        [D, node]=PSDirv(n,1);
        [Xh,Yh,Th] = meshgrid(node(2:n),node(2:n),node(1:n));
        Xh=Xh(:,:,1);Yh=Yh(:,:,1);
        UN = reshape(solAA,n-1,n-1,n);UN=UN(:,:,1);
        figure(2); surf(Xh,Yh,UN);title('Numerical solution, AA-FP','interpreter','latex','fontsize', 12, 'fontweight','bold');
        UN = [zeros(1,n+1); zeros(n-1,1) UN zeros(n-1,1); zeros(1,n+1)];
        [Xh,Yh] = meshgrid(node,node);
        [Xq,Yq] = meshgrid(-1:0.01:1);
        Uq = interp2(Xh,Yh,UN,Xq,Yq,'spline');
        figure(3); contour(Xq,Yq,Uq, min(min(Uq)):0.003:max(max(Uq)));
        title('Stream function','interpreter','latex','fontsize', 12, 'fontweight','bold');colormap(jet);
        plotformat(1.5,6)
    end
    plot1(ITN1,Err1)
elseif ecase==4   %Nonlinear Sch
    n=20;
    epsilon=1;
    [ITN0,Err0,sol]=AA_nonlinear(n,ecase,0,epsilon,[],0);
    [ITN1,Err1,solAA]=AA_nonlinear(n,ecase,1,epsilon,[],4);
    if (IsNAN(sol)==1) && (IsNAN(solAA)==1)
        [D, node]=PSDirv(n,1);
        u=reshape(sol,n,n-1)';
        uAA=reshape(solAA,n,n-1)';
        figure(2), plot(node,[0;uAA(:,1);0],'k-');
        title('Numerical solution, AA-FP','interpreter','latex','fontsize', 12, 'fontweight','bold')
        plotformat(1.5,6)
    end
    plot2(ITN0,Err0, ITN1,Err1)
end

function nan=IsNAN(x)
NN=size(x,1);
if isnan(x)==zeros(NN,1)
    nan=1;
else
    nan=0;
end
end

function plot2(ITN0,Err0, ITN1,Err1)
figure (1), semilogy(1:ITN0,Err0,'k-d','LineWidth',1)
hold on
semilogy(1:ITN1,Err1,'r-o','LineWidth',1)
xlabel('Iteration number','interpreter','latex','fontsize', 14,'fontweight','bold')
ylabel('$\|r_k\|_{\infty}$','interpreter','latex','fontsize', 14,'fontweight','bold')
legend('FP', 'AA-FP','interpreter','latex','fontsize', 12, 'fontweight','bold', 'Location','best')
title('AA-FP vs. FP','interpreter','latex','fontsize', 14, 'fontweight','bold')
plotformat(1.5,6)
axis([0 100, 1e-13 20])
hold off
end

function plot1(ITN1,Err1)
figure (1),  semilogy(1:ITN1,Err1,'r-o','LineWidth',1)
xlabel('Iteration number','interpreter','latex','fontsize', 14,'fontweight','bold')
ylabel('$\|r_k\|_{\infty}$','interpreter','latex','fontsize', 14,'fontweight','bold')
title('AA-FP','interpreter','latex','fontsize', 14, 'fontweight','bold')
plotformat(1.5,6)
axis([0 100, 1e-13 20])
end




