% Description: this code is used to compare error of fixed point (FP) iteration
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
% The results have been proposed in the manuscript:
% Anderson acceleration for nonlinear PDEs discretized by space--time spectral methods
% by  Sarah Nataj(sarah.nataj@gmail.com) and Yunhui He.
% Feb 2023
clear all
close all
ecase=1;

if ecase==1       %nonlinear reaction diffusion
    n=15; w_AA=[3,4,5];
    epsilon=0.001;lambda=0.31;
    for i=1:length(w_AA)
        tic
        [ITN(i),Err,solAA]=AA_nonlinear(n,ecase,1,epsilon,lambda,w_AA(i));
        time(i)=toc;
        Err1(i)=Err(end);
    end
elseif ecase==2   %steady Navier Stokes
    n=40;
    w_AA=[5];
    Re=[1000,1100, 1200];
    for i=1:length(Re)
        tic
        [ITN(i),Err,solAA]=AA_nonlinear(n,ecase,1,Re(i),[],w_AA);
        time(i)=toc;
        Err1(i)=Err(end);
    end
elseif ecase==3  %unsteady Navier Stokes
    n=15;
    w_AA=[3];
    Re=[4000,3000,5000];
    for i=1:length(Re)
        tic
        [ITN(i),Err,solAA]=AA_nonlinear(n,ecase,1,Re(i),[],w_AA);
        time(i)=toc;
        Err1(i)=Err(end);
    end
elseif ecase==4   %Nonlinear Sch
    n=20;
    epsilon=1;
    w_AA=[2,3,4];
    for i=1:length(w_AA)
        tic
        [ITN(i),Err,solAA]=AA_nonlinear(n,ecase,1,epsilon,[],w_AA(i));
        time(i)=toc;
        Err1(i)=Err(end);
    end
end
if (ecase==2) || (ecase==3)
    Results=table(Re',ITN',Err1',time','VariableNames',{'Re', 'Iter','Error','Time'})
    %Results=table(AAw',ITN1',Err',time','VariableNames',{'Re', 'Iter','Error','Time'})
else
    Results=table(w_AA',ITN',Err1',time','VariableNames',{'w', 'Iter','Error','Time'})
end


