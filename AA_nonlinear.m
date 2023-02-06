function  [ITN,Err,sol]=AA_nonlinear(n,ecase,With_AA,var,var2,w)
% This code calculates the fixed point of a nonlinear system with FP and AA--PP with window size w.
if ecase==1, N=n*(n-1);
elseif ecase==2, N=(n-1)*(n-1);
elseif ecase==3, N=(n-1)*(n-1)*n;
elseif ecase==4, N=n*(n-1);
end

x0 =  rand(N,1);   % initial guess
max_Iter=100-w-2;
Max_Iter =max_Iter;
Tol =1;
Iter_num =0;
%======== standard fixed point method ======================
if With_AA==0
   k=1;
    while (Iter_num < Max_Iter) && Tol>1e-12

        if ecase==1, [solu, r]=NonlinHeat(x0,n,var,var2);
        elseif ecase==2, [solu,r]=NavierStokes(x0,n,var);
        elseif ecase==3, [solu,r]=UnSteadyNavierStokes(x0,n,var);
        elseif ecase==4, [solu,r]=NonlinSch(x0,n,var);
        end
        Tol = norm(r,inf);
        %Tol = norm(r);
        tol(k)=Tol;
        % update solution
        x0 = solu;
        Iter_num = Iter_num+1;
        k = k+1;

    end
    %% record the iteration number and infinity-norm
    ITN = k-1;
    Err=tol(1:ITN);
    sol=solu;
end


if With_AA==1
    % store the update variables
    MAT_r = sparse(N,w+1);
    MAT_sol = cell(1,w+1);

    % set up AA(w) with w=1 inital guess
    for kk=1:w+1
        if ecase==1, [solu, r]=NonlinHeat(x0,n,var,var2);
        elseif ecase==2, [solu,r] =NavierStokes(x0,n,var);
        elseif ecase==3, [solu,r]=UnSteadyNavierStokes(x0,n,var);
        elseif ecase==4, [solu,r]=NonlinSch(x0,n,var);
        end

        % update solution
        MAT_sol{1,kk}=solu;
        MAT_r(:,kk)=x0-solu;
        x0 = solu;
        Tol = norm(r,inf);
        %Tol = norm(r);
        tolp(kk)=Tol;
    end

    k=1;
    while Iter_num < Max_Iter+1 && Tol>1e-12
        %% Anderson acceleration with size w
        % R: r_k- r_{k-1}, r_k -r_{k-2}, ....
        % r_k =x_k-q(x_k)
        R = zeros(N,w);
        for jj=1:w
            R(:,jj) =MAT_r(:,w+1)-MAT_r(:,w+1-jj);
        end
        rr = MAT_r(:,w+1);
        Beta = -R\rr;
        dq=zeros(N,1);

        for kk=1:w
            dq =dq+Beta(kk)*(MAT_sol{1,w+1}-MAT_sol{1,w+1-kk});
        end

        % update the solution by AA(w)
        x0 = MAT_sol{1,w+1}+dq;


        %% next loop AA
        %% MAT_sol(1,1:w)=MAT_sol(1,2:w+1);
        for j=1:w
            MAT_sol{1,j}=MAT_sol{1,j+1};
        end

        %% u =G(u)=U-inv(F)\F;
        if ecase==1, [solu, r]=NonlinHeat(x0,n,var,var2);
        elseif ecase==2, [solu,r] =NavierStokes(x0,n,var);
        elseif ecase==3, [solu,r]=UnSteadyNavierStokes(x0,n,var);
        elseif ecase==4, [solu,r]=NonlinSch(x0,n,var);
        end
        Tol = norm(r,inf);
        %Tol = norm(r);
        tol(k)=Tol;

        % update solution
        MAT_sol{1,w+1}=solu;
        r_p = x0-solu;
        MAT_r =[MAT_r(:,2:end),r_p];
        Iter_num = Iter_num+1;
        k = k+1;

    end

    %% record the iteration number and error
    ITN1 = k-1;
    ITN=(w+1)+ITN1;
    Err=[tolp tol(1:ITN1)];
    sol=solu;

end

end

