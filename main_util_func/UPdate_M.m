function [M_new,Time_SOVLEM,iter]=UPdate_M(Cold_prod,Half_prod,ord_Blk,Pas,M_in,Y_r,Lam)
%=================================%
% SPD
%=================================%
Z_r=Y_r;
tsol=tic;
iter=0;
flag=1;
tol_iter=2;
tol=0.1;
max_iter=20;
while flag
    iter=iter+1;
    preo_M=InduceLM(M_in,Cold_prod,Half_prod,ord_Blk,Pas);
    %NAPM=numel(find(isnan(preo_M)))
    Z=preo_M+0.5*Lam*Z_r;
    clear preo_M 
    [U,~,V]=svd(Z,0);
    M_new=U*V';
    % preo_Mnew=InduceLM(M_new,Cold_prod,Half_prod,ord_block,num_clu);
    % Z_new=preo_Mnew+lambda*Y_r;
    % Obj(it)=trace(M_new'*Z_new);
    Obj(iter)=norm(M_new-M_in*(M_in'*M_new),'fro');
    if iter>tol_iter && (abs((Obj(iter)-Obj(iter-1))/Obj(iter-1))<tol ||iter>=max_iter) 
       flag=0;
       clear Z U V Obj M_in preo_Mnew
    else
       M_in=M_new;
       clear Z U V M_new 
    end
end
clear flag
%plot(Obj)
Time_SOVLEM=toc(tsol);
end

