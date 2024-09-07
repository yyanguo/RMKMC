function [Obj,SaveAlpha,res,KM_time,iter,SaveBeta]=main_ReKM_SRK_ViewBA(X,Y,CorMat_func,Pas,La_Prel,para,Niu,Gam)
%===========================================================================%
%Overall Algorithm£ºKM + update M, alpha, beta{i}
%X: cell(1,Nview)
%Y: true label
%CorMat_func: @CorMat_Rsamp to calculate C^(v) 
%Pas: Nsam, Nview, Nclu
%La_Prel: initial label
%para: LaSker, SmSker, Ovsam
%===========================================================================%
global Nbaker
ta=tic;
flag=1;
iter=0;
tol_iter=2;
tol_V=0.01;
max_iter=30;
Vo=sparse(Pas.Nsam,Pas.Nclu); 
%----------kernel correlation matrix 
Corker=cell(1,Pas.Nview);
for i=1:Pas.Nview
    Corker{i} = CorMat_func(X{i},Pas); 
end
%----------Initial Beta: All Beta are 1/Nbaker
Beta_old=cell(1,Pas.Nview);
for i=1:Pas.Nview
    Beta_old{i} = ones(Nbaker,1)/Nbaker; %kernel correlation matrix 
end
Alpha_old=ones(Pas.Nview,1)/Pas.Nview;
while flag
    iter=iter+1;
    fprintf('Itering: %d\n',iter)
    %------------------update V
   [V_Blk,V,D_Blk,A_Blk,Ker_Blk,Sig_Blk,US_Blk,ord_Blk]=KM_update_V_ViewBA(X,Alpha_old,Beta_old,La_Prel,Pas,para);
   clear Alpha_old La_Prel
   %-----------------update Alpha
   VD_Blk=cell(1,Pas.Nclu);
    for j=1:Pas.Nclu
        VD_Blk{j}=D_Blk{j}.*V_Blk{j}; %column vector
    end
    clear V_Blk D_Blk
    Alpha=KM_update_Alpha(A_Blk,Sig_Blk,US_Blk,ord_Blk,VD_Blk,Pas,Niu); %vec
    clear A_Blk Sig_Blk US_Blk 
    SaveAlpha(:,iter)=Alpha;
    %----------------update Beta^(v)
    Beta=KM_update_ViewBeta(Corker,Alpha,Ker_Blk,ord_Blk,VD_Blk,Pas,Beta_old,para,Gam);
    clear VD_Blk Ker_Blk ord_Blk Beta_old 
    SaveBeta(:,:,iter)=cell2mat(Beta);
    %----------------Obj value
    Obj(iter)=norm(V-Vo*(Vo'*V),'fro'); 
    clear Vo 
    %----------------clustering measure
    [label_new,~,~]=Re_litekmeans(V,Pas.Nclu,'MaxIter',500,'Replicates',1);
    if iter>tol_iter && (abs((Obj(iter)-Obj(iter-1))/Obj(iter-1))<tol_V || iter>=max_iter)
       flag=0;
       labelr=bestMap(Y,label_new);
       res=ClusteringMeasure_Ave(Y,labelr);
       clear V labelr label_new Alpha Beta Corker
    else
       La_Prel=label_new;
       Vo=V; 
       Alpha_old=Alpha;
       Beta_old=Beta;
       clear V label_new Alpha Beta
    end
end
clear flag
KM_time=toc(ta);
fprintf('PerTime: %f s\n', KM_time);
fprintf('PerResult: %.4f  %.4f  %.4f  %.4f %.4f %.4f %.4f %.4f \n', res);
end