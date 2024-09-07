function [Obj,SaveAlpha,res,SR_time,iter,SaveBeta]=main_ReSR_SRK_ViewBA(X,Y,CorMat_func,Pas,La_Prel,para,Lam,Niu,Gam)
%=========================================================================%
%Overall Algorithm£ºSR + update M, alpha, beta{i}
%X: cell(1,Nview)
%Y: true label
%CorMat_func: @CorMat_Rsamp to calculate C^(v) 
%Pas: Nsam, Nview, Nclu
%La_Prel: initial label
%para: LaSker, SmSker, Ovsam
%=========================================================================%
global Nbaker
ta=tic;
M_init=orth(randn(Pas.Nsam,Pas.Nclu));
Y_init=idx2pm(La_Prel);
flag=1;
iter=0;
tol_iter=2;
tol_M=0.01;
max_iter=30;
%is_zero=0;
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
%--------for block L£¬Li~=Di^(-1/2)*keri*Di^(-1/2)
[A_Blk,Sig_Blk,US_Blk,Ker_Blk,D_Blk,Cold_prod,Half_prod,D_bent,ord_Blk]=SR_VSelK_AD_ViewBA(X,La_Prel,Beta_old,Alpha_old,Pas,para);
clear Alpha_old La_Prel
%--------Y_r : D^(1/2)*Y*(Y'*D*Y)^(-1/2)
Y_rd=sqrt(D_bent).*Y_init;  %D be a column vector
y=sum(D_bent.*Y_init,1);    %sum(*),row vector
Y_r=(1./sqrt(y)).*Y_rd;
clear D_bent Y_init
%--------update eigvector M
[M,~,~]=UPdate_M(Cold_prod,Half_prod,ord_Blk,Pas,M_init,Y_r,Lam);
clear Cold_prod Half_prod 
clear Y_rd Y_r y
%--------update alpha
Dsen=cat(1,D_Blk{:}); %column vector Nsam*1
VD=Dsen.*M; %n*k %diag*M=column vector.*M
clear Dsen D_Blk
Alpha=SR_update_Alpha(A_Blk,Sig_Blk,US_Blk,ord_Blk,VD,Pas,Niu);
clear A_Blk Sig_Blk US_Blk  
SaveAlpha(:,iter)=Alpha;
%----------------update Beta
Beta=SR_update_ViewBeta(Corker,Alpha,Ker_Blk,ord_Blk,VD,Pas,Beta_old,para,Gam);
clear VD Ker_Blk ord_Blk Beta_old
SaveBeta(:,:,iter)=cell2mat(Beta);
%--------update 0-1 matrix Y,sparse
[Yind,~]=UPdate_Y(M,Pas);
%is_zero=is_zero+flag_zero;
%--------update label
prel_new=vec2ind(Yind')';  %turn n*k to n*1
%--------output or not 
Obj(iter)=norm(M-M_init*(M_init'*M),'fro');
clear M_init
if iter>tol_iter && (abs((Obj(iter)-Obj(iter-1))/Obj(iter-1))<tol_M || iter>=max_iter)
   flag=0;
   labelr=bestMap(Y,prel_new);
   res=ClusteringMeasure_Ave(Y,labelr);
   clear prel_new labelr Yind M Alpha Beta
else
   La_Prel=prel_new;
   M_init=M;
   Y_init=Yind;
   Alpha_old=Alpha;
   Beta_old=Beta;
   clear prel_new Yind M Alpha Beta
end
end
clear flag
%Time_SOVLEM=sum(Time_SOVLEM)
%Is_Zero_Column=is_zero;
SR_time=toc(ta);
fprintf('PerTime: %.4f s\n', SR_time);
fprintf('PerResult: %.4f %.4f %.4f %.4f %.4f %.4f %.4f %.4f \n', res);
end
