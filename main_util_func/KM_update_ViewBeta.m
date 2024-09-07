function Beta=KM_update_ViewBeta(Corker,alpha,Ker_Blk,ord_Blk,VD_Blk,Pas,Beta_old,para,Gam)
%===========================================%
%update Beta for each view£ºNbaker*1
%===========================================%
global Nbaker
Beta=cell(1,Pas.Nview);
Dw=cell(1,Pas.Nview);
for v=1:Pas.Nview   
    Df=zeros(Nbaker,1);    %wi
    for i=1:Nbaker
        if Beta_old{v}(i)>1e-3
            Deta=0;
            for j=1:Pas.Nclu
                if ord_Blk(j)<=1000
                   VKPro=VD_Blk{j}'*Ker_Blk{v}{j}{i}*VD_Blk{j};
                   Deta=Deta+VKPro;
                   clear VKPro
                elseif ord_Blk(j)>1000 && ord_Blk(j)<=50000
                   Nsamker=ceil(para.SmSker*ord_Blk(j));
                   cols = randperm(ord_Blk(j),Nsamker); 
                   ropt=Nsamker-para.Ovsam;   
                   W=Ker_Blk{v}{j}{i}(cols,:);  
                   [Q_ker,R_ker] = qr(Ker_Blk{v}{j}{i},0); 
                   W=R_ker*pinv(W,1e-12)*R_ker'; 
                   W=max(W,W');  
                   [P_W,Sigma,~]=svd(W,0);
                   P_W=P_W(:,1:ropt);
                   Sig=Sigma(1:ropt,1:ropt);
                   US=Q_ker*P_W;
                   A_PKer=sqrt(diag(Sig)).*US'*VD_Blk{j};
                   VKPro=norm(A_PKer)^2;
                   Deta=Deta+VKPro;
                   clear Nsamker cols ropt W Q_ker R_ker P_W Sigma Sig US A_PKer VKPro
                else
                   Nsamker=ceil(para.LaSker*ord_Blk(j));
                   cols = randperm(ord_Blk(j),Nsamker); 
                   ropt=Nsamker-para.Ovsam;   
                   W=Ker_Blk{v}{j}{i}(cols,:);  
                   [Q_ker,R_ker] = qr(Ker_Blk{v}{j}{i},0); 
                   W=R_ker*pinv(W,1e-12)*R_ker'; 
                   W=max(W,W');  
                   [P_W,Sigma,~]=svd(W,0);
                   P_W=P_W(:,1:ropt);
                   Sig=Sigma(1:ropt,1:ropt);
                   US=Q_ker*P_W;
                   A_PKer=sqrt(diag(Sig)).*US'*VD_Blk{j};
                   VKPro=norm(A_PKer)^2;
                   Deta=Deta+VKPro;
                   clear Nsamker cols ropt W Q_ker R_ker P_W Sigma Sig US A_PKer VKPro
                end
             end
             Df(i)=Deta;
             clear Deta
        end
    end
    Dw{v}=alpha(v)*Df;  %Nbaker*1
    clear Df
    Beta{v}=quadprog(Gam*Corker{v},-Dw{v},[],[],ones(1,Nbaker),1,zeros(Nbaker,1),ones(Nbaker,1));
end
clear Dw 
end