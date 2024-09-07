function Beta=SR_update_ViewBeta(Corker,alpha,Ker_Blk,ord_Blk,VD,Pas,Beta_old,para,Gam)
%=======================================%
%update Beta for each view£ºNbaker*1
%=======================================%
global Nbaker
Beta=cell(1,Pas.Nview);
Dw=cell(1,Pas.Nview);
for i=1:Pas.Nview   
    Df=zeros(Nbaker,1);    %wi
    for ll=1:Nbaker
        if Beta_old{i}(ll)>1e-3
            Deta=0;
            preb=0;
            for j=1:Pas.Nclu
                VDD=VD((preb+1):(ord_Blk(j)+preb),:);
                if ord_Blk(j)<=1000
                    VKPro=trace(VDD'*Ker_Blk{i}{j}{ll}*VDD); 
                    Deta=Deta+VKPro;
                    clear VKPro VDD
                elseif ord_Blk(j)>1000 && ord_Blk(j)<=50000
                    Nsamker=ceil(para.SmSker*ord_Blk(j));
                    cols = randperm(ord_Blk(j),Nsamker); 
                    ropt=Nsamker-para.Ovsam;
                    W=Ker_Blk{i}{j}{ll}(cols,:);
                    [Q_ker,R_ker] = qr(Ker_Blk{i}{j}{ll},0); 
                    W=R_ker*pinv(W,1e-12)*R_ker'; 
                    W=max(W,W'); 
                    [P_W,Sigma,~]=svd(W,0);
                    P_W=P_W(:,1:ropt);
                    Sig=Sigma(1:ropt,1:ropt);
                    US=Q_ker*P_W;
                    A_PKer=sqrt(diag(Sig)).*US'*VDD;
                    VKPro=norm(A_PKer,'fro').^2;
                    Deta=Deta+VKPro;
                    clear Nsamker cols ropt W Q_ker R_ker P_W Sigma Sig US A_PKer VKPro VDD
                else
                    Nsamker=ceil(para.LaSker*ord_Blk(j));
                    cols = randperm(ord_Blk(j),Nsamker); 
                    ropt=Nsamker-para.Ovsam;
                    W=Ker_Blk{i}{j}{ll}(cols,:);
                    [Q_ker,R_ker] = qr(Ker_Blk{i}{j}{ll},0); 
                    W=R_ker*pinv(W,1e-12)*R_ker'; 
                    W=max(W,W'); 
                    [P_W,Sigma,~]=svd(W,0);
                    P_W=P_W(:,1:ropt);
                    Sig=Sigma(1:ropt,1:ropt);
                    US=Q_ker*P_W;
                    A_PKer=sqrt(diag(Sig)).*US'*VDD;
                    VKPro=norm(A_PKer,'fro').^2;
                    Deta=Deta+VKPro;
                    clear Nsamker cols ropt W Q_ker R_ker P_W Sigma Sig US A_PKer VKPro VDD
                end
             preb=preb+ord_Blk(j);
             end
             Df(ll)=Deta;
             clear Deta preb
        end
    end
    Dw{i}=alpha(i)*Df;  %Nbaker*1
    clear Df
    Beta{i}=quadprog(Gam*Corker{i},-Dw{i},[],[],ones(1,Nbaker),1,zeros(Nbaker,1),ones(Nbaker,1));
end
clear Dw Corker
end