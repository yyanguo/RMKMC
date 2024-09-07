function alpha=SR_update_Alpha(A_Blk,Sig_Blk,US_Blk,ord_Blk,VD,Pas,Niu)
%==============================================%
%update alpha 
%==============================================%
Df=zeros(Pas.Nview,1);
for i=1:Pas.Nview   
     Deta=0;
     preb=0;
     for j=1:Pas.Nclu
         MB=VD((preb+1):(ord_Blk(j)+preb),:);
         if ord_Blk(j)<=1000
            VKPro=trace(MB'*A_Blk{i}{j}*MB);
            Deta=Deta+VKPro;
            clear VKPro 
         else
             A_PKer=sqrt(diag(Sig_Blk{i}{j})).*US_Blk{i}{j}'*MB;
             VKPro=norm(A_PKer,'fro').^2;
             Deta=Deta+VKPro;
             clear A_PKer VKPro
          end
          preb=preb+ord_Blk(j);
     end
     Df(i)=Deta;
     clear Deta
end
clear MB preb
alpha= EProjSimplex_new(Df/(2*Niu));
clear Df 
end