function alpha=KM_update_Alpha(A_Blk,Sig_Blk,US_Blk,ord_Blk,VD_Blk,Pas,Niu)
%============================================%
%update alpha£¬the combination of the base kernels is inside the A_Blk and US_Blk, Sig_Blk
%only need to compute the vector Df
%=============================================%
Df=zeros(Pas.Nview,1);
for v=1:Pas.Nview 
    Deta=0;
    for j=1:Pas.Nclu
        if ord_Blk(j)<=1000
           VKPro=VD_Blk{j}'*A_Blk{v}{j}*VD_Blk{j};
           Deta=Deta+VKPro;
           clear VKPro
        else 
           A_PKer=sqrt(diag(Sig_Blk{v}{j})).*US_Blk{v}{j}'*VD_Blk{j};
           VKPro=norm(A_PKer)^2;
           Deta=Deta+VKPro;
           clear A_PKer VKPro 
         end
    end
    Df(v)=Deta;
    clear Deta 
end
alpha=EProjSimplex_new(Df/(2*Niu));
clear Df
end