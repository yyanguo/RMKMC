function preo_M=InduceLM(M_in,Cold_prod,Half_prod,ord_Blk,Pas)
preb=0;
preo_M=[];
miu=1e-3;
for i=1:Pas.Nclu
    if ord_Blk(i)<=1000
        P=Cold_prod{i}*M_in((preb+1):(ord_Blk(i)+preb),:)+miu*M_in((preb+1):(ord_Blk(i)+preb),:);
    else
        P=Half_prod{i}*(Half_prod{i}'*M_in((preb+1):(ord_Blk(i)+preb),:))+miu*M_in((preb+1):(ord_Blk(i)+preb),:);
    end
    preo_M=[preo_M;P];
    preb=preb+ord_Blk(i);
end
clear P preb miu
end