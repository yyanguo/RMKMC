function [Y_all,flag_zero]=UPdate_Y(G,Pas)
Y_all=sparse(Pas.Nsam,Pas.Nclu);
for j=1:Pas.Nsam
    [~,ind]=max(abs(G(j,:))); %max in abs
    %[~,ind]=max(G(j,:));%max in algebraic
    Y_all(j,ind)=1;
end
% How to do when there are zero columns?
 r=sum(Y_all);
 [~,inx_r]=find(r==0);
 lr=length(inx_r);
%If have zero columns, revise it; else, return.
 if lr>0 
    flag_zero=1;
    [~,inx_m]=max(r); %index of max r 
    [inx_row,~]=find(Y_all(:,inx_m)==1);  
    %r_row=randperm(length(inx_row));
    lw=length(inx_row);
     if lw>=lr
         for j=1:lr
             Y_all(inx_row(j),inx_r(j))=1;
             Y_all(inx_row(j),inx_m)=0;
         end
     else
         disp('zero column')  
       %  [~,inx_sort]=sort(r,'descend');
         error('error algorithm')
     end
 else
    flag_zero=0;
    return
 end
end