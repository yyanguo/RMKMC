function [V_Blk,V,D_Blk,A_Blk,Ker_Blk,Sig_Blk,US_Blk,ord_Blk]=KM_update_V_ViewBA(X,Alpha_old,Beta_old,La_Prel,Pas,para)
%===========================================%
%X£ºcell(1,Nview)   
%Alpha_old£ºNview*1
%Beta_old: Nview*1
%La_Prel: initial label 
%Pas, para: parameters
%===========================================%
global KernelTypes KernelPostProcessTypes Degrees Nbaker
V_Blk=cell(1,Pas.Nclu);
D_Blk=cell(1,Pas.Nclu);
A_Blk=cell(1,Pas.Nview);
Ker_Blk=cell(1,Pas.Nview);
Sig_Blk=cell(1,Pas.Nview);
US_Blk=cell(1,Pas.Nview);
ord_Blk=zeros(1,Pas.Nclu);
sub_kernel=cell(1,Pas.Nclu);
ind=cell(1,Pas.Nclu);
Sropt=zeros(1,Pas.Nclu);
Lropt=zeros(1,Pas.Nclu);
%the size and corresponding position of each block are the same for each view
for j=1:Pas.Nclu 
   [ind{j},~]=find(La_Prel==j);
   ord_Blk(j)=numel(ind{j}); %order number of this block
end
%--------construct block diagnoal multiple kernel matrix K^(v)
for v=1:Pas.Nview
    for j=1:Pas.Nclu
       view_subfea=X{v}(ind{j},:);
       %K_j^(V) is calculated in four cases
       if ord_Blk(j)==1 %
           A_Blk{v}{j}=1;
           for i=1:Nbaker
             if Beta_old{v}(i)>1e-3
                 Ker_Blk{v}{j}{i}=1;
             end
          end
       elseif ord_Blk(j)>1 && ord_Blk(j)<=1000  %
          sub_ker=zeros(ord_Blk(j),ord_Blk(j)); 
          for i=1:Nbaker
              if Beta_old{v}(i)>1e-3
                  kernel_option = [];
                  switch lower(KernelTypes{i})
                     case lower('Linear')
                         kernel_option.KernelType = 'Linear';
                         for iPost = KernelPostProcessTypes
                             sub_kernel = constructKernel(view_subfea,[], kernel_option);
                             sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                         end
                     case lower('Polynomial')
                         kernel_option.KernelType = 'Polynomial';
                         kernel_option.d = Degrees(i);
                         for iPost = KernelPostProcessTypes
                            sub_kernel = constructKernel(view_subfea,[], kernel_option);   
                            sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                         end
                     case lower('PolyPlus')
                         kernel_option.KernelType = 'PolyPlus';
                         kernel_option.d = Degrees(i);
                         for iPost = KernelPostProcessTypes
                            sub_kernel = constructKernel(view_subfea,[], kernel_option);
                            sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                         end
                     case lower('Gaussian')
                         kernel_option.KernelType = 'Gaussian';
                         kernel_option.t = Degrees(i);
                         for iPost = KernelPostProcessTypes
                            sub_kernel = constructKernel(view_subfea,[], kernel_option);
                            sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                         end
                      otherwise
                         error('KernelType does not exist!');
                   end
                   sub_ker=sub_ker+Beta_old{v}(i)*sub_kernel; 
                   Ker_Blk{v}{j}{i}=sub_kernel;   
                   clear sub_kernel 
              end
          end
          A_Blk{v}{j} = sub_ker;  %a combined kernel matrix for each block per view
          clear view_subfea sub_ker kernel_option
       elseif ord_Blk(j)>1000 && ord_Blk(j)<=50000   %approximate
          Nsamker=ceil(para.SmSker*ord_Blk(j));
          cols = randperm(ord_Blk(j),Nsamker); 
          Rview_fea=view_subfea(cols,:);
          sub_ker=zeros(ord_Blk(j),Nsamker);
          for i=1:Nbaker
              if Beta_old{v}(i)>1e-3
                 kernel_option = [];
                 switch lower(KernelTypes{i})
                     case lower('Linear')
                         kernel_option.KernelType = 'Linear';
                         for iPost = KernelPostProcessTypes
                            sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                            sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                         end
                      case lower('Polynomial')
                         kernel_option.KernelType = 'Polynomial';
                         kernel_option.d = Degrees(i);
                         for iPost = KernelPostProcessTypes
                            sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                            sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                         end
                     case lower('PolyPlus')
                        kernel_option.KernelType = 'PolyPlus';
                        kernel_option.d = Degrees(i);
                        for iPost = KernelPostProcessTypes
                           sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                           sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                        end
                    case lower('Gaussian')
                       kernel_option.KernelType = 'Gaussian';
                       kernel_option.t = Degrees(i);
                       for iPost = KernelPostProcessTypes
                          sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                          sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                       end
                   otherwise
                   error('KernelType does not exist!');
                end
                sub_ker=sub_ker+Beta_old{v}(i)*sub_kernel; %C after combination
                Ker_Blk{v}{j}{i}=sub_kernel;
                clear sub_kernel
             end
          end
         Sropt(j)=Nsamker-para.Ovsam;
         W=sub_ker(cols,:);  
         [Q_ker,R_ker] = qr(sub_ker,0); 
         W=R_ker*pinv(W,1e-12)*R_ker'; 
         W=max(W,W');
         [P_W,Sigma,~]=svd(W,0);
         P_W=P_W(:,1:Sropt(j));
         Sig_Blk{v}{j}=Sigma(1:Sropt(j),1:Sropt(j));
         US_Blk{v}{j}=Q_ker*P_W;
         clear sub_ker view_subfea Rview_fea Nsamker cols ker W Q_ker R_ker P_W Sigma kernel_option
     else 
        Nsamker=ceil(para.LaSker*ord_Blk(j));
        cols = randperm(ord_Blk(j),Nsamker); 
        Rview_fea=view_subfea(cols,:);
        sub_ker=zeros(ord_Blk(j),Nsamker);
        for i=1:Nbaker
            if Beta_old{v}(i)>1e-3
               kernel_option = [];
              switch lower(KernelTypes{i})
                   case lower('Linear')
                        kernel_option.KernelType = 'Linear';
                       for iPost = KernelPostProcessTypes
                           sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                           sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                       end
                   case lower('Polynomial')
                       kernel_option.KernelType = 'Polynomial';
                       kernel_option.d = Degrees(i);
                       for iPost = KernelPostProcessTypes
                           sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option); 
                           sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                       end
                  case lower('PolyPlus')
                      kernel_option.KernelType = 'PolyPlus';
                      kernel_option.d = Degrees(i);
                      for iPost = KernelPostProcessTypes
                          sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                          sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                      end
                 case lower('Gaussian')
                      kernel_option.KernelType = 'Gaussian';
                      kernel_option.t = Degrees(i);
                      for iPost = KernelPostProcessTypes
                          sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option); 
                          sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                      end
                 otherwise
                  error('KernelType does not exist!');
              end
              sub_ker=sub_ker+Beta_old{v}(i)*sub_kernel; %C after combination
              Ker_Blk{v}{j}{i}=sub_kernel;
              clear sub_kernel
            end
        end
        Lropt(j)=Nsamker-para.Ovsam;   
        W=sub_ker(cols,:);  
        [Q_ker,R_ker] = qr(sub_ker,0); 
        W=R_ker*pinv(W,1e-12)*R_ker'; 
        W=max(W,W');  
        [P_W,Sigma,~]=svd(W,0);
        P_W=P_W(:,1:Lropt(j));
        Sig_Blk{v}{j}=Sigma(1:Lropt(j),1:Lropt(j));
        US_Blk{v}{j}=Q_ker*P_W;
        clear sub_ker view_subfea Rview_fea Nsamker cols ker W Q_ker R_ker P_W Sigma kernel_option
       end     
   end 
end 
clear ind
%---------------compute the largest eigenpair for each block                                                                                                                                                                                                                
for j=1:Pas.Nclu  
    if ord_Blk(j)==1   
        D_Blk{j}=1;
        V_Blk{j}=1;
    elseif ord_Blk(j)>1 && ord_Blk(j)<=1000 %combination
       Abent_Blk=zeros(ord_Blk(j));
       for v=1:Pas.Nview
            Abent_Blk=Abent_Blk+Alpha_old(v)*A_Blk{v}{j};
       end
       Dbent_Blk=sum(Abent_Blk,2); %column vector
       for rd=1:ord_Blk(j)
           if Dbent_Blk(rd)<=1e-12
              Dbent_Blk(rd)=1;
           end
       end
       clear rd
       r=1./sqrt(Dbent_Blk);
       D_Blk{j}=r;
       Cold_prod=r.*(r'.*Abent_Blk);
       Cold_prod=max(Cold_prod,Cold_prod');
       [eigvec,eigval]=eig(Cold_prod); 
       [~,inx]=sort(diag(eigval),'descend');
       V_Blk{j}=eigvec(:,inx(1));  %largest eigenvector 
       clear Abent_Blk Dbent_Blk r Cold_prod eigvec eigval inx
    elseif ord_Blk(j)>1000 && ord_Blk(j)<=50000   %approximate
       Ubt=[];
       Dbt=zeros(ord_Blk(j),1);
       SiGb=cell(1,Pas.Nview);
       for v=1:Pas.Nview
           Ubt=[Ubt US_Blk{v}{j}];
           SiGb{v}=Alpha_old(v)*Sig_Blk{v}{j};
           Dbt=Dbt+abs(US_Blk{v}{j}*(SiGb{v}*(US_Blk{v}{j}'*ones(ord_Blk(j),1))));
       end  
       SiGbt=sparse(blkdiag(SiGb{:}));
       for rd=1:ord_Blk(j)
          if Dbt(rd)<=1e-12
             Dbt(rd)=1;
          end
       end
       clear rd SiGb
       UU=1./sqrt(Dbt);
       D_Blk{j}=UU;
       U=UU.*Ubt;
       [Q_U,R_U]=qr(U,0);
       Sropt(j)=Sropt(j)*Pas.Nview;
       if Sropt(j)<=1000
          P_U=R_U*SiGbt*R_U'; 
          P_U=max(P_U,P_U');
          [eigvec,eigval]=eig(P_U);
          eigvec=Q_U*eigvec;   
          [~,inx]=sort(diag(eigval),'descend');
          V_Blk{j}=eigvec(:,inx(1));  %largest one
       else
          [eigvec,~]=eigs(@(x)RSRproduct(R_U,SiGbt,x),Sropt(j),1,'la','Tolerance',1e-3,'MaxIterations',1e+4);
          V_Blk{j}=Q_U*eigvec;
       end  
       clear SiGbt Ubt Dbt U UU R_U Q_U P_U eigvec eigval inx
    else
       Ubt=[];
       Dbt=zeros(ord_Blk(j),1);
       SiGb=cell(1,Pas.Nview);
       for v=1:Pas.Nview
           Ubt=[Ubt US_Blk{v}{j}];
           SiGb{v}=Alpha_old(v)*Sig_Blk{v}{j};
           Dbt=Dbt+abs(US_Blk{v}{j}*(SiGb{v}*(US_Blk{v}{j}'*ones(ord_Blk(j),1))));
       end  
       SiGbt=sparse(blkdiag(SiGb{:}));
       for rd=1:ord_Blk(j)
          if Dbt(rd)<=1e-12
             Dbt(rd)=1;
          end
       end
       clear rd SiGb
       UU=1./sqrt(Dbt);
       D_Blk{j}=UU;
       U=UU.*Ubt;
       [Q_U,R_U]=qr(U,0); 
       Lropt(j)=Lropt(j)*Pas.Nview;
       if Lropt(j)<=1000
          P_U=R_U*SiGbt*R_U'; 
          P_U=max(P_U,P_U');
          [eigvec,eigval]=eig(P_U);
          eigvec=Q_U*eigvec;   
          [~,inx]=sort(diag(eigval),'descend');
          V_Blk{j}=eigvec(:,inx(1));  %largest one
       else
          [eigvec,~]=eigs(@(x)RSRproduct(R_U,SiGbt,x),Lropt(j),1,'la','Tolerance',1e-3,'MaxIterations',1e+4);
          V_Blk{j}=Q_U*eigvec;
        end  
        clear SiGbt Ubt Dbt U UU R_U Q_U P_U eigvec eigval inx
     end
end
clear Sropt Lropt
V=sparse(blkdiag(V_Blk{:}));
end

function y=RSRproduct(D,A,x)
%y=D*A*D'
y=D'*x;
y=A*y;
y=D*y;
end