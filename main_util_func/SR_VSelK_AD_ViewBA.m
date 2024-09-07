function [A_Blk,Sig_Blk,US_Blk,Ker_Blk,D_Blk,Cold_prod,Half_prod,D_block,ord_Blk]=SR_VSelK_AD_ViewBA(X,La_Prel,Beta_old,Alpha_old,Pas,para)
%========================================================%
%X£ºcell(1,Nview) 
%La_Prel: initial label
%Beta_old: Nview*1
%Alpha_old£ºNview*1
%Pas, para: parameters
%SR, select a "representative" kernel matrix for each view 
%========================================================%
global KernelTypes KernelPostProcessTypes Degrees Nbaker
A_Blk=cell(1,Pas.Nview);
Ker_Blk=cell(1,Pas.Nview);
D_Blk=cell(1,Pas.Nclu);
Sig_Blk=cell(1,Pas.Nview);
US_Blk=cell(1,Pas.Nview);
ord_Blk=zeros(1,Pas.Nclu);
Cold_prod=cell(1,Pas.Nclu);
Half_prod=cell(1,Pas.Nclu);
sub_kernel=cell(1,Pas.Nview);
ind=cell(1,Pas.Nclu);
%the size and corresponding position of each block are the same for each view
for j=1:Pas.Nclu 
    [ind{j},~]=find(La_Prel==j);
    ord_Blk(j)=numel(ind{j}); %order number of this block
end
for v=1:Pas.Nview
   for j=1:Pas.Nclu 
       view_subfea=X{v}(ind{j},:);
       %K_j^(V) is calculated in four cases
       if ord_Blk(j)==1 
          A_Blk{v}{j}=1;
          for ii=1:Nbaker
             if Beta_old{v}(ii)>1e-3
                 Ker_Blk{v}{j}{ii}=1;
             end
          end
       elseif ord_Blk(j)>1 && ord_Blk(j)<=1000  
          sub_ker=zeros(ord_Blk(j),ord_Blk(j));
          %form base kernel matrix according to beta^(v)
          for ii=1:Nbaker
             if Beta_old{v}(ii)>1e-3
                kernel_option = [];
               switch lower(KernelTypes{ii})
                  case lower('Linear')
                       kernel_option.KernelType = 'Linear';
                       for iPost = KernelPostProcessTypes
                           sub_kernel = constructKernel(view_subfea,[], kernel_option);
                           sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                       end
                  case lower('Polynomial')
                       kernel_option.KernelType = 'Polynomial';
                       kernel_option.d = Degrees(ii);
                       for iPost = KernelPostProcessTypes
                          sub_kernel = constructKernel(view_subfea,[], kernel_option); 
                          sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                       end
                  case lower('PolyPlus')
                      kernel_option.KernelType = 'PolyPlus';
                      kernel_option.d = Degrees(ii);
                      for iPost = KernelPostProcessTypes
                          sub_kernel = constructKernel(view_subfea,[], kernel_option);
                          sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                      end
                  case lower('Gaussian')
                      kernel_option.KernelType = 'Gaussian';
                      kernel_option.t = Degrees(ii);
                      for iPost = KernelPostProcessTypes
                          sub_kernel = constructKernel(view_subfea,[], kernel_option); 
                          sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                      end
                  otherwise
                   error('KernelType does not exist!');
               end
               sub_ker=sub_ker+Beta_old{v}(ii)*sub_kernel; 
               Ker_Blk{v}{j}{ii}=sub_kernel;
               clear sub_kernel 
             end
          end
          A_Blk{v}{j} = sub_ker;  %ord_block(j)*ord_block(j)
          clear view_subfea sub_ker kernel_option
      elseif ord_Blk(j)>1000 && ord_Blk(j)<=50000  
          %----------------------rank ropt approximation
          Nsamker=ceil(para.SmSker*ord_Blk(j));
          cols = randperm(ord_Blk(j),Nsamker); 
          Rview_fea=view_subfea(cols,:);
          sub_ker=zeros(ord_Blk(j),Nsamker);
          for ii=1:Nbaker
              if Beta_old{v}(ii)>1e-3
                 kernel_option = [];
                 switch lower(KernelTypes{ii})
                     case lower('Linear')
                          kernel_option.KernelType = 'Linear';
                          for iPost = KernelPostProcessTypes
                             sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                             sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                          end
                     case lower('Polynomial')
                          kernel_option.KernelType = 'Polynomial';
                          kernel_option.d = Degrees(ii);
                          for iPost = KernelPostProcessTypes
                             sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                             sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                          end
                     case lower('PolyPlus')
                          kernel_option.KernelType = 'PolyPlus';
                          kernel_option.d = Degrees(ii);
                          for iPost = KernelPostProcessTypes
                             sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                             sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                          end
                     case lower('Gaussian')
                          kernel_option.KernelType = 'Gaussian';
                          kernel_option.t = Degrees(ii);
                          for iPost = KernelPostProcessTypes
                             sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option); 
                             sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                          end
                     otherwise
                     error('KernelType does not exist!');
                 end
                 sub_ker=sub_ker+Beta_old{v}(ii)*sub_kernel; %C after combination
                 Ker_Blk{v}{j}{ii}=sub_kernel;
                 clear sub_kernel
              end
          end
          ropt=Nsamker-para.Ovsam;   
          W=sub_ker(cols,:);  
          [Q_ker,R_ker] = qr(sub_ker,0); 
          W=R_ker*pinv(W,1e-12)*R_ker'; 
          W=max(W,W');  
          [P_W,Sigma,~]=svd(W,0);
          P_W=P_W(:,1:ropt);
          Sig_Blk{v}{j}=Sigma(1:ropt,1:ropt);
          US_Blk{v}{j}=Q_ker*P_W;
          clear sub_ker view_subfea Rview_fea samker_num cols ker ropt W Q_ker R_ker P_W Sigma kernel_option
       else
          %---------------------rank ropt approximation
          Nsamker=ceil(para.LaSker*ord_Blk(j));
          cols = randperm(ord_Blk(j),Nsamker); 
          Rview_fea=view_subfea(cols,:);
          sub_ker=zeros(ord_Blk(j),Nsamker);
          for ii=1:Nbaker
              if Beta_old{v}(ii)>1e-3
                 kernel_option = [];
                 switch lower(KernelTypes{ii})
                     case lower('Linear')
                         kernel_option.KernelType = 'Linear';
                         for iPost = KernelPostProcessTypes
                             sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                             sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                         end
                     case lower('Polynomial')
                         kernel_option.KernelType = 'Polynomial';
                         kernel_option.d = Degrees(ii);
                         for iPost = KernelPostProcessTypes
                             sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option); 
                             sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                         end
                    case lower('PolyPlus')
                        kernel_option.KernelType = 'PolyPlus';
                        kernel_option.d = Degrees(ii);
                        for iPost = KernelPostProcessTypes
                            sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option);
                            sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                        end
                   case lower('Gaussian')
                        kernel_option.KernelType = 'Gaussian';
                        kernel_option.t = Degrees(ii);
                        for iPost = KernelPostProcessTypes
                            sub_kernel = constructKernel(view_subfea,Rview_fea, kernel_option); 
                            sub_kernel = KernelNormalize(sub_kernel, iPost{1});%#ok
                        end
                     otherwise
                     error('KernelType does not exist!');
                 end
                 sub_ker=sub_ker+Beta_old{v}(ii)*sub_kernel; %C after combination
                 Ker_Blk{v}{j}{ii}=sub_kernel;
                 clear sub_kernel
              end
          end
          ropt=Nsamker-para.Ovsam;   
          W=sub_ker(cols,:);  
          [Q_ker,R_ker] = qr(sub_ker,0); 
          W=R_ker*pinv(W,1e-12)*R_ker'; 
          W=max(W,W');  
          [P_W,Sigma,~]=svd(W,0);
          P_W=P_W(:,1:ropt);
          Sig_Blk{v}{j}=Sigma(1:ropt,1:ropt);
          US_Blk{v}{j}=Q_ker*P_W;
          clear sub_ker view_subfea Rview_fea samker_num cols ker ropt W Q_ker R_ker P_W Sigma kernel_option
       end    
   end 
end
clear ind
Dbt=cell(1,Pas.Nview);
for j=1:Pas.Nclu
    %Save U,S for each block
    if ord_Blk(j)==1
        Dbt{j}=1;
        D_Blk{j}=1;
        Cold_prod{j}=1;
    elseif ord_Blk(j)>1 && ord_Blk(j)<=1000 %combination
          Abent_block=zeros(ord_Blk(j));
          for v=1:Pas.Nview
             Abent_block=Abent_block+Alpha_old(v)*A_Blk{v}{j};
          end
             Dbt{j}=sum(Abent_block,2);
          for rd=1:ord_Blk(j)
             if Dbt{j}(rd)<=1e-12
                Dbt{j}(rd)=1;
             end
          end
          r=1./sqrt(Dbt{j});
          D_Blk{j}=r;
          Cold_prod{j}=r.*(r'.*Abent_block);
          Cold_prod{j}=max(Cold_prod{j},Cold_prod{j}');
          clear Abent_block Dbent_block r
    else
          Ubt=[];
          Dbt{j}=zeros(ord_Blk(j),1);
          SiGb=cell(1,Pas.Nview);
          for v=1:Pas.Nview
             Ubt=[Ubt US_Blk{v}{j}];
             SiGb{v}=Alpha_old(v)*Sig_Blk{v}{j};
             Dbt{j}=Dbt{j}+abs(US_Blk{v}{j}*(SiGb{v}*(US_Blk{v}{j}'*ones(ord_Blk(j),1))));
          end  
          SiGbt=sparse(blkdiag(SiGb{:}));
          for rd=1:ord_Blk(j)
            if Dbt{j}(rd)<=1e-12
               Dbt{j}(rd)=1;
            end
          end
          clear rd SiGb
          U=1./sqrt(Dbt{j}); 
          D_Blk{j}=U;
          U=U.*Ubt; 
          DIAG_Sig=sqrt(diag(SiGbt));
          Half_prod{j}=DIAG_Sig'.*U; %D{j}^(-1/2)*U{j}*SiGbt{j}^(1/2) 
          clear Ubt SiGbt U DIAG_Sig
     end
end
%D_block=sum(sparse(blkdiag(Dbt{:})),2);  %Nsam*1
D_block=cat(1,Dbt{:});
clear Dbt 
end