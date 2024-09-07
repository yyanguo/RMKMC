function [Ufac,vecSig]=Ker_svd(kerM,NKMsamker)
%=================================================%
%kerM: s*s small kernel matrix (explicitly formed)
%Ufac: s*? left singular vector of kerM   
%request£ºNsam>>s
%=================================================%
[Ufac,~,~]=svd(kerM,0);
% DSig=diag(Sig);
% DSig=DSig/DSig(1); %relative truncation
% vecSig=length(DSig(DSig>1e-12)); %truncation svd
vecSig=ceil(0.1*NKMsamker); 
%select 0.1*s columns 
Ufac=Ufac(:,1:vecSig);
clear DSig 
end