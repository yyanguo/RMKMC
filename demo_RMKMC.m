%%===============================================================%%
% Randomly initial label + two restarted algorithms 
% Overall framework: 
%        1. KM + update alpha + update view beta
%        2. SR + update alpha + update view beta
% Written by Guo yongyan
%%===============================================================%%
clc; clear; close all; warning off
trials=2;   %mean of 10 times
para.LaSker=0.01; para.SmSker=0.05;   para.Ovsam=5;
Lam=0.8; Niu=1e-2;  Gam=2^-1;
global KernelTypes KernelPostProcessTypes Degrees Nbaker
KernelTypes = {'Linear', 'PolyPlus', 'PolyPlus', 'PolyPlus','Polynomial','Polynomial','Polynomial',...
    'Gaussian','Gaussian', 'Gaussian', 'Gaussian'};
KernelPostProcessTypes = {'MAX'};
Degrees = [1, 2, 4, 6, 2, 4, 6, 0.01, 0.1, 1, 10]; % 11 base kernels
Nbaker=length(Degrees);
file_path = '.\Data\';
files = dir(strcat(file_path,'\','*.mat'));
datanum = length(files);
mkdir('.\Res_RMKMC\');
Begin_Time = ['Running at *** ' datestr(now)]
for r=1:datanum
   dataName=[files(r).name]
   load([file_path,'\',dataName],'X','Y','Pas');
   Dimsion=[num2str(Pas.Nsam) '*' num2str(Pas.Nclu) '*' num2str(Pas.Nview)]
   addpath('.\main_util_func\');
   addpath('.\clustering_measure\');
   Tprel=zeros(1,trials);
   for tr = 1:trials
      fprintf('\n\t TRIAL %d\n', tr);
      tp=tic;
      La_Prel=crossvalind('Kfold',Pas.Nsam,Pas.Nclu);
      if length(unique(La_Prel))~=Pas.Nclu
         error('incorrect cluster')
      end
      Tprel(tr)=toc(tp);
      %=================Alg.1===================%
      [KV.Obj{tr},KV.SaveAlpha{tr},KV.res(tr,:),KV.time(tr),KV.iter(tr),KV.SaveBeta{tr}]=main_ReKM_SRK_ViewBA(X,Y,@CorMat_Rsamp,Pas,La_Prel,para,Niu,Gam);
      disp(['KV--',num2str(tr),'--Over']);
      %=================Alg.2====================%
      [RV.Obj{tr},RV.SaveAlpha{tr},RV.res(tr,:),RV.time(tr),RV.iter(tr),RV.SaveBeta{tr}]=main_ReSR_SRK_ViewBA(X,Y,@CorMat_Rsamp,Pas,La_Prel,para,Lam,Niu,Gam);
      disp(['RV--',num2str(tr),'--Over']);
   end
   rmpath('.\main_util_func\');
   rmpath('.\clustering_measure\');
   %=================Alg.1=====================%
   KV.mean_Res=mean(KV.res,1);   KV.mean_Time=mean(KV.time)+mean(Tprel);  KV.RESULT=[KV.mean_Res KV.mean_Time];
   %=================Alg.2=====================%
   RV.mean_Res=mean(RV.res,1);   RV.mean_Time=mean(RV.time)+mean(Tprel);  RV.RESULT=[RV.mean_Res  RV.mean_Time];
   
   OUTPUT_inx=[KV.RESULT;RV.RESULT]; 
   %======save data========%
   save(['.\Res_RMKMC\',dataName],'KV','RV','OUTPUT_inx');
   clear X Y Pas La_Prel KV RV OUTPUT_inx 
   disp('========Results have been saved==========') 
end
clear global KernelTypes KernelPostProcessTypes Degrees Nbaker para Lam Niu Gam
End_Time = ['Running over *** ' datestr(now)]
