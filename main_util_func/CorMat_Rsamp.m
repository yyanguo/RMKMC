function Corker = CorMat_Rsamp(Dfea,Pas)
%=============================================%
%for data Nsam*d£¬sample s rows randomly£¬s*d to replace Nsam*d, svd, cos
%Dfea£ºdataset
%Corker£ºkernel similarity matrix
%=============================================%
if Pas.Nsam<=5000
    sam_prob=0.1;
elseif Pas.Nsam>5000 && Pas.Nsam<=20000
    sam_prob=0.05;
elseif Pas.Nsam>20000 && Pas.Nsam<=50000
    sam_prob=0.01;
else
    sam_prob=0.005;
end
NKMsamker=ceil(sam_prob*Pas.Nsam);  %s
cols = randperm(Pas.Nsam,NKMsamker); 
Cent=Dfea(cols,:); %s*d
Corker=CalMat_Corker(Cent,NKMsamker);
clear Dfea cols Cent NKMsamker sam_prob
end