function Corker=CalMat_Corker(Cent,NKMsamker)
%=====================================================%
%for X{i}, perform svd to get "Corker"
%Cent£ºs*d 
%Corker£ºNbaker*Nbaker 
%==================================================%
global KernelTypes KernelPostProcessTypes Degrees Nbaker
%--------------------compute sub_kernel: s*s*Nbaker
Ufac=cell(1,Nbaker); 
for i =1:Nbaker  
     kernel_option = [];
     switch lower(KernelTypes{i})
        case lower('Linear')
            kernel_option.KernelType = 'Linear';
            for iPost = KernelPostProcessTypes
                    Ker_sub = constructKernel(Cent,[], kernel_option);
                    Ker_sub = KernelNormalize(Ker_sub, iPost{1});%#ok
                    [Ufac{i},~]=Ker_svd(Ker_sub,NKMsamker);
            end
        case lower('Polynomial')
            kernel_option.KernelType = 'Polynomial';
                kernel_option.d = Degrees(i);
                for iPost = KernelPostProcessTypes
                    Ker_sub = constructKernel(Cent,[], kernel_option);
                    Ker_sub = KernelNormalize(Ker_sub, iPost{1});%#ok
                    [Ufac{i},~]=Ker_svd(Ker_sub,NKMsamker);
                end
        case lower('PolyPlus')
            kernel_option.KernelType = 'PolyPlus';
                kernel_option.d = Degrees(i);
                for iPost = KernelPostProcessTypes
                    Ker_sub = constructKernel(Cent,[], kernel_option);
                    Ker_sub = KernelNormalize(Ker_sub, iPost{1});%#ok
                    [Ufac{i},~]=Ker_svd(Ker_sub,NKMsamker); 
                end
        case lower('Gaussian')
            kernel_option.KernelType = 'Gaussian';
                kernel_option.t = Degrees(i);
                for iPost = KernelPostProcessTypes
                    Ker_sub = constructKernel(Cent,[], kernel_option);
                    Ker_sub = KernelNormalize(Ker_sub, iPost{1});%#ok
                    [Ufac{i},~]=Ker_svd(Ker_sub,NKMsamker);
                end
        otherwise
            error('KernelType does not exist!');
    end
end 
clear Ker_sub kernel_option
%--------------------compute correlevant matrix
Corker = zeros(Nbaker);
for ii=1:(Nbaker-1)
    for jj=(ii+1):Nbaker
        Corker(ii,jj) = norm(Ufac{jj}'*Ufac{ii});
    end
end 
Corker = (Corker+Corker')+eye(Nbaker);
clear Ufac Cent 
end