function f=NLmeansGlb(J,I,patchSigma,patchSize,filtSigma)

%% SCRIPT: SAMPLE_KERNEL
%
% Sample usage of GPU kernel through MATLAB

% DEPENDENCIES
%
%  nlmGlobal.cu
%

 % clear variables;

  %% PARAMETERS 
  % filter sigma value
  %filtSigma = 0.02;
  %patchSize = 5;
  %patchSigma = 5/3;
  %patchsize =[5 5];
 
  
  threadsPerBlock = [32 32];
  [m,n]=size(J);
  


  %% (BEGIN)

 % fprintf('...begin %s...\n',mfilename);  
  %% gaussian patch
  H = fspecial('gaussian',patchSize, patchSigma);
  H = H(:) ./ max(H(:));
 
  H=reshape(H,patchSize);
  %% KERNEL
  
%  k = parallel.gpu.CUDAKernel( '../cuda/sampleAddKernel.ptx', ...
 %                            '../cuda/sampleAddKernel.cu','meanCalc');
  p =  parallel.gpu.CUDAKernel( '../cuda/nlmGlobal.ptx', ...
                                    '../cuda/nlmGlobal.cu','Zcalc');
      
  fi=  parallel.gpu.CUDAKernel( '../cuda/nlmGlobal.ptx', ...
                              '../cuda/nlmGlobal.cu','fCalc');
                          
  numberOfBlocks  = ceil( [m n] ./ threadsPerBlock );
  

  p.ThreadBlockSize = threadsPerBlock;
  p.GridSize        = numberOfBlocks;
  p.SharedMemorySize= 16000;
  fi.ThreadBlockSize = threadsPerBlock;
  fi.GridSize        = numberOfBlocks;
  fi.SharedMemorySize= 16000;
  %% DATA
  
  I=J;
  J=padarray(J,(patchSize-1)./2,'symmetric','both');
  [m ,n]=size(J);
  J=gpuArray(single(J));
 
  Z= zeros([m n], 'gpuArray');
  Z=single(Z);
  f= zeros([m n],'gpuArray');
  f=single(f);
  H=gpuArray(single(H));
fileID=fopen('NLMGLOBAL','a+');  
  
  tic;

  Z=gather(feval(p,J,Z,H,patchSize(1),patchSigma,filtSigma,m,n));
time1=toc;
tic; 
  f=gather(feval(fi,J,Z,H,f,patchSize(1),patchSigma,filtSigma,m,n));
time2=toc;
  f=f(1+(patchSize(1)-1)/2:end-(patchSize(1)-1)/2,1+(patchSize(2)-1)/2:end-(patchSize(1)-1)/2);
I=single(I);
fprintf(fileID,'\n\nGlobal memory Image: %dx%d patchSize: %dx%d \n',size(I),patchSize);
fprintf(fileID,'Z kernel : %f ,  f kernel: %f  \n',time1,time2);
[peakpsnr,snr]=psnr(f,I,1);
  fprintf(fileID,'\n The peak -SNR value is %f',peakpsnr);
  fprintf(fileID,'\n The SNR value is %f \n',snr);
end
  
  %% SANITY CHECK
  
 % fprintf('Error: %e\n', norm( B - (A+1), 'fro' ) );
  
  %% (END)

  %fprintf('...end %s...\n',mfilename);

%%%------------------------------------------------------------
%
% AUTHORS
%
%   Panagiotis Kasparidis                 pankasgeo@e.ece.auth.gr
%
% ------------------------------------------------------------
