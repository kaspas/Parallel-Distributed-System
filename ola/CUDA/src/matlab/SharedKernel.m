function f=SharedKernel(J,I,patchSigma,patchSize,filtSigma)

%% SCRIPT: SharedKernel
%
% Sample usage of GPU kernel through MATLAB

% DEPENDENCIES
%
%  Shared.cu
%

  
  [m,n]=size(J);
  


  %% (BEGIN)

  %% gaussian patch
  H = fspecial('gaussian',patchSize, patchSigma);
  H = H(:) ./ max(H(:));
  H=gpuArray(H);
  H=single(H);
  %% KERNEL
  threadsPerBlock=[16 16];  
  p =  parallel.gpu.CUDAKernel( '../cuda/Shared.ptx', ...
                                    '../cuda/Shared.cu','Zcalc');
      
  fi=  parallel.gpu.CUDAKernel( '../cuda/Shared.ptx', ...
                              '../cuda/Shared.cu','fCalc');
                          
  numberOfBlocks  = ceil( [m n] ./ threadsPerBlock );
  
  p.ThreadBlockSize = threadsPerBlock;
  p.GridSize        = numberOfBlocks;
  p.SharedMemorySize= 16000;
  setConstantMemory(p,'s_H',H);
  fi.ThreadBlockSize = threadsPerBlock;
  setConstantMemory(fi,'s_H',H);
  fi.GridSize        = numberOfBlocks;
  fi.SharedMemorySize= 16000;
  %% DATA
  
  %A = rand([m n], 'gpuArray');
  %B = zeros([m n], 'gpuArray');
  %[w,h]=size(J);
I=single(I);
  J=padarray(J,(patchSize-1)./2,'symmetric','both');
  [m ,n]=size(J);
  J=gpuArray(single(J));
 
  Z= zeros([m n], 'gpuArray');
  Z=single(Z);
  f= zeros([m n],'gpuArray');
  f=single(f);
fileID=fopen('SharedKernelNumerics','a+');  
tic;
  Z= gather(feval(p,J,Z,patchSize(1),filtSigma,m,n));  
time1=toc;
tic;
  f=gather(feval(fi,J,Z,f,patchSize(1),filtSigma,m,n));
time2=toc;
fprintf(fileID,'\n\nShared Image: %dx%d patchSize: %dx%d \n',size(I),patchSize);
   f=f(1+(patchSize(1)-1)/2:end-(patchSize(1)-1)/2,1+(patchSize(2)-1)/2:end-(patchSize(2)-1)/2);
[peakpsnr,snr]=psnr(f,I,1);
fprintf(fileID,'Z kernel : %f ,  f kernel: %f  \n',time1,time2);
  fprintf(fileID,'\n The peak -SNR value is %f',peakpsnr);
  fprintf(fileID,'\n The SNR value is %f \n',snr);
fclose(fileID);


  end



%%%------------------------------------------------------------
%
% AUTHORS
%
%   Panagiotis Kasparidis                         pankasgeo@ece.auth.gr
%
%
% ------------------------------------------------------------
