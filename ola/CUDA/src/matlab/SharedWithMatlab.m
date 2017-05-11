function f = SharedWithMatlab(I,J, patchSize, filtSigma, patchSigma)
%% SCRIPT: SAMPLE_KERNEL
%
% Sample usage of GPU kernel through MATLAB
%
% DEPENDENCIES
%
%  SharedWithMatlab.cu
%
  
  %clear variables;

  %% PARAMETERS
  [w,h]=size(I);
  threadsPerBlock = [32 32];
 

  %% (BEGIN)

    %% GUASSIAN PATCH
    H = fspecial('gaussian',patchSize, patchSigma);
    
    H = H(:) ./ max(H(:));
    H=gpuArray(H);
    H=single(H);
    
  %% KERNEL
  
    p = parallel.gpu.CUDAKernel( '../cuda/SharedWithMatlab.ptx', ...
                               '../cuda/SharedWithMatlab.cu','Zcalc');
    fi = parallel.gpu.CUDAKernel( '../cuda/SharedWithMatlab.ptx', ...
                               '../cuda/SharedWithMatlab.cu','fCalc');

  numberOfBlocks  = ceil( [w h] ./ threadsPerBlock );
  p.SharedMemorySize = 16000;
  p.ThreadBlockSize = threadsPerBlock;
  p.GridSize        = numberOfBlocks;
  setConstantMemory(p,'s_H',H);
    fi.SharedMemorySize = 16000;
  fi.ThreadBlockSize = threadsPerBlock;
  fi.GridSize        = numberOfBlocks;
  setConstantMemory(fi,'s_H',H);

  %% DATA
  I=padarray(I,(patchSize-1)/2,'symmetric','both');

    I=single(I);
    A = gpuArray(I);
    [m,n] =size(I);

  Z = zeros([m n], 'gpuArray');
  f = zeros([m n], 'gpuArray');
  Z=single(Z);
  f=single(f);

  H=single(H);
  fileID=fopen('SWM','a+');  

tic;
  for j=[1 threadsPerBlock(1):threadsPerBlock(1):w-threadsPerBlock(1)]
	
     for i= [1 threadsPerBlock(2):threadsPerBlock(2):h-threadsPerBlock(2)]
        dim_x=i:i+threadsPerBlock(1)+patchSize(1)-1;
        dim_y=j:j+threadsPerBlock(2)+patchSize(2)-1;
	
         Z = (feval(p,A,A(dim_x,dim_y),Z,H,patchSize(1,1),filtSigma,m,n));
	wait(gpuDevice);

     end
  end
  Z=gather(Z);
time1=toc;
%----------------------------------------------------------------------------------------------
tic;
  for j=[1 threadsPerBlock(1):threadsPerBlock(1):w-threadsPerBlock(1)]
	
     for i= [1 threadsPerBlock(2):threadsPerBlock(2):h-threadsPerBlock(2)]
        dim_x=i:i+threadsPerBlock(1)+patchSize(1)-1;
        dim_y=j:j+threadsPerBlock(2)+patchSize(2)-1;
	
         f = (feval(fi,A,A(dim_x,dim_y),Z,H,f,patchSize(1,1),filtSigma, m, n));
	wait(gpuDevice);
	
     end
  end
  f=gather(f);
time2=toc;  

f=f(1+(patchSize(1)-1)/2:end-(patchSize(1)-1)/2,1+(patchSize(2)-1)/2:end-(patchSize(2)-1)/2);

 J=single(J);
fprintf(fileID,'\n\nShared With Matlab Image: %dx%d patchSize: %dx%d \n',size(I),patchSize);

[peakpsnr,snr]=psnr(f,J,1);
fprintf(fileID,'Z kernel : %f ,  f kernel: %f  \n',time1,time2);
  fprintf(fileID,'\n The peak -SNR value is %f',peakpsnr);
  fprintf(fileID,'\n The SNR value is %f \n',snr);


end

%%%------------------------------------------------------------
%
% AUTHORS
%
%   Panagiotis Kasparidis                         pankasgeo@ece.auth.gr
%
%
% ------------------------------------------------------------
