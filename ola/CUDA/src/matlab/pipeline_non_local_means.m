    %% SCRIPT: PIPELINE_NON_LOCAL_MEANS
%
% Pipeline for non local means algorithm as described in [1].
%-
% The code thus far is implemented in CPU.
%
% DEPENDENCIES
%
% [1] Antoni Buades, Bartomeu Coll, and J-M Morel. A non-local
%     algorithm for image denoising. In 2005 IEEE Computer Society
%     Conferen ce on Computer Vision and Pattern Recognition (CVPR’05),
%      volume 2, pages 60–65. IEEE, 2005.
%
  
  clear all %#ok
  close all

  %% PARAMETERS
  
  % input image
  pathImg   = '../data/lena128.mat';
  strImgVar = 'lena128';
  
  % noise
  noiseParams = {'gaussian', ...
                 0,...
                 0.001};
  
  % filter sigma value
  filtSigma = 0.02;
  patchSize = [7 7];
  patchSigma = 5/3;
  
  %% USEFUL FUNCTIONS

  % image normalizer
  normImg = @(I) (I - min(I(:))) ./ max(I(:) - min(I(:)));
  
  %% (BEGIN)

  fprintf('...begin %s...\n',mfilename);  
  
  %% INPUT DATA
  
  fprintf('...loading input data...\n')
  
  ioImg = matfile( pathImg );
  I     = ioImg.(strImgVar);
  
  %% PREPROCESS
  
  fprintf(' - normalizing image...\n')
  I = normImg( I );
  
   figure('Name','Original Image');
   imagesc(I); axis image;
   colormap gray;
  
  %% NOISE
  
  fprintf(' - applying noise...\n')
  J = imnoise( I, noiseParams{:} );
%   figure('Name','Noisy-Input Image');
%   imagesc(J); axis image;
%   colormap gray;
  %% NON LOCAL MEANS
  
%  tic;
%  If = nonLocalMeans( J,I, patchSize, filtSigma, patchSigma );
%  toc
%  [peakpsnr,snr]=psnr(If,I,1);
%  fprintf('\n The peak -SNR value is %0.04f',peakpsnr);
%  fprintf('\n The SN value is %f \n',snr);

    %% VISUALIZE RESULT
  
%   figure('Name', 'Filtered image');
%   imagesc(If); axis image;
%   colormap gray;
  %% NON LOCAL MEANS CUDA

  f= NLmeansGlb(J,I,patchSigma,patchSize,filtSigma);

  imwrite(f,'../outputs/Global/lena128(7x7).jpg');
  %% VISUALIZE RESULT
  
%   figure('Name', 'My Filtered image');
%   imagesc(f); axis image;
%   colormap gray;

  
  %% VISUALIZE RESULT
  
% figure('Name', 'My Filtered image');
%  imagesc(Sf); axis image;
%  colormap gray;
  
   %% NON LOCAL MEANS CUDA SHARED
fprintf('... Cuda lastOne.cu times : \n'); 
% tic;
Sf= SharedWithMatlab(J,I,patchSize,filtSigma,patchSigma);
%toc;
imwrite(Sf,'../outputs/lastOne/lena128(7x7).jpg');
%tic;
fprintf('... Cuda Shared.cu times: \n');
S2f=SharedKernel(J,I,patchSigma,patchSize,filtSigma);
%toc;
imwrite(S2f,'../outputs/Shared/lena128(7x7).jpg');
  %% VISUALIZE RESULT
%
   figure('Name', 'My Filtered image');
   imagesc(S2f); axis image;
   colormap gray;
% 
 % figure('Name', 'Residual');
 % imagesc(S2f-J); axis image;
 % colormap gray;
  
  %% (END)

  fprintf('...end %s...\n',mfilename);


%%------------------------------------------------------------
%
% AUTHORS
%
%   Dimitris Floros                         fcdimitr@auth.gr
%
% VERSION
%
%   0.1 - December 28, 2016
%
% CHANGELOG
%
%   0.1 (Dec 28, 2016) - Dimitris
%       * initial implementation
%
% ------------------------------------------------------------
