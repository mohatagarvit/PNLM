%% Denoising demo using Pruned Non-Local Means (PNLM)
%  Comparing the results of NLM and PNLM.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%  img      : clean grayscale image
%  h        : width of Gaussian
%  P        : half-size of patch 
%  S        : half-search window 
%
%  Author   : S. Ghosh and Kunal N. Chaudhury, Indian Institute of Science.
%  Date     : Feb. 17, 2017
%
%  Reference: 
%  S. Ghosh and K. N. Chaudhury, "Pruned Non-Local Means", IET Image Processing.  
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; close all force; clear all;

% Clean image
% img = double(imread('lena512.png'));
% img = double(imread('barbara512.png'));
% img = double(imread('cameraman256.png'));
img    = double(imread('peppers256.png'));
[m, n] = size(img);
peak  = 255;

% Add noise
sigma     =  20;
imgNoisy  =  img  +  sigma * randn(m,n);

% NLM parameters
S  = 10;            % half-search window
K  = 3;             % half-size of patch
h  = 10 * sigma;    % width of Gaussian

% Additional PNLM parameters
param.alpha = 200;
param.epsilon = 10^-3;
param.itr_max = 50;
param.R = double((sqrt(5)-1)/2);

%% Denoising
% NLM
disp('Non-Local Means is running.......');
tic;
[imgOut_NLM,  W_nlm, M_nsy, Y_swp] = nlm(imgNoisy, S, K, h);
t_NLM = toc;

% PNLM
disp('Pruned NLM is running.......');
l_est = 4.3 *10^-7 * sigma^3 - 1.1 * 10^-4 *sigma^2 + 9.2 * 10^-3 *sigma^1 + 0.039
l_min = max(0, (l_est - 0.10));
l_max = l_est + 0.10;
% tic;
[imgOut_PNLM, MSE_sure1] = pnlm(W_nlm, imgNoisy, S, K, h, sigma, M_nsy, Y_swp,l_min, l_max, param) ;
t_PNLM = toc;

%% Results
PSNR_Nsy = 10 * log10(m * n * peak^2 / sum(sum((imgNoisy - img).^2)) )
PSNR_NLM = 10 * log10(m * n * peak^2 / sum(sum((imgOut_NLM - img).^2)) )
PSNR_PNLM = 10 * log10(m * n * peak^2 / sum(sum((imgOut_PNLM - img).^2)) )

SSIM_Nsy = 100*ssim(img, imgNoisy)
SSIM_NLM = 100*ssim(img, imgOut_NLM)
SSIM_PNLM = 100*ssim(img, imgOut_PNLM)

t_NLM
t_PNLM

%% Displaying images
figure('Units','normalized','Position',[0 0 1 1]);
colormap gray,
subplot(2,2,1); imshow(uint8(img)); title('Clean');
subplot(2,2,2); imshow(uint8(imgNoisy)); 
title([ 'Noisy, ', num2str(PSNR_Nsy, '%.2f'), 'dB, ', num2str(SSIM_Nsy, '%.2f')] , 'FontSize', 10),
subplot(2,2,3); imshow(uint8(imgOut_NLM)); 
title([ 'NLM, ', num2str(PSNR_NLM, '%.2f'), 'dB, ', num2str(SSIM_NLM, '%.2f')] , 'FontSize', 10),
subplot(2,2,4); imshow(uint8(imgOut_PNLM)); 
title([ 'PNLM, ', num2str(PSNR_PNLM, '%.2f'), 'dB, ', num2str(SSIM_PNLM, '%.2f')] , 'FontSize', 10),
