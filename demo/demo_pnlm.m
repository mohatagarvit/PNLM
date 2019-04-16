
%% Pruned Non-Local Means (PNLM)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  INPUTS:
%  x     : grayscale input image .
%  S     : search window parameter. Size of search window: [(2S+1), (2S+1)].
%  K     : patch width parameter. Size of each patch: [(2K+1), (2K+1)].
%  Sigma : noise level of the noisy image; typically [0 - 100].
%
%  h     : 10 * sigma, smoothing parameter in NLM.
%  alpha : slope of the sigmoidal function, used for smoothing.
% 
%  
%  CALLED FUNCTION:
%  GUI_mex.c: mex file that performs the pnlm computation.
%
%  Authors     : Garvit Mohata, Sanjay Ghosh, and Kunal Narayan Chaudhury.
%
%  Reference   : S. Ghosh, A. K. Mandal, and K. N. Chaudhury, "Pruned Non-Local Means", 
%               IET Image Processing, vol. 11, no. 5, pp. 317-323, April 2017.
%
%  Date: June 28, 2017
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear; clc; close all force;



x = double(imread('./images/peppers.png')); 
[m, n] =  size(x);

S      =  10;
K      =  3;
sigma  =  50;
h      =  10 * sigma;
alpha  =  100;


% Adding Gaussian Noise
y      =  x + sigma .* randn(m,n);

y_padded = padarray(y,[S+K S+K],'symmetric');

% PNLM filtering
tic;  
[x_hat, diff_x_hat_with_y] = GUI_mex (y_padded,sigma,S,K,alpha,h);
toc;                

% Results and display
peak       =  255;
PSNR_pnlm  =  10 * log10(m * n * peak^2 / sum(sum((x_hat - x).^2)) );
PSNR_noisy =  10 * log10(m * n * peak^2 / sum(sum((y - x).^2)) );
SSIM       =  ssim(x,x_hat);
str_noisy  =  sprintf('Noisy  (PSNR = %.2f dB)',PSNR_noisy);
str_pnlm   =  sprintf('PNLM (PSNR = %.2f dB)',PSNR_pnlm);
figure;
subplot(2,2,1); imshow(uint8(x),[0 255]);       title('Clean Image','fontsize',8,'fontweight','bold'), axis('image','off');
subplot(2,2,2); imshow(uint8(y),[0 255]);       title(str_noisy,'fontsize',8,'fontweight','bold'), axis('image','off');
subplot(2,2,3); imshow(uint8(x_hat),[0 255]);   title(str_pnlm,'fontsize',8,'fontweight','bold'), axis('image','off');
subplot(2,2,4); imshow(uint8(y-x_hat),[0 255]); title('Method Noise','fontsize',8,'fontweight','bold'), axis('image','off');
