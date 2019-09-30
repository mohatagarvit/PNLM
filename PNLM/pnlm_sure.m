function [fNLM, SURE] = pnlm_sure(A, f, S, K, h, phi_m, alpha, sig, M_nsy, Y_swp)

%--------------------------------------------------------------------------
% nlm_sure : computes truncated NLM image and SURE of the MSE;
%
% ------------------ OUTPUT ----------------------------------------
% fNLM     : truncated and smoothed-NLM-denoised-image 
% SURE     : Stein's Unbiased Risk Estimate for the Mean Square Error
%
%------------------ INPUT -----------------------------------------
% A        : NLM actual weight matrix
% f        : Noisy image
% S        : Search window
% K        : Patch size
% h        : Smoothing parameter in NLM
% phi_m    : Truncation coefficient or mean of the sigmf function 
% alpha    : Slope of the sigmf function 
% sigma    : Square-root of the additive gaussian noise variance
% M_nsy    : Matrix, which stores noisy image arranged to multiply with corresponding NLM weights. Size = ( N^2, (2*S+1)^2 ).
% Y_swp    : Matrix, which stores the noisy image pixel arranged in a way required for SURE MSE computation. Size = ( N^2, (2*S+1)^2 ).
%
%--------------------------------------------------------------------------

[M,N] = size(f);

E = ones(N^2, (2*S+1)^2); % E: All 1 matrix of size ((M N), (2*S+1)^2)
ONE = ones((2*S+1)^2, 1); % ONE: All 1 column vector of size (2*S+1)^2.
one = ones(N^2, 1);       % one: All 1 column vector of size N^2.

W_out = A./( E + exp( - alpha.*(A - phi_m*E)) ); %W_out: truncated weight matrix

%%%%%%%%%%%% SURE using matrix multiplication approach  %%%%%%%%%%%%%%%%%%
%-------------------- Image Denoising -----------------------------------%
Img_out = W_out.*M_nsy;
vec_out = sum(Img_out,2);
den = sum(W_out,2);
v_out = vec_out./den;           % NLM output image of size (N^2, 1), AFTER normalizing.
fNLM = reshape(v_out,N, N)' ;   % NLM output image of size (N, N).     

%------------------ SURE computation ------------------------------------%
x_cap = repmat(v_out, 1, (2*S +1)^2);
y_temp = reshape(f',N^2, 1);
y_i = repmat(y_temp, 1, (2*S+1)^2);
W2 = ( A.*( E + (E + alpha.*A).*exp( -alpha.*(A -phi_m.*E)) ) )./(0.5 *h^2 .* (E + exp( - alpha.*(A - phi_m.*E))).^2 );

% Weights required for 3-rd term of divergence in SURE formula %
W3 = zeros ( N^2, (2*S+1)^2 );
for l=1:N        % Looping over each column
   for i=1:M     % Looping over each row.
       
       for k=-K:K     %looping over search-window column-wise
            for j=-K:K %looping over seach-window row-wise.
                
                row = (1+j+S);
                col = (1+k+S);                 
                W3( (i-1)*N + l, (col-1)*(2*S+1) + row ) = W2( (i-1)*N + l, (col-1)*(2*S+1) + row );               
            end
       end
   end
end      
       
C_i = W_out*ONE;
%tmp = psi_func(1, phi_m, alpha); 
tmp = 1/(1 + exp( -alpha*(1 - phi_m) )); % 1-st term of divergence in SURE formula

nu = W2.*(M_nsy - x_cap).*(M_nsy - y_i) + W3 .*(M_nsy - x_cap).*(Y_swp - y_i);
nu1 = tmp*one + sum(nu,2);
fr = nu1./C_i;
div = sum(fr);           
       
SURE =(   ( sum(sum((fNLM - f).^2))/(N^2) ) - (sig^2) + ( (2* (sig^2) * div )/(N^2))   );
end