function[f_out, MSE_sure] = pnlm(W_nlm, imgNoisy, S, K, h, sigma, M_nsy, Y_swp, l_min, l_max, param)  

%--------------------------------------------------------------------------
% pnlm     : computes truncated NLM image and SURE of the MSE;
%
% ------------------ OUTPUT ----------------------------------------
% fout     : truncated and smoothed-NLM-denoised-image 
% MSE_sure : Stein's Unbiased Risk Estimate for the Mean Square Error
%
%------------------ INPUT -----------------------------------------
% W_nlm    : NLM weight-matrix of size ( N^2, (2*S+1)^2 ).
% imgNoisy : Noisy image.
% S        : Search window
% K        : Patch size
% h        : Smoothing parameter in NLM
% sigma    : Square-root of the additive gaussian noise variance
% M_nsy    : Matrix, which stores noisy image arranged to multiply with corresponding NLM weights. Size = ( N^2, (2*S+1)^2 ).
% Y_swp    : Matrix, which stores the noisy image pixel arranged in a way required for SURE MSE computation. Size = ( N^2, (2*S+1)^2 ).
% l_min    : Lower-range of lambda
% l_max    : Upper-range of lambda
% param    : Additional Parameters of PNLM
%
%           --- Called Function ------
% nlm_sure : computes truncated NLM image and SURE of the MSE;
%--------------------------------------------------------------------------

itr_max = param.itr_max;
alpha = param.alpha;
ep = param.epsilon;
R = param.R;

x2 = l_max - R*(l_max - l_min);
x1 = l_min + R*(l_max - l_min);

%[fNLM, SURE_mat] = NLM_org_wt_sure(A, f, S, K, h, phi_m, alpha, sig, M_nsy, Y_swp)
[f_out1, MSE_sure1] = pnlm_sure(W_nlm, imgNoisy, S, K, h, x1, alpha, sigma, M_nsy, Y_swp);   
[f_out2, MSE_sure2] = pnlm_sure(W_nlm, imgNoisy, S, K, h, x2, alpha, sigma, M_nsy, Y_swp);   

itr = 0; 
while ( (abs(l_max-l_min)> ep) && (itr<itr_max) )
    itr = itr+1;

    if (MSE_sure1 < MSE_sure2)
        l_min = x2;
        x2 = x1;
        %l_max = l_max;
        x1= l_min +R*(l_max - l_min);            

        [f_out1, MSE_sure1] = pnlm_sure(W_nlm, imgNoisy, S, K, h, x1, alpha, sigma, M_nsy, Y_swp);                         
        [f_out2, MSE_sure2] = pnlm_sure(W_nlm, imgNoisy, S, K, h, x2, alpha, sigma, M_nsy, Y_swp);

    else
        l_max = x1;
        x1 = x2;
        x2=l_max - R*(l_max - l_min);            

        [f_out1, MSE_sure1] = pnlm_sure(W_nlm, imgNoisy, S, K, h, x1, alpha, sigma, M_nsy, Y_swp);                      
        [f_out2, MSE_sure2] = pnlm_sure(W_nlm, imgNoisy, S, K, h, x2, alpha, sigma, M_nsy, Y_swp);     

    end   

end

lambda_str = 0.5*(l_max + l_min);
[f_out, MSE_sure] = pnlm_sure(W_nlm, imgNoisy, S, K, h, lambda_str, alpha, sigma, M_nsy, Y_swp); 