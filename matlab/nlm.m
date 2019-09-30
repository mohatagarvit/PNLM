function [fNLM, Weight, M_nsy, Y_swp] = nlm(f, S, K, h)
% --------------------------- ---------------------------------------------
% nlm      : computes original NLM image and the actual NLM weights
%
% -------------- OUTPUT --------------------
% fNLM     : original NLM de-noised image
% Weight   : NLM weight-matrix of size ( N^2, (2*S+1)^2 ).
% M_nsy    : Matrix, which stores noisy image arranged to multiply with corresponding NLM weights. Size = ( N^2, (2*S+1)^2 ).
% Y_swp    : Matrix, which stores noisy image pixel arranged in a way required for SURE MSE computation. Size = ( N^2, (2*S+1)^2 ).
%
% -------------- INPUT --------------------
% l_min    : Lower range of lambda
% S        : search window
% K        : patch 
% h        : smoothing parameter in NLM
%--------------------------------------------------------------------------

[M,N] = size(f);
A = zeros ( N^2, (2*S+1)^2 );       %zero initialization
M_nsy = zeros ( N^2, (2*S+1)^2 );   %zero initialization 
Y_swp = zeros ( N^2, (2*S+1)^2 );   %zero initialization

B=padarray(f,[(K+S),(K+S)],'symmetric'); %zero-oadding
for l=1:N        % Looping over each column
   for i=1:M     % Looping over each row. 

       Pi=B((S+K+i)-K:(S+K+i)+K,(S+K+l)-K:(S+K+l)+K);

        for k=-S:S           %looping over search-window column-wise
            for j=-S:S       %looping over seach-window row-wise.

                Pj=B((S+K+i+j)-K:(S+K+i+j)+K,(S+K+l+k)-K:(S+K+l+k)+K);
                w= (exp(-( (norm( (Pj-Pi),'fro' ))^2 )/(h^2)) ) ;                           

                row = (1+j+S);
                col = (1+k+S);

                A( (i-1)*N + l, (col-1)*(2*S+1) + row ) = w;
                M_nsy( (i-1)*N + l, (col-1)*(2*S+1) + row ) = B((S+K+i+j),(S+K+l+k));
                Y_swp( (i-1)*N + l, (col-1)*(2*S+1) + row ) = B((S+K+i-j),(S+K+l-k));
            end
        end
   end
end
Weight = A;
M_out = A.*M_nsy;
vec_out = sum(M_out,2);
den = sum(A,2);
v_out = vec_out./den;           % NLM output image of size (N^2, 1), AFTER normalizing.
fNLM = reshape(v_out,N, N)' ;  % NLM output image of size (N, N).
end