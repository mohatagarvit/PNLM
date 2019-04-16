#include"mex.h"
#include"matrix.h"
#include"string.h"
#include"math.h"

#define NSY_IMG_IN         prhs[0]
#define SIGMA_IN           prhs[1]
#define S_IN               prhs[2]
#define K_IN               prhs[3]
#define alpha_IN           prhs[4]
#define h_IN               prhs[5]
#define IMAGE_OUT          plhs[0]
#define DIFF_X_HAT_WITH_Y  plhs[1]

/* ========================================================================
 * Mex code for performing computation of denoised image and its partial 
 * derivative w.r.t. y(i)
 *
 * INPUTS:
 *
 * y_padded : noisy image array,y symmetrically padded to [S+K,S+K] dimensions
 * sigma    : noise level of the noisy image
 * S        : Width/extension of the Search window of pixel i in any direction from it.
 *            The dimensions of its search window are [2*S+1,2*S+1]. 
 * K        : Width/extension of the patch of pixel i in any direction from it. The
 *            dimensions of its patch are [2*K+1,2*K+1].
 * alpha    : the slope of the sigmoid function, phi at the threshold weight, lambda.
 * h        : a smoothing parameter, which is a function of sigma.
 *
 * OUTPUTS:
 *
 * x_hat             : the pnlm-processed denoised image
 * diff_x_hat_with_y : sum of the partial derivative terms of x_hat with 
 *                     y(i) for every i
 * 
 * OTHER Variables:
 *
 * lambda: optimal value of a threshold for weight that is being used in calculation of denoised image
 *
 *======================================================================*/

/*The PHI function*/
double phi(double weight,double alpha,double optimal_lambda)     
{
    double out=1/(1+exp(-alpha*(weight-optimal_lambda)));
    return out;
}

/*MAIN MEX FUNCTION*/
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    /* Variable declarations */
    
        double sigma,S,K,alpha,h;   
        double *y_padded,*x_hat,**wght; 
        double *sumofwght;
        double  wghtd_sum=0;
        mwSize mrows,ncols,buf_len; 
        int i1,i2,j1,j2,k1,k2,j_dash1,j_dash2,i1_dash,i2_dash,j1_wght,j2_wght;
        int i,j;
        double a,sum,weight;
        int m,n,counter=0,s;
        double *diff_x_hat_with_y;
        double final=0,xi1=0;
        double *x_hat_p,*x_hat_q;
        
                
        /* Checking the inputs*/
    
       /*Number of i/o arguments*/
        if (nrhs !=6 )
            mexErrMsgTxt("Six inputs required: The input padded-noisy image,sigma,S,K,alpha,h");
        else if (nlhs !=2)
            mexErrMsgTxt("Too many output arguments. Only two required. X-hat and diff_x_hat_with_y.");
       
       /*Datatype consistency*/
        if (!mxIsNumeric(prhs[0]))   
            mexErrMsgTxt("Input image must be a padded noisy image.");
        if (!mxIsDouble(prhs[1])||mxGetNumberOfElements(prhs[1])!=1||mxIsComplex(prhs[1])||(sigma = (double)mxGetScalar(SIGMA_IN)) <= 0)   
            mexErrMsgTxt("Sigma must be a positive scalar double.");
        if (!mxIsDouble(prhs[2])||mxGetNumberOfElements(prhs[2])!=1||mxIsComplex(prhs[2])||(S = (double)mxGetScalar(S_IN)) <= 0||S!=(int)S)
            mexErrMsgTxt("S must be a positive scalar integer."); 
        if (!mxIsDouble(prhs[3])||mxGetNumberOfElements(prhs[3])!=1||mxIsComplex(prhs[3])||(K = (double)mxGetScalar(K_IN)) <= 0||K!=(int)K)
            mexErrMsgTxt("K must be a positive scalar integer.");   
        if (!mxIsDouble(prhs[4])||mxGetNumberOfElements(prhs[4])!=1||mxIsComplex(prhs[4])||(alpha = (double)mxGetScalar(alpha_IN)) <= 0)
            mexErrMsgTxt("alpha must be a positive scalar double.");   
        if (!mxIsDouble(prhs[5])||mxGetNumberOfElements(prhs[5])!=1||mxIsComplex(prhs[5])||(h = (double)mxGetScalar(h_IN)) <= 0)
            mexErrMsgTxt("h must be a positive scalar double."); 
            
        /* Memory Allocation and Retrieving input arguments */
        
        mrows=mxGetM(prhs[0]);
        ncols=mxGetN(prhs[0]);
        m=mrows-2*(S+K); 
        n=ncols-2*(S+K);
        buf_len=mrows*ncols;

        y_padded=mxCalloc(buf_len,sizeof(double)); 
        y_padded=mxGetPr(prhs[0]);
        sigma=mxGetScalar(prhs[1]); 
        S=mxGetScalar(prhs[2]);
        K=mxGetScalar(prhs[3]);
        alpha=mxGetScalar(prhs[4]);
        h=mxGetScalar(prhs[5]);
        
        IMAGE_OUT = mxCreateDoubleMatrix(m, n, mxREAL);
        x_hat=mxCalloc( m*n ,sizeof(double));
        x_hat=mxGetPr(plhs[0]);
        
        x_hat_p=mxCalloc( m*n ,sizeof(double));
        x_hat_q=mxCalloc( m*n ,sizeof(double));
        
        plhs[1]=mxCreateDoubleScalar(0);
        diff_x_hat_with_y=mxCalloc(1,sizeof(double));
        diff_x_hat_with_y= mxGetPr(plhs[1]);
        *diff_x_hat_with_y=0;
        
        /* weight vector allocation */
        wght=mxCalloc(m*n,sizeof(double));
        for(i1=S+K;i1<S+K+m;i1++)   
        {   
            for(i2=S+K;i2<S+K+n;i2++)
            {
                i1_dash=i1-(S+K);/*0-m*/
                i2_dash=i2-(S+K);/*0-n*/
                *(wght+i1_dash*n+i2_dash)=mxCalloc((2*S+1)*(2*S+1),sizeof(double));
            }
        }
        
        sumofwght=mxCalloc(m*n,sizeof(double));/*calloc initialises to 0*/
        s=2*S+1;
        
      /* Weight Computation*/
                                                                           
     double xi,value2,value3;   
     
     for(i1=S+K;i1<S+K+m;i1++)   
     {   
        for(i2=S+K;i2<S+K+n;i2++)
        {
            i1_dash=i1-(S+K);
            i2_dash=i2-(S+K);
            
            for(j_dash1=-S;j_dash1<=S;j_dash1++)
            {
                for(j_dash2=-S;j_dash2<=S;j_dash2++)
                {
                    j1=i1+j_dash1;
                    j2=i2+j_dash2;
                    j1_wght=j_dash1+S;
                    j2_wght=j_dash2+S;
                    
                    sum=0;
                    for(k1=-K;k1<=K;k1++)
                    {   
                        for(k2=-K;k2<=K;k2++)
                        {
                            a=*(y_padded+(i1+k1)*ncols+(i2+k2))-*(y_padded+(j1+k1)*ncols+(j2+k2)); 
                            a=a*a;
                            sum+=a;
                        }                        
                    }
                    *(*(wght+i1_dash*n+i2_dash)+j1_wght*s+j2_wght)=exp(-sum/(h*h)); 
                }
            }
        }
    } 
    
    /* GOLDEN SECTION SEARCH  */
     
     /*Variable declarartion of variables being used in the search*/
     double lambda0,ro,delta,low_lim,up_lim,p,q,  e,f ,sum_p,sum_q,lambda,exp_value,diff_xi;
     double diff_x_hat_with_y_p,diff_x_hat_with_y_q,value2_p,value3_p,value2_q,value3_q,SURE_p,SURE_q;
     int flag=0,loop_count=0, k_dash1,k_dash2;

     lambda0   =  4.3*(1e-7)*pow(sigma,3) - 1.1*(1e-4)*pow(sigma,2) + 9.2*(1e-3)*sigma + 0.039;
     ro        =  0.618;
     delta     =  0.05;
     loop_count=  0;
            
     low_lim   =  lambda0 - delta;   
     up_lim    =  lambda0 + delta;
            
     p         =  up_lim - ro*(up_lim - low_lim);   
     q         =  low_lim + ro*(up_lim - low_lim);

     while (up_lim - low_lim > 1e-2)
     {   
          diff_x_hat_with_y_p=0;
          diff_x_hat_with_y_q=0;
    
          /*For p*/
          
          if(flag==0 || loop_count==0) 
          {  
                 
             /*Calculating x_hat_p*/
                 
              for(i1=S+K;i1<S+K+m;i1++)   
              {   
                  for(i2=S+K;i2<S+K+n;i2++)
                  {
                      i1_dash=i1-(S+K);
                      i2_dash=i2-(S+K);
                      wghtd_sum=0;
                      for(j_dash1=-S;j_dash1<=S;j_dash1++)
                      {
                          for(j_dash2=-S;j_dash2<=S;j_dash2++)
                          {
                              j1      =  i1+j_dash1;
                              j2      =  i2+j_dash2;
                              j1_wght =  j_dash1+S;
                              j2_wght =  j_dash2+S;
                                
                              weight=*(*(wght+i1_dash*n+i2_dash)+j1_wght*s+j2_wght);
                              xi    = weight  *  phi(weight,alpha,p);

                              *(sumofwght+i1_dash*n+i2_dash)+=xi;
                              wghtd_sum+=xi**(y_padded+j1*ncols+j2);                
                          }
                      }
                      *(x_hat_p+i1_dash*n+i2_dash)=wghtd_sum/(*(sumofwght+i1_dash*n+i2_dash));
                  }
              } 
                
              sum_p=0; 
              for(i1=S+K;i1<S+K+m;i1++)   
              {   
                  for(i2=S+K;i2<S+K+n;i2++)
                  {
                      i1_dash=i1-(S+K);
                      i2_dash=i2-(S+K);
                      e     = *(x_hat_p+i1_dash*n+i2_dash)-*(y_padded + ncols*i1 + i2); 
                      e     = e*e;
                      sum_p = sum_p + e;
                  }
              } 
              /*Calculating value2_p and value3_p*/
                
               for(i1=S+K;i1<S+K+m;i1++)   
               {   
                   for(i2=S+K;i2<S+K+n;i2++)
                   {
                       i1_dash=i1-(S+K);
                       i2_dash=i2-(S+K);
                       value2_p=0;
                       value3_p=0;
                       
                       for(j_dash1=-S ; j_dash1<=S ; j_dash1++)
                       {
                           for(j_dash2=-S ; j_dash2<=S ; j_dash2++)
                           {
                               j1      =  i1+j_dash1;
                               j2      =  i2+j_dash2;
                               j1_wght =  j_dash1+S;
                               j2_wght =  j_dash2+S;

                               weight     =  *(*(wght+i1_dash*n+i2_dash)+j1_wght*s+j2_wght);
                               exp_value  =  exp(-alpha*(weight-p));
                               diff_xi    =  (1+(1+alpha*weight )*exp_value)/((1+exp_value)*(1+exp_value));
                               value2_p   =  value2_p+ 2/(h*h)*diff_xi * weight * (*(y_padded+j1*ncols+j2)-*(y_padded + i1*ncols+i2))*(*(y_padded+j1*ncols+j2)-*(x_hat_p + i1_dash*n+i2_dash));

                               k_dash1=i1-j1; k_dash2=i2-j2;
                               if(k_dash1<=K && k_dash1>=-K && k_dash2>=-K && k_dash2<=K)
                               {
                                   value3_p =  value3_p + 2/(h*h) *diff_xi* weight  *(*(y_padded+(i1+k_dash1)*ncols+(i2+k_dash2))-*(y_padded + i1*ncols+i2))*(*(y_padded+(i1-k_dash1)*ncols+(i2-k_dash2))-*(x_hat_p + i1_dash*n+i2_dash));
                               }
                           }
                       }
                       xi1  = phi(1,alpha,p);
                       diff_x_hat_with_y_p =  diff_x_hat_with_y_p + (xi1 + value2_p + value3_p )/(*(sumofwght+i1_dash*n+i2_dash));
                   }
               }
               SURE_p=sum_p/(m*n) - sigma*sigma + 2*sigma*sigma*diff_x_hat_with_y_p/(m*n) ;
               
                
               /*Making sumofwght again 0*/
                
               for(i1=S+K;i1<S+K+m;i1++)   
               {   
                   for(i2=S+K;i2<S+K+n;i2++)
                   {
                       i1_dash=i1-(S+K);/*0-m*/
                       i2_dash=i2-(S+K);/*0-n*/
                       *(sumofwght+i1_dash*n+i2_dash)=0;
                   }
               }
          }
                 
          /*For q*/
          
          if(flag==1 || loop_count==0)     
          {    
              /*Calculating x_hat_q*/
                 
              for(i1=S+K;i1<S+K+m;i1++)   
              {   
                  for(i2=S+K;i2<S+K+n;i2++)
                  {
                      i1_dash=i1-(S+K);
                      i2_dash=i2-(S+K);
                      wghtd_sum=0;
                      for(j_dash1=-S;j_dash1<=S;j_dash1++)
                      {
                          for(j_dash2=-S;j_dash2<=S;j_dash2++)
                          {
                              j1      =  i1+j_dash1;
                              j2      =  i2+j_dash2;
                              j1_wght =  j_dash1+S;
                              j2_wght =  j_dash2+S;
                                
                              weight=*(*(wght+i1_dash*n+i2_dash)+j1_wght*s+j2_wght);
                              xi  = weight * phi(weight,alpha,q);

                              *(sumofwght+i1_dash*n+i2_dash)+=xi;
                              wghtd_sum+=xi**(y_padded+j1*ncols+j2);                
                          }
                      }
                      *(x_hat_q+i1_dash*n+i2_dash)=wghtd_sum/(*(sumofwght+i1_dash*n+i2_dash));
                  }
              } 
              
              sum_q=0;
              for(i1=S+K;i1<S+K+m;i1++)   
              {   
                  for(i2=S+K;i2<S+K+n;i2++)
                  {
                      i1_dash=i1-(S+K);
                      i2_dash=i2-(S+K);
                      f     = *(x_hat_q+i1_dash*n+i2_dash)-*(y_padded + ncols*i1 + i2);
                      f     = f*f;
                      sum_q = sum_q + f;
                  }
              } 
              /*Calculating value2_q and value3_q*/
              
              for(i1=S+K;i1<S+K+m;i1++)   
              {   
                  for(i2=S+K;i2<S+K+n;i2++)
                  {
                      i1_dash=i1-(S+K);
                      i2_dash=i2-(S+K);
                      value2_q=0;
                      value3_q=0;
                      for(j_dash1=-S ; j_dash1<=S ; j_dash1++)
                      {
                          for(j_dash2=-S ; j_dash2<=S ; j_dash2++)
                          {
                              j1      =  i1+j_dash1;
                              j2      =  i2+j_dash2;
                              j1_wght =  j_dash1+S;
                              j2_wght =  j_dash2+S;

                              weight     =  *(*(wght+i1_dash*n+i2_dash)+j1_wght*s+j2_wght);
                              exp_value  =  exp(-alpha*(weight-q));
                              diff_xi    =  (1+(1+alpha*weight )*exp_value)/((1+exp_value)*(1+exp_value));
                              value2_q   =  value2_q+ 2/(h*h)*diff_xi * weight * (*(y_padded+j1*ncols+j2)-*(y_padded + i1*ncols+i2))*(*(y_padded+j1*ncols+j2)-*(x_hat_q + i1_dash*n+i2_dash));

                              k_dash1=i1-j1; k_dash2=i2-j2;
                              if(k_dash1<=K && k_dash1>=-K && k_dash2>=-K && k_dash2<=K)
                              {
                                 value3_q =  value3_q + 2/(h*h) *diff_xi* weight *(*(y_padded+(i1+k_dash1)*ncols+(i2+k_dash2))-*(y_padded + i1*ncols+i2))*(*(y_padded+(i1-k_dash1)*ncols+(i2-k_dash2))-*(x_hat_q + i1_dash*n+i2_dash));
                              }
                          }
                      }
                      xi1  = phi(1,alpha,q);
                      diff_x_hat_with_y_q =  diff_x_hat_with_y_q + (xi1 + value2_q + value3_q )/(*(sumofwght+i1_dash*n+i2_dash));
                  }
              }
              SURE_q=sum_q/(m*n)  - sigma*sigma + 2*sigma*sigma*diff_x_hat_with_y_q/(m*n);
                
              /*Making sumofwght=0*/
                
              for(i1=S+K;i1<S+K+m;i1++)   
              {   
                  for(i2=S+K;i2<S+K+n;i2++)
                  {
                      i1_dash=i1-(S+K);/*0-m*/
                      i2_dash=i2-(S+K);/*0-n*/
                      *(sumofwght+i1_dash*n+i2_dash)=0;
                  }
              }
          }
              
          if ( SURE_p > SURE_q ) 
          {    low_lim=p;
               p=q;
               q=low_lim + ro*(up_lim - low_lim);
               SURE_p=SURE_q;
               flag=1;
          }     
          else 
          {    up_lim=q;
               q=p;
               p=up_lim - ro*(up_lim - low_lim);
               SURE_q=SURE_p;
               flag=0;
          }
          loop_count=loop_count+1;
     }
     lambda=(up_lim+low_lim)/2;
     
    /*Final Computation*/  
     
    for(i1=S+K;i1<S+K+m;i1++)   
    {   
        for(i2=S+K;i2<S+K+n;i2++)
        {
            i1_dash=i1-(S+K);
            i2_dash=i2-(S+K);
            wghtd_sum=0;
            for(j_dash1=-S;j_dash1<=S;j_dash1++)
            {
                for(j_dash2=-S;j_dash2<=S;j_dash2++)
                {
                    j1      =  i1+j_dash1;
                    j2      =  i2+j_dash2;
                    j1_wght =  j_dash1+S;
                    j2_wght =  j_dash2+S;
                    
                    weight     =  *(*(wght+i1_dash*n+i2_dash)+j1_wght*s+j2_wght);
                    xi  = *(*(wght+i1_dash*n+i2_dash)+j1_wght*s+j2_wght)  *  phi(weight,alpha,lambda);
                    
                    *(sumofwght+i1_dash*n+i2_dash)+=xi;
                    wghtd_sum+=xi**(y_padded+j1*ncols+j2);                
                }
            }
            *(x_hat+i1_dash*n+i2_dash)=wghtd_sum/(*(sumofwght+i1_dash*n+i2_dash)); 
        }
    } 
            
   for(i1=S+K;i1<S+K+m;i1++)   
   {   
       for(i2=S+K;i2<S+K+n;i2++)
       {
           i1_dash=i1-(S+K);
           i2_dash=i2-(S+K);
           value2=0;
           value3=0;
           for(j_dash1=-S ; j_dash1<=S ; j_dash1++)
           {
               for(j_dash2=-S ; j_dash2<=S ; j_dash2++)
               {
                   j1      =  i1+j_dash1;
                   j2      =  i2+j_dash2;
                   j1_wght =  j_dash1+S;
                   j2_wght =  j_dash2+S;
                   
                   weight     =  *(*(wght+i1_dash*n+i2_dash)+j1_wght*s+j2_wght);                    
                   exp_value  =  exp(-alpha* (weight-lambda));
                   diff_xi    =  (1+(1+alpha*weight )*exp_value)/((1+exp_value)*(1+exp_value));
                   value2     =  value2+ 2/(h*h)*diff_xi * weight * (*(y_padded+j1*ncols+j2)-*(y_padded + i1*ncols+i2))*(*(y_padded+j1*ncols+j2)-*(x_hat + i1_dash*n+i2_dash));
                    
                   k_dash1=i1-j1; k_dash2=i2-j2;
                   if(k_dash1<=K && k_dash1>=-K && k_dash2>=-K && k_dash2<=K)
                   {
                       value3 = value3 + 2/(h*h) *diff_xi* weight *(*(y_padded+(i1+k_dash1)*ncols+(i2+k_dash2))-*(y_padded + i1*ncols+i2))*(*(y_padded+(i1-k_dash1)*ncols+(i2-k_dash2))-*(x_hat + i1_dash*n+i2_dash));
                   }
                }
            }
            xi1  = phi(1,alpha,lambda);
            final =  final + (xi1 + value2 + value3 )/(*(sumofwght+i1_dash*n+i2_dash));
       }
   }
   *diff_x_hat_with_y=final;
}