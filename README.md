# Pruned Non-Local Means (PNLM)
*This is my [Science Academies' Summer Research Fellowship Programme (SRFP) 2017](https://web-japps.ias.ac.in:8443/fellowship2017/lists/selectedList.jsp) project done under [Dr. K. N. Chaudhury](https://sites.google.com/site/kunalnchaudhury/home), Dept. of Electrical Engineering, IISc Bengaluru based on the [Pruned Non-Local Means (PNLM)](https://ieeexplore.ieee.org/document/7932291) paper.*  

Pruned Non-Local Means (PNLM) is a denoising algorithm that discards small weights below a certain threshold (computed using Golden Section Search) in the Non-local Means computation. This works well experimentally and the demo attached is a testament to that. The demo has been written in mex code and thus is much faster than the MATLAB code.  

The details for using the demo have been attached in the [readme](demo/readme.txt) file inside the demo folder.

### About mex code
Mex codes are a mixture of C and MATLAB code. Here the variables are called from MATLAB but the operations are performed in C making it faster than its MATLAB version.  

### Running mex file
First, ensure that a MinGW compiler compatible with your MATLAB version is installed.  
- To install the compiler, use the Add-Ons menu described in Get Add-Ons. Search for MinGW or select from Features.  
- To choose between multiple C or C++ compilers, use `mex -setup` to choose MinGW.  
For more MinGW setup details, refer [this](https://in.mathworks.com/help/matlab/matlab_external/install-mingw-support-package.html).

Finally, compile the mex file. `mex GUI_mex.c`
