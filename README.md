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

Finally, compile the mex file. `mex GUI_mex.c`. It would look like this.
![Sample MATLAB GUI view](https://github.com/mohatagarvit/PNLM/blob/master/demo/Screenshot.png)

### Links of the Complete Intern Report
- [AuthorCafe](https://edu.authorcafe.com/academies/6714/study-and-implementation-of-non-local-means-and-its-variants)
- [Drive](/)

<!----
- Non-Local Means (NLM) is a standard denoising technique in image processing. Pruned Non-Local Means (PNLM) is a variant of NLM that discards the small neighbourhood weights as they are primarily responsible for the noise. Separable Non-Local Means (SNLM) is a NLM variant that aims to achieve comparable performance to the 2-D NLM technique using 1-D filters.
- The project involved a comprehensive study as well as MATLAB and mex implementation of NLM and PNLM. A GUI with complete abstraction was also developed for demo of PNLM's improvement against NLM for the same parameters.
- After completing the designated tasks for my internship, I was involved in active research in the work of Separable Non-Local Means (SNLM) during the last weeks of my internship. It was when we were trying to remove horizontal and vertical line artifacts that my intern got over. These artifacts were later removed using bilateral filtering and the work published.
---->
