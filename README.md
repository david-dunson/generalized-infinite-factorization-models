<a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by-nc/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by-nc/4.0/">Creative Commons Attribution-NonCommercial 4.0 International License</a>.

# GIF_SIS
Bayesian inference of Generalized Infinite Factorization models with Structured Increasing Shrinkage prior

The repository includes R functions to perform Bayesian inference under the Structured Increasing Shrinkage model presented in the following paper:
Schiavon, L., Canale, A., Dunson, D.B. (2022) Generalized infinite factorization models, Biometrika 109 (3), 817-835.
<a rel="paper" href="https://arxiv.org/abs/2103.10333">Link to the pre-print version of the paper</a>
  
File Schiavon_SISGaussian.R refers to the application of the structured increasing shrinkage prior in case of a Gaussian data matrix.

File Schiavon_SIScovariatesRegression.R refers to the application of the structured increasing shrinkage prior when we are interested in factorize a covariate matrix in a linear regression model. 
In particular, it refers to the paper:
Schiavon, L., Canale, A. (2021) Bayesian regularized regression of football tracking data through structured factor models, in Book of Short Papers SIS 2021 (Editors: Pernal, C., Salvati, N. and Schirippa Spagnolo, F.), ISBN: 9788891927361.
<a rel="paper" href="https://air.unimi.it/bitstream/2434/851706/7/pearson-sis-book-2021-parte-1.pdf#page=535">Link to the paper</a>


All files include a function for the MCMC sampler and a function to estimate a meaningful posterior summary of the loadings matrix according to what discussed in Section 3.3 of Schiavon, Canale and Dunson (in press).
