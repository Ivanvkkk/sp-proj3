#Name: Yifan Wu   UUN:s2316499
#Smoothing and function estimation play important parts in applied statistics
#and data science. One approach combines basis expansions and penalized 
#regression. In this project, write R functions for smoothing x, y data.
library(MASS) #loaded library for data

#define x and y data from mcycle dataset
x <- mcycle$times
y <- mcycle$accel

pspline <- function(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100){
  '
  pspline(x,y,k=20,logsp=c(-5,5),bord=3,pord=2,ngrid=100)
    This function is for writing a routine to fit p-spline smoothers 
    to x,y data with GCV smoothing parameter selection
  
  input:
    x,y(float): the vectors of x, y data to smooth
    k(int): the number of basis functions to use
    logsp(float):  the ends of the interval over 
    which to search for the smoothing parameter (log λ scale)
    bord(int): the B-spline order to use
    pord(int): the order of difference to use in the penalty
    ngrid(int): the number of smoothing parameter values to try
    
  output: a list included below elements
    coef(float): the vector of coefficient of function
    fitted(float): the fitted value from function
    sig2(float): the varience
    edf(float): the effective degrees of freedom
    gcv(float):  the smallest generalized cross validation criterion
    r2(float): the coefficient of determination
    lambda(float): smoothing parameter
    V(float):  covariance matrix for the coefficients
    x,y(float): the vectors of x, y data to smooth
  '
  
  dk <- diff(range(x))/(k-bord) ## knot spacing
  #use knots vector used to set up the B-spline basis
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  #use splineDesign to set up the basis (X matrix)
  X <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  #D matrix is obtained by finding the element differences for qr decompose
  D <- diff(diag(k),differences=pord)
  #count the length of y data
  n <- length(y)
  
  #QR decomposition
  Q <- qr.Q(qr(X)) #q matrix
  R <- qr.R(qr(X)) #r matrix
  #get the matrix for eigen-decomposition
  des_matrix <- t(solve(R))%*%t(D)%*%D%*%solve(R)
  
  #eigen-decomposition
  lam <- diag(eigen(des_matrix)$values) #diagonal matrix of eigen-value
  U <- eigen(des_matrix)$vectors #matrix about eigen-vectors
  
  #create a series of arrays to store the details of smooth function
  #at different lambda's
  kapa_list <- c() #store the edf value
  coef_list <- c() #store the coefficient
  fitted_list <- c() #store the fitted value 
  sig2_list <- c() #store the variance
  GCV_list <- c() #store the gcv value
  para_range <- c() #store the lambda value 
  
  #find the range of values for lambda based on the interval
  interval <- (logsp[2]-logsp[1])/(ngrid-1) #Find the range of exponent
  for (i in 0:(ngrid-1)){
    para_range<- append(para_range,exp(logsp[1]+i*interval))
  }
  
  #use the loop to minimize the generalized cross validation criterion and 
  #find the corresponding details of smooth function
  for (lambda in para_range){
    I <- diag(dim(lam)[1]) #create a unit matrix 
    #calculate the effective degrees of freedom
    kapa <- sum(diag((solve(I + lambda*lam)))) 
    #calculate the coefficient of function with QR decomposition
    beta <- solve(R) %*% U %*% solve(I + lambda*lam) %*% t(U) %*% t(Q) %*% y
    #beta <- solve(crossprod(X)+lambda*crossprod(D)) %*% t(X) %*% y 
    #calculate the fitted value
    mi <- X %*% beta
    #calculate the variance
    sig2 <- sum((y - mi)^2) /(n - kapa)
    #calculate the the smallest generalized cross validation criterion
    GCV <- sig2/(n - kapa)
    
    # append each details into the list created above
    kapa_list <- append(kapa_list, kapa)
    coef_list <- append(coef_list, beta)
    fitted_list <- append(fitted_list, mi)
    sig2_list <- append(sig2_list, sig2)
    GCV_list <- append(GCV_list, GCV)
  }
  #find the index of the smallest gcv value
  min_index <- which.min(GCV_list)
  #use the index to get corresponding details of smooth function
  edf <- kapa_list[min_index]
  coef <- coef_list[((min_index-1)*k+1):(min_index*k)]
  fitted <- fitted_list[((min_index-1)*length(y)+1):(min_index*length(y))]
  sig2 <- sig2_list[min_index]
  gcv <- GCV_list[min_index]
  lambda <- para_range[min_index]
  
  #calculate the the coefficient of determination
  r2 <- 1-(n-1)*sig2/sum((y-mean(y))^2)
  #get the residuals standard deviation
  sig <- sqrt(sig2) 
  #get the covariance matrix for the coefficients
  V <- solve(crossprod(X)+lambda*crossprod(D)) *sig2
  #define the class for this function
  class(pspline) <- 'pspline'
  #return the list about details of smooth function
  list(coef=coef,fitted=fitted,sig2=sig2,edf=edf,gcv=gcv,r2=r2,lambda,V,x,y)
  
}

#pspline(x,y)

print.pspline <- function(m) {
  '
  print.pspline(m)
   This function is for printing some statistics data about p-spline smoothers
  
  input:
   m(pspline): the fit p-spline smoothers to x,y data
    
  output: 
   the sentence includes below variable:
    rsd(float): residual standard deviation
    coef(int): length of coefficient
    gcv(float):  the smallest generalized cross validation criterion
    edf(float): the effective degrees of freedom
    r2(float): the coefficient of determination
   silently output list includes:
    gcv, edf, r2(explained above)
  '
  #output the details about the smooth function
  cat(paste('Order 3 p-spline with order 2 penalty\n'))
  cat(paste('Effective degrees of freedom:', unlist(m[4]),' '))
  cat(paste('Coefficients:', length(unlist(m[1])),'\n'))
  cat(paste('residual std dev:', sqrt(unlist(m[3])),' '))
  cat(paste('r-squared:', unlist(m[6]),' '))
  cat(paste('GCV:', unlist(m[5]),'\n'))
  #silently return the list about some statistic data
  invisible(list(gcv=unlist(m[5]),edf=unlist(m[4]),r2=unlist(m[6])))
}

#print.pspline(pspline(x,y))

predict.pspline <- function(m,x,se=TRUE) {
  '
  predict.pspline(m,x,se=True)
   This function is for predicting new x data by using the p-spline smoothers 
   and get the corresponding standard errors from covariance matrix for the 
   coefficients
  
  input:
   m(pspline): the fit p-spline smoothers to x,y data
   x(float): the new value for predicting
   se(bool): indicate if output the standard errors
    
  output: 
   fit(float): the fiited value for new x data
   se(float): the corresponding standard errors
  '
  # define the inital value explained above
  bord <- 3
  pord <- 2
  k <- 20
  dk <- diff(range(x))/(k-bord) ## knot spacing
  #use knots vector used to set up the B-spline basis
  knots <- seq(min(x)-dk*bord,by=dk,length=k+bord+1)
  #use splineDesign to set up the basis (Xp matrix)
  Xp <- splines::splineDesign(knots,x,ord=bord+1,outer.ok=TRUE)
  #D matrix is obtained by finding the element differences for qr decompose
  D <- diff(diag(k),differences=pord)
  #use the coefficients got above
  beta <- matrix(unlist(m[1]))
  #multiple the new x data and coefficients to obtain the new fitted values
  fit <- Xp %*% beta
  
  if (se == FALSE){
    #return the vector of predictions
    fit
  }else{
    #use the covariance matrix for the coefficients got above
    V <-matrix(unlist(m[8]), nrow = 20, ncol = 20)
    #calculate the standard errors by using the row sum of matrix
    se <- rowSums(Xp*(Xp %*% V))^.5  
    #return a named list include fit and se value
    list(fit= fit,se=se)
  }
}

#predict.pspline(pspline(x,y),c(1,2,3))
#predict.pspline(pspline(x,y),c(1,2,3),se=FALSE)

plot.pspline <- function(m) {
  '
  plot.pspline(m)
   This function is for ploting some graph about data, smooth function and 
   residuals, which is for evaluation of fitted models and residual analysis
  
  input:
   m(pspline): the fit p-spline smoothers to x,y data
    
  output: 
   three plots include:
    plot1: original x,y data, estimated smooth function,  approximate 95% 
    credible intervals for the smooths
    plot2: model residuals against fitted values
    plot3: qqplot of the residuals
   silently returns list includes:
    ul(float): lower confidence limits
    ll(float): upper confidence limits
    x(float): corresponding x values
  '
  #define the x,y data
  x = unlist(m[9])
  y = unlist(m[10])
  
  #plot the original x,y data
  plot(x,y,main='Data with Smooth Function')
  #define the upper confidence limits from standard error 
  ul <- predict.pspline(m,x,se=FALSE)+1.96*sqrt(unlist(m[3]))/sqrt(length(y))
  #define the lower confidence limits from standard error 
  ll <- predict.pspline(m,x,se=FALSE)-1.96*sqrt(unlist(m[3]))/sqrt(length(y))
  #plot the smooth function
  lines(x,predict.pspline(m,x,se=FALSE),type="l",col='blue')
  #plot approximate 95% credible intervals for the smooth(between red lines)
  lines(x,ul,type="l",col='red')
  lines(x,ll,type="l",col='red')
  
  #plot the fitted value vs residuals
  fitted_value <- unlist(m[2])
  #subtract the fitted values from the true values to get the residuals 
  residuals <- y - fitted_value 
  plot(fitted_value, residuals, main = 'Fitted Value vs Residuals')
  
  #plot the qqplot 
  qqnorm(residuals)
  #silently returns the list about confidence limits and x
  invisible(list(ll=ll,ul=ul,x=x))
}
#plot.pspline(pspline(x,y))

