# # 10 stations n=10
# # the four methods include zidec, cressie, cressie_corrected, cressie_covariance matrix
# library(Matrix)
# library(fields)
# library(matrixcalc)
# library(pracma)
# library(data.table)
# library(rdist)
# # library(matrixcalc)

#' The Isotropic Semivariogram Model.
#'
#' Define a semivariogram function (no nugget effect).
#' Herein, we provide three options, "Exponential", "Power","Gaussian".
#' @param phi Two dimensional vector-valued parameters, \eqn{\phi_1, \phi_2 \geq 0}.
#' @param x The Euclidean distance of two locations, \eqn{x \geq 0}.
#' @param model The semivariogram model, it has three options,
#' "Exponential", "Power", "Gaussian" (default:"Exponential").
#' @return A scalar representing the spatial isotropic variogram estimate for x.
#' @export
#' @references Shanshan Qin, Bin Sun, Yuehua Wu, Yuejiao Fu.(2020)
#' Generalized Least-Squares in Dimension Expansion Method for Nonstationary Processes.
Gamma_phiFun = function(phi, x, model="Exponential"){
  if (model == "Exponential"){
    output = phi[1]*(1 - exp(-x/phi[2]))
  } else if (model == "Power"){
    output = phi[1]*(x^phi[2])
  } else if (model == "Gaussian"){
    output = phi[1]*(1-exp(-x^2/phi[2]))
  }
  return(output)
}

#' Generalized Least-squares (GLS) objective function.
#'
#' The GLS fitting procedure is to minimize the objective function,
#' \deqn{f_{GLS}
#' = \{vec(U^{T})-\gamma_{\phi}W([X,Z])\}^{T}\hat{\Sigma}^{-1}
#'  \{vec(U^{T})-\gamma_{\phi}W([X,Z])\}
#' + \lambda\sum\limits_{k=1}^{p}\|Z_{.k}\|_{1}}
#' @param X The records of location coordinates, a \eqn{n \times d} matrix,
#' where \eqn{n} is the number of locations, and \eqn{d} is the dimension of the record location.
#' Usually, \eqn{d =2} (latitude and longitude coordinates).
#' @param para A vector-valued parameters \eqn{(\phi,Z)} in the objective function, dimension (Maxd-2)n + 2.
#' @param lambda A prespecified tuning parameter, \eqn{\lambda >0}.
#' @param Maxd A prespecified maximum coordinate dimension of each location.
#' @param EmpSemivariogram The empirical semivariogram that is computed by \code{\link{EmpSemivariogramFit}}.
#' @param model The semivariogram model, it has three options, "Exponential", "Power","Gaussian" (default:"Exponential").
#' @return An objective value
#' @export
#' @references Shanshan Qin, Bin Sun, Yuehua Wu, Yuejiao Fu.(2020)
#' Generalized Least-Squares in Dimension Expansion Method for Nonstationary Processes.
ObjectiveGLSFun = function(X, para, lambda, Maxd, EmpSemivariogram, model){
  n = dim(X)[1]
  phi = para[1:2]
  Z = as.matrix(para[-(1:2)],row = n,col = Maxd - 2)
  TempDis = dist(cbind(X,Z))
  TempSemivar = Gamma_phiFun(phi,as.matrix(TempDis),model)
  ######covariance matrix formula (4) #######
  m = choose(n,2)
  CovSemivar=matrix(c(rep(0,m^2)),ncol=m)
  index<-which(upper.tri(diag(1:n)),arr.ind=T)
  index<-as.matrix(data.table::data.table(index, key="row"))
  cressij=function(i,j,ip,jp){((TempSemivar[j,ip]+TempSemivar[i,jp]-TempSemivar[i,ip]-TempSemivar[j,jp])^2)/2}
  for(k in 1:m){
    tempc<-sapply(1:m,function(x) cressij(index[k,1],index[k,2],index[x,1],index[x,2]))
    CovSemivar[k,]<-tempc
  }
  CovSemivarA = pracma::nearest_spd(CovSemivar)
  CovSemivarA= Matrix::nearPD(CovSemivar,corr=FALSE,do2eigen=TRUE)

  u1 = EmpSemivariogram - Gamma_phiFun(phi,as.vector(TempDis),model) # an m-dim vector
  u2 = solve(CovSemivar)    # inverse of the covariance matrix of EmpSemivariogram
  output = t(u1)%*%u2%*%u1 + lambda * sum(sqrt(apply(Z^2,2, function(x) sum(x))))
  #output = t(u1)%*%u2%*%u1 + lambda * sum(apply(Z,2, function(x) sum(abs(x))))
  return(output)
}


#' Weighted Least-squares (WLS) objective function.
#'
#' The WLS fitting procedure is to minimize the objective function,
#' \deqn{f_{WLS}
#' = \sum\limits_{j<i}\frac{1}{\gamma_{\phi}^{2} d_{i,j}([X,Z])}
#'  \{\hat{\gamma}_{i,j}-\gamma_{\phi} d_{i,j}([X,Z]) \}^{2}
#' + \lambda\sum\limits_{k=1}^{p}\|Z_{.k}\|_{1}}
#' @param X The records of location coordinates, a \eqn{n \times d} matrix,
#' where \eqn{n} is the number of locations, and \eqn{d} is the dimension of the record location.
#' Usually, \eqn{d =2} (latitude and longitude coordinates).
#' @param para A vector-valued parameters \eqn{(\phi,Z)} in the objective function, dimension (Maxd-2)n + 2.
#' @param lambda A prespecified tuning parameter, \eqn{\lambda >0}.
#' @param Maxd A prespecified maximum coordinate dimension of each location.
#' @param EmpSemivariogram The empirical semivariogram that is computed by \code{\link{EmpSemivariogramFit}}.
#' @param model The semivariogram model, it has three options, "Exponential", "Power","Gaussian" (default:"Exponential").
#' @return An objective value
#' @export
#' @references Shanshan Qin, Bin Sun, Yuehua Wu, Yuejiao Fu.(2020)
#' Generalized Least-Squares in Dimension Expansion Method for Nonstationary Processes.
ObjectiveWLSFun = function(X, para, lambda, Maxd, EmpSemivariogram,model){
  n = dim(X)[1]
  phi = para[1:2]
  Z = as.matrix(para[-(1:2)],row = n,col = Maxd - 2)
  TempDis = dist(cbind(X,Z))
  TempSemivar = Gamma_phiFun(phi,as.matrix(TempDis),model)
  weight = TempSemivar[lower.tri(TempSemivar)]^(-2)
  u1 = EmpSemivariogram - Gamma_phiFun(phi,as.vector(TempDis),model)
  output = 2*t(weight) %*% (u1^2) + lambda * sum(sqrt(apply(Z^2,2, function(x) sum(x))))
  #output = t(weight) %*% (u1^2) + lambda * sum(apply(Z,2, function(x) sum(abs(x))))
  return(output)
}

#' The objective function in Bornn's paper.
#'
#' The OLS fitting procedure is to minimize the objective function,
#' \deqn{f_{OLS}
#' = \sum\limits_{j<i} \{ \hat{\gamma}_{i,j}-\gamma_{\phi} d_{i,j}([X,Z])\}^{2}
#' + \lambda\sum\limits_{k=1}^{p}\|Z_{.k}\|_{1} }
#' @param X The records of location coordinates, a \eqn{n \times d} matrix,
#' where \eqn{n} is the number of locations, and \eqn{d} is the dimension of the record location.
#' Usually, \eqn{d =2} (latitude and longitude coordinates).
#' @param para A vector-valued parameters \eqn{(\phi,Z)} in the objective function, dimension (Maxd-2)n + 2.
#' @param lambda A prespecified tuning parameter, \eqn{\lambda >0}.
#' @param Maxd A prespecified maximum coordinate dimension of each location.
#' @param EmpSemivariogram The empirical semivariogram that is computed by \code{\link{EmpSemivariogramFit}}.
#' @param model The semivariogram model, it has three options, "Exponential", "Power","Gaussian" (default:"Exponential").
#' @return An objective value
#' @export
#' @references Bornn, L.,Shaddick, G., and Zidek, J.V.(2012). Modeling non-stationary processes through dimension expansion.
#'    Journal of American Statistical Association, 107:281–289.
ObjectiveFunBorrn = function(X, para, lambda, Maxd, EmpSemivariogram,model){
  n = dim(X)[1]
  phi = para[1:2]
  Z = as.matrix(para[-(1:2)],row = n,col = Maxd - 2)
  TempDis = dist(cbind(X, Z))
  u1 = EmpSemivariogram - Gamma_phiFun(phi,as.vector(TempDis),model)
  output = 2* t(u1)%*%u1 + lambda * sum(sqrt(apply(Z^2,2, function(x) sum(x))))
  #output = t(u1)%*%u1 + lambda * sum(apply(Z,2, function(x) sum(abs(x)))) # Lasso
  return(output)
}

#' Estimate the empirical semivariogram function.
#'
#' This is to estimate the empirical semivariogram function.
#'
#' @param Y Realizations of the spatial random process \eqn{Y(X)};
#' @return A vector formed by concatenating the columns of upper triangular of the empirical semivariogram.
#' @export
#' @references Shanshan Qin, Bin Sun, Yuehua Wu, Yuejiao Fu.(2020)
#' Generalized Least-Squares in Dimension Expansion Method for Nonstationary Processes.
EmpSemivariogramFit = function(Y){
  if(dim(Y)[1]>1){
    n = dim(Y)[2]
    CovYObs = cov(as.matrix(Y))
    EmpSemivariogram = matrix(0, ncol=n, nrow=n)
    for(i in 1:n){
      for(j in 1:n){
        EmpSemivariogram[i,j] = CovYObs[i,i] + CovYObs[j,j] - 2*CovYObs[i,j]
      }
    }
    EmpSemivariogram = EmpSemivariogram[lower.tri(EmpSemivariogram)] # a vector
  } else {
    EmpSemivariogram = as.vector(dist(Y)^2)/2
  }
  return(EmpSemivariogram)
}

#' Generate three dimensional geographical data.
#'
#' The spatial random process \eqn{Y(X)} follows a Gaussian process.
#' The locations coordinates are simulated on a three-dimensional half-ellipsoid
#' centered at (0,0,0) and the projection of the first two dimensions
#' is a disk centered at the origin.
#' @param n The number of locations.
#' @param SimRep The number of realizations of the Gaussian process \eqn{Y(X)}.
#' @param trueTheta A given positive value that controls the correlation of the spatial random process.
#' @return
#' \itemize{
#'  \item LocIndex: three dimensional location coordinates;
#'  \item YObs: realizations of the Gaussian process \eqn{Y(X)};
#'  \item SigmaObs: covariance matrix of the random process \eqn{Y(X)}.
#' }
#' @export
#' @references
#' \itemize{
#'  \item Bornn, L.,Shaddick, G., and Zidek, J.V.(2012).
#'    Modeling non-stationary processes through dimension expansion.
#'    Journal of American Statistical Association, 107:281–289.
#'  \item Shanshan Qin, Bin Sun, Yuehua Wu, Yuejiao Fu.(2020).
#'    Generalized Least-Squares in Dimension Expansion Method for Nonstationary Processes.
#' }
GenerateData = function(n,SimRep,trueTheta){
  # n: no.of location obs; m: repeated times; d: dimension
  L = matrix(0, ncol=3, nrow=n)
  for(i in 1:n){
    L[i,3] = runif(1)
    ang = 2*pi*runif(1)
    r = sqrt(1-L[i,3]^2)
    L[i,1] = r * cos(ang) # generate x-axis
    L[i,2] = r * sin(ang) # generate y-axis
    L[i,3] = 10 * L[i,3] # z-axis
  }
  SigmaObs = exp(-fields::rdist(L)/trueTheta)
  YObs =  matrix(rnorm(n*SimRep),SimRep,n) %*% chol(SigmaObs)
  return(list(LocIndex = L,YObs = YObs,SigmaObs = SigmaObs))
}

###### Main Functions of Three methods ############
#' Main fitting procedure.
#'
#' The fitting procedure is to minimize one of the objective functions,
#' \eqn{f_{OLS},f_{WLS},f_{GLS}}.
#' Its purpose is to estimate the semivariogram function and the coordinates of the learned dimensions.
#' The method applied in each iteration is BFGS (Broyden,1979).
#' @param X The records of location coordinates, a \eqn{n \times d} matrix,
#' where \eqn{n} is the number of locations, and \eqn{d} is the dimension of the record location.
#' Usually, \eqn{d =2} (latitude and longitude coordinates).
#' @param Y Realizations of the Gaussian process \eqn{Y(X)}.
#' @param lambda A prespecified tuning parameter, \eqn{\lambda >0}.
#' @param Maxd A prespecified maximum coordinate dimension of each location.
#' @param EmpSemivariogram The empirical semivariogram that is computed by \code{\link{EmpSemivariogramFit}}.
#' @param para.initial Initial parameters of \eqn{(\phi,Z)} in the objective functions, \eqn{f_{OLS},f_{WLS},f_{GLS}}.
#' @param Method Three optional methods, Method = c("OLS","WLS","GLS"), which corresponds to
#' the Least-squares, weighted least-squares, and generalized least-squares fitting procedures (default:OLS).
#' @param model Semivariogram model, it has three options, "Exponential", "Power","Gaussian" (default: "Exponential").
#' @return Return the estimates of \eqn{\phi} and locations coordinates \eqn{Z} of the learned dimensions.
#' @examples n = 10  # number of locations.
#' Maxd = 3
#' SimRep = 1000 # number of observations for each location
#' lambda = 0.01
#' trueTheta = 10
#' set.seed(123)
#' Data = GenerateData(n=n, SimRep = SimRep,trueTheta = trueTheta)
#' Y = Data$YObs
#' X = Data$LocIndex[,1:2]
#' EmpSemivariogram = EmpSemivariogramFit(Y)
#' para.initial = c(1,1,Data$LocIndex[,-c(1:2)])
#' DimExpansion(X,Y,lambda,Maxd,EmpSemivariogram, para.initial, Method="OLS",model = "Exponential")
#' DimExpansion(X,Y,lambda,Maxd,EmpSemivariogram, para.initial, Method="GLS",model = "Exponential")
#' DimExpansion(X,Y,lambda,Maxd,EmpSemivariogram, para.initial, Method="WLS",model = "Exponential")
#' @export
#' @references
#' \itemize{
#'  \item Bornn, L.,Shaddick, G., and Zidek, J.V. (2012). Modeling non-stationary processes through dimension expansion.
#'    Journal of American Statistical Association, 107:281–289.
#'  \item Broyden, C.G.(1979). The convergence of a class of double-rank minimization algorithms.
#'    Journal of the Institute of Mathematics and Its Application,6:76-90.
#'  \item Shanshan Qin, Bin Sun, Yuehua Wu, Yuejiao Fu. (2020)
#' Generalized Least-Squares in Dimension Expansion Method for Nonstationary Processes.
#' }
DimExpansion = function(X,Y,lambda,Maxd,EmpSemivariogram,para.initial,Method="OLS",model="Exponential"){
  if (Method == "OLS") {
    ObjFun = ObjectiveFunBorrn
  } else if (Method == "GLS") {
    ObjFun = ObjectiveGLSFun
  } else if (Method == "WLS") {
    ObjFun = ObjectiveWLSFun
  }
  n = dim(X)[1]
  OptimFit = optim(
    par=para.initial, fn=ObjFun, gr=NULL,
    lambda = lambda, X=X, Maxd=Maxd,
    EmpSemivariogram = EmpSemivariogram,
    model = model,
    method="BFGS",
    control=list(maxit=5000, fnscale=.01))
    phi = OptimFit$par[(1:2)]
    LatentLocation = as.matrix(OptimFit$par[-(1:2)],n,Maxd-2)
    list(phi=phi, LatentLocation = LatentLocation)
}


