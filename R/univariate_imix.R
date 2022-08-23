#' @title Fit Univariate IMIX model
#' @description This function assigns initial values for the input for univariate IMIX. 
#'
#' @param input1 This is the output 'IMIX_Input_Zscores' after using 'CreateSpatialIMIXObejct' function on data type 1. An n x d data frame or matrix of the summary statistics z score, n is the nubmer of genes, d is the number of labels. Each row is a gene, each column is a label within data type 1. 
#' @param input2 This is the output 'IMIX_Input_Zscores' after using 'CreateSpatialIMIXObejct' function on data type 2. An n x d data frame or matrix of the summary statistics z score, n is the nubmer of genes, d is the number of labels. Each row is a gene, each column is a label within data type 2. 
#' @param maxiter The maximum number of iteration, default is 1000
#' @param seed0 set.seed, default is 10
#'
#' @return 
#' @export
#' @examples
#' # A toy example
#' test_uni=fit_uni(input1=imix_object_datatype1$IMIX_Input_Zscores,input2=imix_object_datatype2$IMIX_Input_Zscores)

#Assign initial values for input for univariate IMIX
fit_uni=function(input1,
                 input2,
                 maxiter=1000,
                 seed0=10
                 ){
  if(identical(rownames(input1),rownames(input2)) == FALSE) {cat(crayon::red("Error: The gene names of two data types don't match!")); return(1)}
  
  res=list()
  for(i in 1:dim(input1)[2]){
    res[[i]] = mixtools::normalmixEM(input1[, i], lambda=c(0.9,0.1), mu=c(0,3), sigma=c(1,1), k=2, maxit = maxiter)
    
  }
  for(i in 1:dim(input2)[2]){
    res[[i+3]] = mixtools::normalmixEM(input2[, i], lambda=c(0.9,0.1), mu=c(0,3), sigma=c(1,1), k=2, maxit = maxiter)
    
  }
  
  #########################################
  #Initial values based on single EM
  #########################################
  id_ini1=order(res[[1]]$mu)
  id_ini2=order(res[[2]]$mu)
  id_ini3=order(res[[3]]$mu)
  
  mu1 = res[[1]]$mu[id_ini1]
  mu2 = res[[2]]$mu[id_ini2]
  mu3 = res[[3]]$mu[id_ini3]
  mu_vec = list()
  mu_vec[[1]] = c(mu1[1], mu2[1], mu3[1])
  mu_vec[[2]] = c(mu1[2], mu2[1], mu3[1])
  mu_vec[[3]] = c(mu1[1], mu2[2], mu3[1])
  mu_vec[[4]] = c(mu1[1], mu2[1], mu3[2])
  mu_vec[[5]] = c(mu1[2], mu2[2], mu3[1])
  mu_vec[[6]] = c(mu1[2], mu2[1], mu3[2])
  mu_vec[[7]] = c(mu1[1], mu2[2], mu3[2])
  mu_vec[[8]] = c(mu1[2], mu2[2], mu3[2])
  
  sigma1 = res[[1]]$sigma[id_ini1]
  sigma2 = res[[2]]$sigma[id_ini2]
  sigma3 = res[[3]]$sigma[id_ini3]
  
  cov = list()
  cov[[1]] = diag(x = c(sigma1[1], sigma2[1], sigma3[1]),
                  nrow = 3)
  cov[[2]] = diag(x = c(sigma1[2], sigma2[1], sigma3[1]),
                  nrow = 3)
  cov[[3]] = diag(x = c(sigma1[1], sigma2[2], sigma3[1]),
                  nrow = 3)
  cov[[4]] = diag(x = c(sigma1[1], sigma2[1], sigma3[2]),
                  nrow = 3)
  cov[[5]] = diag(x = c(sigma1[2], sigma2[2], sigma3[1]),
                  nrow = 3)
  cov[[6]] = diag(x = c(sigma1[2], sigma2[1], sigma3[2]),
                  nrow = 3)
  cov[[7]] = diag(x = c(sigma1[1], sigma2[2], sigma3[2]),
                  nrow = 3)
  cov[[8]] = diag(x = c(sigma1[2], sigma2[2], sigma3[2]),
                  nrow = 3)
  
  
  p1 = res[[1]]$lambda[id_ini1]
  p2 = res[[2]]$lambda[id_ini2]
  p3 = res[[3]]$lambda[id_ini3]
  p = c(
    p1[1] * p2[1] * p3[1],
    p1[2] * p2[1] * p3[1],
    p1[1] * p2[2] * p3[1],
    p1[1] * p2[1] * p3[2],
    p1[2] * p2[2] * p3[1],
    p1[2] * p2[1] * p3[2],
    p1[1] * p2[2] * p3[2],
    p1[2] * p2[2] * p3[2]
  )
  
  #for data type 1
  mu_data1=mu_vec
  cov_data1=cov
  p_data1=p
  
  
  
  #########################################
  #Initial values based on single EM
  #########################################
  id_ini1=order(res[[4]]$mu)
  id_ini2=order(res[[5]]$mu)
  id_ini3=order(res[[6]]$mu)
  
  mu1 = res[[4]]$mu[id_ini1]
  mu2 = res[[5]]$mu[id_ini2]
  mu3 = res[[6]]$mu[id_ini3]
  mu_vec = list()
  mu_vec[[1]] = c(mu1[1], mu2[1], mu3[1])
  mu_vec[[2]] = c(mu1[2], mu2[1], mu3[1])
  mu_vec[[3]] = c(mu1[1], mu2[2], mu3[1])
  mu_vec[[4]] = c(mu1[1], mu2[1], mu3[2])
  mu_vec[[5]] = c(mu1[2], mu2[2], mu3[1])
  mu_vec[[6]] = c(mu1[2], mu2[1], mu3[2])
  mu_vec[[7]] = c(mu1[1], mu2[2], mu3[2])
  mu_vec[[8]] = c(mu1[2], mu2[2], mu3[2])
  
  sigma1 = res[[4]]$sigma[id_ini1]
  sigma2 = res[[5]]$sigma[id_ini2]
  sigma3 = res[[6]]$sigma[id_ini3]
  
  cov = list()
  cov[[1]] = diag(x = c(sigma1[1], sigma2[1], sigma3[1]),
                  nrow = 3)
  cov[[2]] = diag(x = c(sigma1[2], sigma2[1], sigma3[1]),
                  nrow = 3)
  cov[[3]] = diag(x = c(sigma1[1], sigma2[2], sigma3[1]),
                  nrow = 3)
  cov[[4]] = diag(x = c(sigma1[1], sigma2[1], sigma3[2]),
                  nrow = 3)
  cov[[5]] = diag(x = c(sigma1[2], sigma2[2], sigma3[1]),
                  nrow = 3)
  cov[[6]] = diag(x = c(sigma1[2], sigma2[1], sigma3[2]),
                  nrow = 3)
  cov[[7]] = diag(x = c(sigma1[1], sigma2[2], sigma3[2]),
                  nrow = 3)
  cov[[8]] = diag(x = c(sigma1[2], sigma2[2], sigma3[2]),
                  nrow = 3)
  
  
  p1 = res[[4]]$lambda[id_ini1]
  p2 = res[[5]]$lambda[id_ini2]
  p3 = res[[6]]$lambda[id_ini3]
  p = c(
    p1[1] * p2[1] * p3[1],
    p1[2] * p2[1] * p3[1],
    p1[1] * p2[2] * p3[1],
    p1[1] * p2[1] * p3[2],
    p1[2] * p2[2] * p3[1],
    p1[2] * p2[1] * p3[2],
    p1[1] * p2[2] * p3[2],
    p1[2] * p2[2] * p3[2]
  )
  
  #for data type 2
  mu_data2=mu_vec
  cov_data2=cov
  p_data2=p
  
  out=list("modelfit"=res,"mu_data1"=mu_data1,"mu_data2"=mu_data2,"cov_data1"=cov_data1,"cov_data2"=cov_data2,"p_data1"=p_data1,"p_data2"=p_data2)
  return(out)
}


