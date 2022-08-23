#' @title Fit IMIX Spatial Model
#' @description This function fits the spatial IMIX model
#'
#' @param input1 This is the output 'IMIX_Input_Zscores' after using 'CreateSpatialIMIXObejct' function on data type 1. An n x d data frame or matrix of the summary statistics z score, n is the nubmer of genes, d is the number of labels. Each row is a gene, each column is a label within data type 1. 
#' @param input2 This is the output 'IMIX_Input_Zscores' after using 'CreateSpatialIMIXObejct' function on data type 2. An n x d data frame or matrix of the summary statistics z score, n is the nubmer of genes, d is the number of labels. Each row is a gene, each column is a label within data type 2. 
#' @param maxiter The maximum number of iteration, default is 1000
#' @param seed0 set.seed, default is 10
#'
#' @return 
#' @export
#' @examples
#' fit_imix_spatial=imix_spatial(input1=imix_object_datatype1$IMIX_Input_Zscores,input2=imix_object_datatype2$IMIX_Input_Zscores,model_type = "univariate",input_initial_model = test_uni)
#' fit_imix_spatial_multi=imix_spatial(input1=imix_object_datatype1$IMIX_Input_Zscores,input2=imix_object_datatype2$IMIX_Input_Zscores,model_type = "multivariate",input_initial_model = test_multi)

#To integrate two data types
#Here, mu1 and mu2 are lists of vectors, cov1 and cov2 are lists of matrices
imix_spatial=function(input1,
                      input2,
                      model_type=c("univariate","multivariate"), #Whether to use univariate IMIX or multivariate IMIX, default is univaraite IMIX
                      input_initial_model, #If model_type is 'univariate', this is the output from 'fit_uni'; if model_type is 'multivariate', this is the output from 'fit_multi'.
                      g_spatial=64, #The number of components, default is 64
                      tol0=1e-6, #The convergence criterion. Convergence is declared when the change in the observed data log-likelihood increases by less than epsilon.
                      maxiter0=1000, #The maximum number of iteration, default is 1000
                      seed0=10,#set.seed, default is 10
                      verbose0=FALSE){
  #Combine the two input z scores, need to make sure that the row names (gene names) are the same
  if(identical(rownames(input1),rownames(input2)) == FALSE) {cat(crayon::red("Error: The gene names of two data types don't match!")); return(1)}
  
  model_type <- match.arg(model_type)
 
  mu1=input_initial_model$mu_data1
  mu2=input_initial_model$mu_data2
  cov1=input_initial_model$cov_data1
  cov2=input_initial_model$cov_data2
  p1=input_initial_model$p_data1
  p2=input_initial_model$p_data2 
  
  
  input_tmp=cbind(input1,input2)
  mu_tmp=list()
  cov_tmp=list()
  p_tmp=0
  for(i in 1:length(mu1)){
    for(j in 1:length(mu2)){
      #Put the mean vectors together
      mu_tmp[[(i-1)*length(mu1)+j]]=c(mu1[[i]],mu2[[j]])
      
      #Bind the two covariance matrices together in block diagnol
      cov_tmp[[(i-1)*length(mu1)+j]]=rbind(cbind(cov1[[i]], matrix(0, nrow=nrow(cov1[[i]]), ncol=ncol(cov2[[j]]))),
                                           cbind(matrix(0, nrow=nrow(cov2[[j]]), ncol=ncol(cov1[[i]])), cov2[[j]]))
      
      #Calculate the initial pi
      p_tmp[(i-1)*length(mu1)+j]=p1[i]*p2[j]
    }
  }
  
  fit=IMIX_pi(input_tmp,data_type="z",g=g_spatial,mu_vec=mu_tmp,cov=cov_tmp,p=p_tmp,tol=tol0,maxiter=maxiter0,seed=seed0,verbose=verbose0)
  return(fit)
}
