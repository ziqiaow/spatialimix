#' @title Fit Multivariate IMIX model
#' @description This function assigns initial values for the input for multivariate IMIX.
#'
#' @param input1 This is the output 'IMIX_Input_Zscores' after using 'CreateSpatialIMIXObejct' function on data type 1. An n x d data frame or matrix of the summary statistics z score, n is the nubmer of genes, d is the number of labels. Each row is a gene, each column is a label within data type 1. 
#' @param input2 This is the output 'IMIX_Input_Zscores' after using 'CreateSpatialIMIXObejct' function on data type 2. An n x d data frame or matrix of the summary statistics z score, n is the nubmer of genes, d is the number of labels. Each row is a gene, each column is a label within data type 2. 
#' @param maxiter The maximum number of iteration, default is 1000
#' @param seed0 set.seed, default is 10
#'
#' @return 
#' @export
#' @examples
#' test_multi=fit_multi(input1=imix_object_datatype1$IMIX_Input_Zscores,input2=imix_object_datatype2$IMIX_Input_Zscores)

#Assign initial values for input for univariate IMIX
fit_multi=function(input1,
                   input2,
                   maxiter=1000,
                   seed0=10
){
  if(identical(rownames(input1),rownames(input2)) == FALSE) {cat(crayon::red("Error: The gene names of two data types don't match!")); return(1)}
  
  fit1=IMIX::IMIX(input1,data_type = "z",model = "IMIX_cor_twostep")
  fit2=IMIX::IMIX(input2,data_type = "z",model = "IMIX_cor_twostep")
  #for data type 1
  mu_data1=fit1$IMIX_cor_twostep$mu
  cov_data1=fit1$IMIX_cor_twostep$cov
  p_data1=fit1$IMIX_cor_twostep$pi
  #for data type 2
  mu_data2=fit2$IMIX_cor_twostep$mu
  cov_data2=fit2$IMIX_cor_twostep$cov
  p_data2=fit2$IMIX_cor_twostep$pi
  
  res=list(fit1,fit2)
  
  out=list("modelfit"=res,"mu_data1"=mu_data1,"mu_data2"=mu_data2,"cov_data1"=cov_data1,"cov_data2"=cov_data2,"p_data1"=p_data1,"p_data2"=p_data2)
  return(out)
}


