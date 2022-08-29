#' @title Summarize the results after fitting IMIX Spatial Model
#' @description This function maps the results after fitting the spatial IMIX model for the two data types to field effect genes
#'
#' @param input_spatialmodel This is the output 'IMIX_Input_Zscores' after using 'imix_spatial' function 
#' @param fdr The desired FDR threshold for the results, default is 0.1
#'
#' @return 
#' @export
#' @examples
#' fit_imix_spatial=imix_spatial(input1=imix_object_datatype1$IMIX_Input_Zscores,input2=imix_object_datatype2$IMIX_Input_Zscores,model_type = "univariate",input_initial_model = test_uni)
#' fit_imix_spatial_multi=imix_spatial(input1=imix_object_datatype1$IMIX_Input_Zscores,input2=imix_object_datatype2$IMIX_Input_Zscores,model_type = "multivariate",input_initial_model = test_multi)

#To integrate two data types
#Here, mu1 and mu2 are lists of vectors, cov1 and cov2 are lists of matrices
summary_spatial=function(input_spatialmodel, #The saved object after fitting "imix_spatial"
                         fdr=0.1 #The desired FDR threshold for the results, default is 0.1
                         ){
  #Combine the two input z scores, need to make sure that the row names (gene names) are the same
  if(identical(rownames(input1),rownames(input2)) == FALSE) {cat(crayon::red("Error: The gene names of two data types don't match!")); return(1)}
  
  model_type <- match.arg(model_type)
  
  post.prob=int_imix_res$`posterior prob`
  
  load("//q1prpfs04/workspace/zwang21/Aim3/approach6_spark/real_data/map19/nlme_corExp_stat.RData")
  stat.res.matrix2=stat.res.matrix[match(rownames(post.prob),rownames(stat.res.matrix)),]
  load("//q1prpfs04/workspace/zwang21/Aim3/approach6_spark/real_data/rnaseq/map19/nlme_corExp_stat.RData")
  stat.res.matrix1=stat.res.matrix[match(rownames(post.prob),rownames(stat.res.matrix)),]
  
  
  #Field effect for both data types
  id=c((2-1)*8+c(2,5,6,8),(5-1)*8+c(2,5,6,8),(6-1)*8+c(2,5,6,8),(8-1)*8+c(2,5,6,8)) # 16 groups in field effect
  lfdr_field=1-rowSums(post.prob[,id])
  names(lfdr_field)=rownames(post.prob)
  field_fdr=FDR_control_adaptive.test(lfdr_field,alpha=threshold)$significant_genes_with_FDRcontrol
  field_name=names(field_fdr)[which(field_fdr==1)]
  
  
  tmp1=apply(post.prob,1,which.max)
  
  id_field=match(field_name,names(tmp1))
  id_filter=which(tmp1==id[1] | (tmp1==id[2] & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,1])) | #first consider data type 1 is comp2; data type 2 is 2,5,6,8
                    (tmp1==id[3] & sign(stat.res.matrix2[,1])==sign(stat.res.matrix2[,3])) | (tmp1==id[4] & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,3]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,1])) |
                    #second consider data type 1 is comp5; data type 2 is 2,5,6,8
                    (tmp1==id[5]& sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,1])) | (tmp1==id[6] & sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,1]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,1])) |
                    (tmp1==id[7] & sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,1]) & sign(stat.res.matrix2[,1])==sign(stat.res.matrix2[,3])) | (tmp1==id[8] & sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,1]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,3]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,1])) |
                    #third consider data type 1 is comp6; data type 2 is 2,5,6,8
                    (tmp1==id[9] & sign(stat.res.matrix1[,1])==sign(stat.res.matrix1[,3])) | (tmp1==id[10] & sign(stat.res.matrix1[,1])==sign(stat.res.matrix1[,3])& sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,1])) |
                    (tmp1==id[11] & sign(stat.res.matrix1[,1])==sign(stat.res.matrix1[,3])& sign(stat.res.matrix2[,1])==sign(stat.res.matrix2[,3])) | (tmp1==id[12] & sign(stat.res.matrix1[,1])==sign(stat.res.matrix1[,3]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,3]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,1])) |
                    #last consider data type 1 is comp8; data type 2 is 2,5,6,8
                    (tmp1==id[13]& sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,3]) & sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,1])) | (tmp1==id[14]& sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,3]) & sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,1]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,1])) | 
                    (tmp1==id[15]& sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,3]) & sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,1]) & sign(stat.res.matrix2[,1])==sign(stat.res.matrix2[,3])) | (tmp1==id[16]& sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,3]) & sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,1]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,3]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,1])) 
  )
  
  
  id2=intersect(id_field,id_filter)
  
  
  #UC for both data types
  id=c((4-1)*8+4) # 1 group in UC: comp4 for both data
  lfdr_uc=1-post.prob[,id]
  names(lfdr_uc)=rownames(post.prob)
  uc_fdr=FDR_control_adaptive.test(lfdr_uc,alpha=threshold)$significant_genes_with_FDRcontrol
  uc_name=names(uc_fdr)[which(uc_fdr==1)]
  id_uc=match(uc_name,names(tmp1))
  id_filter=which(tmp1==id)
  id4=intersect(id_uc,id_filter)
  
  
  #HG&UC for both data types
  id=c((7-1)*8+7) # 1 group in hguc: comp7 for both data
  lfdr_hguc=1-post.prob[,id]
  names(lfdr_hguc)=rownames(post.prob)
  hguc_fdr=FDR_control_adaptive.test(lfdr_hguc,alpha=threshold)$significant_genes_with_FDRcontrol
  hguc_name=names(hguc_fdr)[which(hguc_fdr==1)]
  id_hguc=match(hguc_name,names(tmp1))
  id_filter=which(tmp1==id & sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,3]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,3]))
  id3=intersect(id_hguc,id_filter)
  
  res6=rep(NA,dim(post.prob)[1])
  #no effect genes
  id1=which(tmp1==1)
  res6[id1]=1
  res6[id2]=2
  res6[id3]=3
  res6[id4]=4
  
  int_imix_group=res6
  int_imix_group[which(is.na(int_imix_group))]=5
  names(int_imix_group)=names(tmp1)
  table(int_imix_group)
  #   1    2    5 
  #197 4167 9324 
  
  id=intersect(names(int_imix_group)[which(int_imix_group==2)],names(uni_imix_group)[which(uni_imix_group==2)])
  length(id) #[1] 3391
  
  summary0=list()
return(summary0)
}
