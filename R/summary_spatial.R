#' @title Summarize the results after fitting IMIX Spatial Model
#' @description This function maps the results after fitting the spatial IMIX model for the two data types to field effect genes
#'
#' @param input_spatialmodel This is the output 'IMIX_Input_Zscores' after using 'imix_spatial' function 
#' @param threshold The desired FDR threshold for the results, default is 0.1
#'
#' @return 
#' @export
#' @examples
#' res=summary_spatial(fit_imix_spatial,imix_object_datatype1 = imix_object_datatype1,imix_object_datatype2 = imix_object_datatype2)


#To integrate two data types
#Here, mu1 and mu2 are lists of vectors, cov1 and cov2 are lists of matrices
summary_spatial=function(input_spatialmodel, #The saved object after fitting "imix_spatial"
                         threshold=0.1, #The desired FDR threshold for the results, default is 0.1
                         imix_object_datatype1, #The output using 'CreateSpatialIMIXObejct' function on data type 1
                         imix_object_datatype2 # #The output using 'CreateSpatialIMIXObejct' function on data type 2
                         ){
  
  post.prob=input_spatialmodel$`posterior prob`
  stat.res.matrix1=imix_object_datatype1$spatial_mixed_model_result[match(rownames(post.prob),rownames(imix_object_datatype1$spatial_mixed_model_result)),]
  stat.res.matrix2=imix_object_datatype2$spatial_mixed_model_result[match(rownames(post.prob),rownames(imix_object_datatype2$spatial_mixed_model_result)),]
  
  
  #Field effect for both data types
  id=c((2-1)*8+c(2,5,6,8),(5-1)*8+c(2,5,6,8),(6-1)*8+c(2,5,6,8),(8-1)*8+c(2,5,6,8)) # 16 groups in field effect
  lfdr_field=1-rowSums(post.prob[,id])
  names(lfdr_field)=rownames(post.prob)
  field_fdr=IMIX::FDR_control_adaptive(lfdr_field,alpha=threshold)$significant_genes_with_FDRcontrol
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
  uc_fdr=IMIX::FDR_control_adaptive(lfdr_uc,alpha=threshold)$significant_genes_with_FDRcontrol
  uc_name=names(uc_fdr)[which(uc_fdr==1)]
  id_uc=match(uc_name,names(tmp1))
  id_filter=which(tmp1==id)
  id4=intersect(id_uc,id_filter)
  
  
  #HG&UC for both data types
  id=c((7-1)*8+7) # 1 group in hguc: comp7 for both data
  lfdr_hguc=1-post.prob[,id]
  names(lfdr_hguc)=rownames(post.prob)
  hguc_fdr=IMIX::FDR_control_adaptive(lfdr_hguc,alpha=threshold)$significant_genes_with_FDRcontrol
  hguc_name=names(hguc_fdr)[which(hguc_fdr==1)]
  id_hguc=match(hguc_name,names(tmp1))
  id_filter=which(tmp1==id & sign(stat.res.matrix1[,2])==sign(stat.res.matrix1[,3]) & sign(stat.res.matrix2[,2])==sign(stat.res.matrix2[,3]))
  id3=intersect(id_hguc,id_filter)
  
  res6=rep(NA,dim(post.prob)[1])
  #no effect genes
  id1=which(tmp1==1)
  res6[id1]=1 #no field effect genes, others
  res6[id2]=2 #field effect genes
  res6[id3]=3 #HG&UC genes
  res6[id4]=4 #UC only genes
  
  int_imix_group=res6
  int_imix_group[which(is.na(int_imix_group))]=1 #other genes
  names(int_imix_group)=names(tmp1)
  results=data.frame(int_imix_group)
  results$label=0
  results$label[which(results$int_imix_group==1)]="Not significant"
  results$label[which(results$int_imix_group==2)]="Field effect genes"
  results$label[which(results$int_imix_group==3)]="HG&UC genes"
  results$label[which(results$int_imix_group==4)]="UC only genes"
  label_order=c("Field effect genes","HG&UC genes","UC only genes","Not significant")
  results=results[order(factor(results$label, levels = label_order)),]
  summary0=list(results=results,FDR=threshold)
return(summary0)
}
