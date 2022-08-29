#' @title CreateSpatialIMIXObject
#' @description Create Data Input for IMIX model for one data type
#'
#' @param ratio An n x d data frame or matrix of the log2ratio between the samples and the controls
#' @param label A vector of lenghth d corresponding to the sample subtypes of the samples
#' @param location Sample corrdinates to compute the kernel matrix, d x 2 data matrix (d is the number of samples, the columns correspond to the x and y coordinates of the sample on a 2-dimentional space)
#' @param reference_label The reference label for input in the spatial mixed model, here the example is LG (low-grade intraurothelial neoplasia/normal-appearing urothelial) 
#' @param second_label The second layer label for input in the spatial mixed model, here the example is HG (high-grade intraurothelial neoplasia) 
#' @param top_label The top layer label for input in the spatial mixed model, here the example is UC (urothelial carcinoma) 
#' @return The data input ready for spatial IMIX
#'  \item{IMIX_Input_Zscores}{Z scores input for IMIX model, this is inverse normal transformed from the p-values using LRT test on the spatial mixed model results}
#'  \item{spatial_mixed_model_result}{A data matrix of results fitting the spatial mixed model, including the estimated values, SE, DF, t-value, p-value for each gene of LG, HG, UC. The output LG, HG, UC correspond to reference_label, second_label, top_label respectively}
#'  \item{full_model}{List of results fitting the full spatial mixed model with variables including HG, UC and the intercept, these are used to construct the LRT test statistics and p-values}
#'  \item{model_withouthg}{List of results fitting the spatial mixed model with the intercept and UC but without HG, these are used to construct the LRT test statistics and p-values}
#'  \item{model_withoutintercept}{List of results fitting the spatial mixed model with HG and UC but without the intercept, these are used to construct the LRT test statistics and p-values}
#'  \item{model_withoutuc}{List of results fitting the spatial mixed model with the intercept and HG but without UC, these are used to construct the LRT test statistics and p-values}
#' @importFrom stats qnorm
#' @export
#' 
#' @examples 
#' # A toy example
#' data("example_spatial_IMIX")
#' imix_object_datatype1 <- CreateSpatialIMIXObject(ratio=ratio1,label=label,location=location)
#' imix_object_datatype2 <- CreateSpatialIMIXObject(ratio=ratio2,label=label,location=location)
CreateSpatialIMIXObject <- function(ratio,
                                    location,
                                    label,
                                    reference_label="LG",
                                    second_label="HG",
                                    top_label="UC"
                                    ){
reference_label <- match.arg(reference_label)
second_label <- match.arg(second_label)
top_label <- match.arg(top_label)

info <- data.frame(cbind(label,location))
colnames(info) <- c("group","x","y")
if (dim(ratio)[2] != dim(info)[1] ) {
  cat(crayon::red("Error: The number of samples is different in the data matrix and the label/location!"))
  return(1)
}
rownames(info) <- colnames(ratio)
info$group=as.character(info$group)
info$x=as.numeric(as.character(info$x))
info$y=as.numeric(as.character(info$y))

#Change the label names to LG, HG, UC
info$group[which(info$group == reference_label)] = "LG"
info$group[which(info$group == second_label)] = "HG"
info$group[which(info$group == top_label)] = "UC"
if (length(which(info$group %in% c("HG","LG","UC") == "FALSE")) > 0) {
  cat(crayon::red("Error: Labels do not match with the names specified in the function input!"))
  return(1)
}

#create dummy variables
info$LG <- ifelse(info$group == 'LG', 1, 0)
info$HG <- ifelse(info$group == 'HG', 1, 0)
info$UC <- ifelse(info$group == 'UC', 1, 0)

# fit the model
# prepare the dataset
spdata=data.frame(t(ratio))
spdata=cbind(spdata,info)
dummy <- rep(1, dim(spdata)[1]) 
spdata <- cbind(spdata, dummy) 

#Use exponential correlation structure
#First fit a lme function using corExp
#Create list to store updated models
mod.list=list()
#Build model
colnames(spdata)[1]="variable"
variable <- spdata[,1]
mod=lme(fixed = variable ~ HG+UC, data = spdata, random = ~ 1 | dummy,correlation = corExp(1, form = ~ x + y), method = "ML") 

#Create a data list
spdata.input=list()
for(i in 1:dim(ratio)[1]){
  tmp=data.frame(cbind(spdata[,i],spdata[,c("HG","UC","x","y","dummy")]))
  colnames(tmp)[1]="variable"
  spdata.input[[i]]=tmp
}

#Update models using lapply and store in a list
mod.list=lapply(seq_along(spdata.input),function(i) {
  mod2=try(update(mod,data=spdata.input[[i]]),TRUE)
  if(isTRUE(class(mod2)=="try-error")) { return(NULL) } else { return(mod2) } } )

#Build model without intercept
variable <- spdata[,1]
mod=lme(fixed = variable ~ -1+HG+UC, data = spdata, random = ~ 1 | dummy,correlation = corExp(1, form = ~ x + y), method = "ML") 

#Update models using lapply and store in a list
mod.list.int=lapply(seq_along(spdata.input),function(i) {
  mod2=try(update(mod,data=spdata.input[[i]]),TRUE)
  if(isTRUE(class(mod2)=="try-error")) { return(NULL) } else { return(mod2) } } )

#Build model without HG
variable <- spdata[,1]
mod=lme(fixed = variable ~ UC, data = spdata, random = ~ 1 | dummy,correlation = corExp(1, form = ~ x + y), method = "ML") 

#Update models using lapply and store in a list
mod.list.hg=lapply(seq_along(spdata.input),function(i) {
  mod2=try(update(mod,data=spdata.input[[i]]),TRUE)
  if(isTRUE(class(mod2)=="try-error")) { return(NULL) } else { return(mod2) } } )

#Build model without UC
variable <- spdata[,1]
mod=lme(fixed = variable ~ HG, data = spdata, random = ~ 1 | dummy,correlation = corExp(1, form = ~ x + y), method = "ML") 

#Update models using lapply and store in a list
mod.list.uc=lapply(seq_along(spdata.input),function(i) {
  mod2=try(update(mod,data=spdata.input[[i]]),TRUE)
  if(isTRUE(class(mod2)=="try-error")) { return(NULL) } else { return(mod2) } } )

id=which(lengths(mod.list)==0) 
id1=which(lengths(mod.list.int)==0) 
id2=which(lengths(mod.list.hg)==0) 
id3=which(lengths(mod.list.uc)==0) 
#if(length(id)==0){mod.list.complete=mod.list} else {mod.list.complete=mod.list[-id]}
stat.res.nlme=array(0,c(length(mod.list),3))
for(i in 1:length(mod.list)){
  if(is.element(i,id)==TRUE | is.element(i,id1)==TRUE){stat.res.nlme[i,1]=NA} else{
    G2 = 2 * logLik(mod.list[[i]]) - 2 * logLik(mod.list.int[[i]])
    stat.res.nlme[i,1]=pchisq(as.numeric(G2), df=1, lower.tail=F)
  }
  
  if(is.element(i,id)==TRUE | is.element(i,id2)==TRUE){stat.res.nlme[i,2]=NA} else{
    G2 = 2 * logLik(mod.list[[i]]) - 2 * logLik(mod.list.hg[[i]])
    stat.res.nlme[i,2]=pchisq(as.numeric(G2), df=1, lower.tail=F)
  }
  
  if(is.element(i,id)==TRUE | is.element(i,id3)==TRUE){stat.res.nlme[i,3]=NA} else{
    G2 = 2 * logLik(mod.list[[i]]) - 2 * logLik(mod.list.uc[[i]])
    stat.res.nlme[i,3]=pchisq(as.numeric(G2), df=1, lower.tail=F)
  }
}
rownames(stat.res.nlme)=rownames(ratio)
colnames(stat.res.nlme)=c("LG","HG","UC")

#Next extract the fixed effect estimate
#Summarize the results
#check which gene returned error
id=which(lengths(mod.list)==0)  #24
#if(length(id)==0){mod.list.complete=mod.list} else {mod.list.complete=mod.list[-id]}
stat.res.matrix=array(0,c(length(mod.list),15))
for(i in 1:length(mod.list)){
  if(is.element(i,id)==TRUE){stat.res.matrix[i,]=NA} else{
    stat.res=summary(mod.list[[i]])$tTable
    rownames(stat.res)[1]="LG"
    stat.res.matrix[i,]=gdata::unmatrix(stat.res)}
}
colnames(stat.res.matrix)=names(gdata::unmatrix(stat.res))
rownames(stat.res.matrix)=rownames(ratio)

#Next prepare the z score to fit IMIX
#First, remove NA valued genes
stat.res.nlme.complete=stat.res.nlme[complete.cases(stat.res.nlme),] #remove the missing values

#change p-value 0 or 1 to 0.00001 or 0.99999 so that the z scores will not be infinity
stat.res.nlme.complete[stat.res.nlme.complete==1]=0.99999
stat.res.nlme.complete[stat.res.nlme.complete==0]=0.00001
zscores=apply(stat.res.nlme.complete,2,function(x) qnorm(x,lower.tail = F))
res=list(IMIX_Input_Zscores = zscores,
         spatial_mixed_model_result=stat.res.matrix,
         full_model=mod.list,
         model_withouthg=mod.list.hg,
         model_withoutintercept=mod.list.int,
         model_withoutuc=mod.list.uc)
}
