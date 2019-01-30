args <- commandArgs(TRUE)
input <- args[[1]]
ref_path <- args[[2]]
norm_method <- args[[3]]
assembly <- args[[4]]


age_transform <- function(age,cutoff)
{
	age <- age
	cutoff <- cutoff

	age_young <- age[age <= cutoff]
	age_young <- log(age_young+1)-log(cutoff+1)

	age_adult <- age[age > cutoff]
	age_adult <- (age_adult-cutoff)/(cutoff+1)

	age[age <= cutoff] <- age_young
	age[age > cutoff] <- age_adult
	return(age)
}

inv_age_transform <- function(age,cutoff)
{
	age <- age
	cutoff <- cutoff

	age_young_index <- age <= 0
	age_young <- age[age_young_index]
	age_young <- (cutoff+1)*exp(age_young)-1

	age_adult_index <- age > 0
	age_adult <- age[age_adult_index]
	age_adult <- (cutoff+1)*age_adult+cutoff

	age[age_young_index] <- age_young
	age[age_adult_index] <- age_adult
	return(age)	
}

get_beta_grange <- function(beta)
{
	beta <- beta
	id <- rownames(beta)
	pos <- unlist(gregexpr("_",id))
	chr_id <- substr(id,1,pos-1)
	start_pos <- as.numeric(substr(id,pos+1,99))
	chr_id <- chr_id[!is.na(start_pos)]
	start_pos <- start_pos[!is.na(start_pos)]
	results <- GRanges(chr_id,IRanges(start=start_pos,end=start_pos))

	return(results)
}

remove_missing_beta <- function(beta)
{
	beta <- beta
	beta_mean <- apply(beta,1,mean,na.rm=T)
	indc_beta_miss <- which(is.na(beta),arr.ind=T)
	beta[indc_beta_miss] <- beta_mean[indc_beta_miss[,1]]

	return(beta)
}


beta_process <- function(beta,beta_train,beta_train_grange,miss_cutoff,zero_cutoff)
{
	beta <- beta
	beta_train <- beta_train
	beta_train_grange <- beta_train_grange
	miss_cutoff <- miss_cutoff
	zero_cutoff <- zero_cutoff

	num_cpg999 <- apply(beta==999,1,sum)
	beta <- beta[num_cpg999 < ncol(beta)*miss_cutoff,]
	num_cpg0 <- apply(beta==0,1,sum)
	beta <- beta[num_cpg0 < ncol(beta)*zero_cutoff,]

	indc_beta_miss <- which(beta==999,arr.ind=T)
	if(nrow(indc_beta_miss) > 0)
	{
	beta[indc_beta_miss] <- NA
	beta <- remove_missing_beta(beta)
	}

	beta_grange <- get_beta_grange(beta)
	index <- findOverlaps(query=beta_train_grange,subject=beta_grange)
	index <- cbind(queryHits(index),subjectHits(index))
	index2 <- as.numeric(levels(factor(index[,1])))
	index2 <- cbind(index2,index[match(index2,index[,1]),2])

	beta <- beta[index2[,2],]
	beta_train <- beta_train[,index2[,1]]

	results <- vector("list",2)
	results[[1]] <- beta
	results[[2]] <- beta_train
	return(results)
}

simple_norm <- function(beta,beta_train)
{
	beta <- beta
	beta_train <- beta_train

	x <- apply(beta,2,quantile,probs=seq(0,1,by=0.01))
	x_mean <- apply(x,1,mean)
	y <- apply(beta_train,1,quantile,probs=seq(0,1,by=0.01))
	y_mean <- apply(y,1,mean)

	fit <- lm(y_mean ~ poly(x_mean,3,raw=T))
	beta_norm <- beta
	for(i in 1:ncol(beta))
	{beta_norm[,i] <- predict(fit,data.frame(x_mean=beta[,i]))}

	return(beta_norm)
	
}

DNAm_age_predict <- function(beta,beta_train,age)
{
	beta <- beta
	beta_train <- beta_train
	age <- age

	model <- mixed.solve(y=age,Z=beta_train)
	names(model$u) <- colnames(beta_train)
	age_predict <- as.numeric(model$beta)+as.numeric(t(beta)%*%as.numeric(model$u))
	age_predict <- inv_age_transform(age_predict,20)

	results <- vector("list",2)
	results[[1]] <- age_predict
	results[[2]] <- model
	return(results)
}


library(methylKit)
library(rrBLUP)
library(preprocessCore)


cpg_train_anno <- read.csv(paste(ref_path,"/anno_cpg_illumina_array_train_450k.csv",sep=""),row.names=1)

if(assembly=="GRCh37")
{
cpg_train_grange <- GRanges(sub("chr","",cpg_train_anno[,1]),IRanges(start=cpg_train_anno[,2]-1,end=cpg_train_anno[,2]+1))
}

if(assembly=="GRCh38")
{
index_grch38 <- !is.na(cpg_train_anno[,3])
cpg_train_anno <- cpg_train_anno[index_grch38,]
cpg_train_grange <- GRanges(sub("chr","",cpg_train_anno[,3]),IRanges(start=cpg_train_anno[,4],end=cpg_train_anno[,5]))
}


sample <- list.files(input,pattern = ".txt")
loc <- as.list(paste(input,sample,sep=""))
if(length(grep(".txt.gz",sample)) > 0)
{sample.id <- as.list(sub(".txt.gz","",sample))}
if(length(grep(".txt.gz",sample)) == 0)
{sample.id <- as.list(sub(".txt","",sample))}
sampleNum <- length(sample)
caseNum <- round(sampleNum/2)
design <- c(rep(1,caseNum),rep(0,sampleNum-caseNum))


load(paste(ref_path,"/train_data_450k.Rdata",sep=""))
da_train_450k <- results[[1]]
if(assembly=="GRCh38")
{da_train_450k <- da_train_450k[,index_grch38]}
age_train <- results[[2]]
age_train_adj <- age_transform(age_train,20)
age <- read.csv(paste(input,"/age.csv",sep=""),row.names=1)


cpgNum <- c()
age_file <- c()
cpgNum_indc <- 20001
cov_min <- 1

while(cpgNum_indc > 20000 && cov_min <= 10)
{
da_cpg <- methRead(loc,sample.id=sample.id,assembly=assembly,treatment=design,context="CpG",mincov=cov_min)

####################
## get union beta ##
####################

beta <- getData(da_cpg[[1]])
beta[,1] <- sub("chr","",beta[,1])
beta_grange <- GRanges(beta[,1],IRanges(start=beta[,2],end=beta[,2]))
index <- findOverlaps(query=beta_grange,subject=cpg_train_grange)
beta <- beta[queryHits(index),]
beta_grange <- beta_grange[queryHits(index)]

n <- rep(1,nrow(beta))

for(i in 2:length(da_cpg))
{
	temp <- getData(da_cpg[[i]])
	temp[,1] <- sub("chr","",temp[,1])
	temp_grange <- GRanges(temp[,1],IRanges(start=temp[,2],end=temp[,2]))
	index <- findOverlaps(query=temp_grange,subject=cpg_train_grange)
	temp <- temp[queryHits(index),]
	temp_grange <- temp_grange[queryHits(index)]

	index <- findOverlaps(query=beta_grange,subject=temp_grange)

	temp_grange2 <- temp_grange[-subjectHits(index)]

	n[queryHits(index)] <- n[queryHits(index)]+1
	n <- c(n,rep(1,length(temp_grange2)))
	beta_grange <- c(beta_grange,temp_grange2)
}

miss_cutoff <- 0.2
beta_grange <- beta_grange[n >= (1-miss_cutoff)*length(da_cpg)]


beta <- matrix(999,ncol=length(sample),nrow=length(beta_grange))

for(i in 1:length(da_cpg))
{
	temp <- getData(da_cpg[[i]])
	temp[,1] <- sub("chr","",temp[,1])
	temp_grange <- GRanges(temp[,1],IRanges(start=temp[,2],end=temp[,2]))
	index <- findOverlaps(query=beta_grange,subject=temp_grange)

	temp_beta <- temp[subjectHits(index),]
	beta[queryHits(index),i] <- round(temp_beta[,6]/temp_beta[,5],4)
}

meth_id <- paste(seqnames(beta_grange),"_",start(beta_grange),sep="")
rownames(beta) <- meth_id
colnames(beta) <- unlist(sample.id)
#save(beta,file=paste(output,"beta_union_",assembly,"_bismark.Rdata",sep=""))
age <- age[match(colnames(beta),rownames(age)),]



##################
## get DNAm age ##
##################

## 1. process beta
temp <- beta_process(beta,da_train_450k,cpg_train_grange,miss_cutoff=0.2,zero_cutoff=0.5)
beta <- temp[[1]]
da_train <- temp[[2]]
cpgNum_indc <- nrow(beta)

if(cpgNum_indc >= 20000)
{
output <- paste(input,"/age_pred/cov_",cov_min,"/",sep="")
dir.create(output,recursive=T)

cpgNum <- c(cpgNum,nrow(beta))
write.csv(cpgNum_indc,file=paste(output,"/cpgNum.csv",sep=""),row.names=F)

## 2. normalize beta

if(norm_method=="simplified_norm")
{beta <- simple_norm(beta,da_train)}
if(norm_method=="quantile_norm")
{beta <- normalize.quantiles.use.target(beta,apply(da_train,2,mean))}

## 3. prediction

age_pred <- DNAm_age_predict(beta,da_train,age_train_adj)
if(assembly=="GRCh37")
{age_adj <- age_pred[[1]]+age$age*0.51916-0.002074*age$age^2-0.732037}
if(assembly=="GRCh38")
{age_adj <- age_pred[[1]]+age$age*0.2560243+0.0008109*age$age^2+2.4097794}
age_adj[age$age==0] <- age_pred[[1]][age$age==0]

age_file <- cbind(age_file,age_pred[[1]],age_adj)
age_file2 <- cbind(age$age,age_pred[[1]],age_adj)
colnames(age_file2) <- c("Biological age","Methylation age unadj","Methylation age adj")
rownames(age_file2) <- rownames(age)
cpgEffect <- c(age_pred[[2]]$beta,age_pred[[2]]$u)
names(cpgEffect) <- c("beta",names(age_pred[[2]]$u))
write.csv(round(age_file2,2),file=paste(output,"/DNAm_age.csv",sep=""),row.names=T)
write.csv(round(cpgEffect,5),file=paste(output,"/cpgEffect.csv",sep=""),row.names=T)

cov_min <- cov_min+1
}
}

age_unadj <- age_file[,seq(1,ncol(age_file),by=2)] 
age_adj <- age_file[,seq(2,ncol(age_file),by=2)]

cov_indc <- 1:length(cpgNum)

age_unadj_mean <- apply(age_unadj[,cov_indc],1,mean)
age_adj_mean <- apply(age_adj[,cov_indc],1,mean)

age_file <- cbind(age$age,age_unadj,age_adj,age_unadj_mean,age_adj_mean)

colnames(age_file) <- c("Biological age",paste("Methylation age unadj cov",cov_indc,sep="_"),paste("Methylation age adj cov",cov_indc,sep="_"),"Methylation age unadj","Methylation age adj")
rownames(age_file) <- rownames(age)

write.csv(round(age_file,2),file=paste(input,"/age_pred/DNAm_age_all.csv",sep=""),row.names=T)







