args <- commandArgs(TRUE)
input <- args[[1]]
ref_path <- args[[2]]


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


library(rrBLUP)
library(preprocessCore)

file <- list.files(input,pattern="beta")
load(paste(input,"/",file,sep=""))
load(paste(ref_path,"train_data_model.Rdata",sep=""))

beta_mean <- apply(beta,1,mean,na.rm=T)
indc_beta_miss <- which(is.na(beta),arr.ind=T)
beta[indc_beta_miss] <- beta_mean[indc_beta_miss[,1]]

beta <- t(beta)

index <- match(colnames(results[[1]]),colnames(beta))
if(sum(is.na(index))==0)
{
	age <- as.numeric(beta[,index]%*%results[[3]]$u)+as.numeric(results[[3]]$beta)
	age <- inv_age_transform(age,20)
}

if(sum(is.na(index)) > 0)
{
	age_train <- results[[2]]
	age_train_adj <- age_transform(age_train,20)
	
	da_train <- results[[1]][,!is.na(index)]
	model <- mixed.solve(y=age_train_adj,Z=da_train)
	
	age <- as.numeric(beta[,index[!is.na(index)]]%*%model$u)+as.numeric(model$beta)
	age <- inv_age_transform(age,20)
}

names(age) <- rownames(beta)
write.csv(round(age,2),file=paste(input,"/DNAm_age.csv",sep=""),row.names=T)








