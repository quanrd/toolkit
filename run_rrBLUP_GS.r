library(rrBLUP)
library(ggplot2)

args <- commandArgs(TRUE)
geno.file <- args[1]
pheno.file <- args[2]
out.name <- args[3]

# read datas
geno.all <- read.delim(geno.file)
pheno.all <- read.delim(pheno.file)

# impute missing markers
all.impute <- A.mat(geno.all,max.missing=0.05,impute.method="mean",return.imputed=T)
geno.all.impute <- all.impute$imputed

pheno.names <- colnames(pheno.all)
sample.num <- dim(pheno.all)[1]
pheno.num <- dim(pheno.all)[2]

# set train and valid dataset
train.rank <- as.matrix(sample(1:sample.num,sample.num/2))
valid.rank <- setdiff(1:sample.num,train)
geno.train <- geno.all.impute[train.rank,]
geno.valid <- geno.all.impute[valid.rank,]
pheno.train <- pheno.all[train.rank,]
pheno.valid <- pheno.all[valid.rank,]

train.sample.num <- length(train.rank)
train.out.m <- matrix(nrow=train.sample.num,ncol=pheno.num-1)
rownames(train.out.m) <- as.character(pheno.train[,1])
colnames(train.out.m) <- pheno.names[2:pheno.num]
train.accuracy.m <- matrix(nrow=pheno.num-1)
rownames(train.accuracy.m) <- colnames(pheno.train[,2:pheno.num])
colnames(train.accuracy.m) <- c("accuracy")

valid.sample.num <- length(valid.rank)
valid.out.m <- matrix(nrow=valid.sample.num,ncol=pheno.num-1)
rownames(valid.out.m) <- as.character(pheno.valid[,1])
colnames(valid.out.m) <- pheno.names[2:pheno.num]
valid.accuracy.m <- matrix(nrow=pheno.num-1)
rownames(valid.accuracy.m) <- colnames(pheno.valid[,2:pheno.num])
colnames(valid.accuracy.m) <- c("accuracy")

for(i in 2:pheno.num){
	
	pheno.name <- pheno.names[i]
	
	# training
	train.out <- mixed.solve(pheno.train[,i],Z=geno.train,K=NULL,SE=FALSE,return.Hinv=F)

	# validation
	train.out.u <- train.out$u
	e <- as.matrix(train.out.u)
	
	train.pred <- geno.train %*% e
	train.pred.value <- (train.pred[,1]) + train.out$beta
	train.accuracy <- cor(train.pred.value,pheno.train[,i],use="complete")
	train.out.m[,i-1] <- train.pred.value
	train.accuracy.m[i-1,1] <- train.accuracy
	
	valid.pred <- geno.valid %*% e
	valid.pred.value <- (valid.pred[,1]) + train.out$beta
	valid.accuracy <- cor(valid.pred.value,pheno.valid[,i],use="complete")
	valid.out.m[,i-1] <- valid.pred.value
	valid.accuracy.m[i-1,1] <- valid.accuracy
	
	out.png <- paste(out.name,pheno.name,"train.png",sep=".")
	png(out.png)
	print(ggplot() + geom_point(aes(x=train.pred.value,y=pheno.train[,i])) + geom_abline(slope=1) + xlab("prediction") + ylab("real") + ggtitle(pheno.name))
	dev.off()
	
	out.png <- paste(out.name,pheno.name,"valid.png",sep=".")
	png(out.png)
	print(ggplot() + geom_point(aes(x=valid.pred.value,y=pheno.valid[,i])) + geom_abline(slope=1) + xlab("prediction") + ylab("real") + ggtitle(pheno.name))
	dev.off()
}

out.acc.file <- paste(out.name,"accuracy.valid.txt",sep=".")
write.table(valid.accuracy.m,file=out.acc.file,sep="\t",quote=F)
out.acc.file <- paste(out.name,"accuracy.train.txt",sep=".")
write.table(train.accuracy.m,file=out.acc.file,sep="\t",quote=F)

out.value.file <- paste(out.name,"value.valid.txt",sep=".")
write.table(valid.out.m,file=out.value.file,sep="\t",quote=F)
out.value.file <- paste(out.name,"value.train.txt",sep=".")
write.table(train.out.m,file=out.value.file,sep="\t",quote=F)