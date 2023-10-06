###Covariate selection workflow
###Melanie Garrett July 2023
###This is meant to be a guideline as the process is iterative and specific to your data and situation.

setwd("your_working_dir")

#log2 transformed expression matrix - summed by sample
expr <- read.table("sum_counts/expr_sampsum.txt",header=T,row.names=1)

#remove lowly expressed genes (0 in >80% of samples - expr in 20+ samples)
#n=24 samples
expr$noexpr <- rowSums(expr == "0")
expr.goodexpr <- expr[which(expr$noexpr<4),]
expr.goodexpr <- expr.goodexpr[,c(-25)] #remove rowsum

#scale data by gene (mean 0, sd 1)
expr.t <- t(expr.goodexpr)
expr.scaled <- scale(expr.t)

#pca
pca.expr <- prcomp((expr.scaled),scale=T,retx=T,center=T)

#calculate proportion of variance explained by each PC
summary(pca.expr)
ev <- pca.expr$sdev^2
prop.var <- ev/sum(ev)
png(file="prop_var_explained.png")
par(mar=c(4,4,2,2))
barplot(prop.var,ylim=c(0,0.75),xlab="PCs",ylab="proportion of variance explained")
dev.off()

pcs <- pca.expr$x
pcs <- pcs[,c(1:20),drop=FALSE] #keep only first 20 PCs

#metadata - file containing technical variables to test for association with expression data
meta <- read.csv("metadata_scrnaseq.csv",header=T,row.names=1)
all(rownames(meta)==rownames(pcs))
all <- cbind(pcs,meta)

all.new <- all[,!(names(all) %in% c("SampleID"))]

#association of metadata with first PC1
pc1 <- all$PC1
start=21 #this is where the metadata variables begin
end=ncol(all.new)
nvar=end-start+1

var=rep(NA,nvar)
beta=rep(NA,nvar)
se=rep(NA,nvar)
pval=rep(NA,nvar)

number=1

for (i in start:end){
        vars=colnames(all.new)[i]
        model <- glm(pc1 ~ get(vars), family="gaussian", data=all.new)
        beta[number]=as.numeric(coef(summary(model))[2,1])
        se[number]=as.numeric(coef(summary(model))[2,2])
        pval[number]=as.numeric(coef(summary(model))[2,4])
        var[number]=vars

        number=number+1
        }

results <- data.frame(var,beta,se,pval)
results.pc1 <- results[,c(1,4)]

rownames(results.pc1) <- results.pc1$var
results.pc1$neglog10p <- -(log10(results.pc1$pval))
res2 <- results.pc1[,c(3),drop=FALSE]
res3 <- res2[order(res2$neglog10p),,drop=FALSE]

png(file="expr_pc1_assoc1.png",units='px',width=3000,height=3000,res=300)
par(mar=c(4,12,2,2))
barplot(t(res3),horiz=T,xlim=c(0,12),xlab="-log10(p)",las=1,cex.names=0.75)
abline(v=1.3) #nominal association (corresponding to p=0.05)
abline(v=2.88) #Bonferroni significance threshold (p=0.05/# of vars tested)
dev.off()

#If no variables are associated with PC1 (at Bonferroni threshold), then run the code above again for PC2, then PC3, etc.  Continue until proportion of variance explained for the PC in question is below 5%.

#If variables are associated with PC1, use the code below to regress out the effect of that variable from all expression data and repeat the process until no variables are associated with PCs explaining more than 5% of the data.

#regress out most significant variable that meets Bonferroni threshold
all(rownames(meta)==rownames(expr.scaled))
round1 <- cbind(meta,expr.scaled)

associated.var <- as.numeric(round1[,"seq_saturation"])
start=40 #this is where the expression data starts
end=ncol(round1)
nvar=end-start+1

resid=matrix(rep(NA,nvar),ncol=nvar,nrow=24) #24 samples

number=1

for (i in start:end){
        outcome <- as.numeric(round1[,(colnames(round1)[i])])
        model <- glm(outcome ~ associated.var, family="gaussian")
        resid[,number]=residuals(model)

        number=number+1
        }

rownames(resid) <- rownames(expr.scaled)
colnames(resid) <- colnames(expr.scaled)

#now run PCA on residuals
pca.resid <- prcomp(resid,scale=T,retx=T,center=T)

#calculate proportion of variance explained by each PC
summary(pca.resid)
ev <- pca.resid$sdev^2
prop.var <- ev/sum(ev)
png(file="prop_var_explained_resid.png")
barplot(prop.var,ylim=c(0,0.5),xlab="PCs",main="proportion of variance explained")
dev.off()

pcs <- pca.resid$x
pcs <- pcs[,c(1:20),drop=FALSE]

all(rownames(meta)==rownames(pcs))
all <- cbind(pcs,meta)
all.new <- all[,!(names(all) %in% c("SampleID","associated.var"))] #remove these variables because they don't need to be tested

#association of metadata with first PC
pc1 <- all.new$PC1
start=21 #this is where the meta data starts
end=ncol(all.new)
nvar=end-start+1

var=rep(NA,nvar)
beta=rep(NA,nvar)
se=rep(NA,nvar)
pval=rep(NA,nvar)

number=1

for (i in start:end){
        vars=colnames(all.new)[i]
        model <- glm(pc1 ~ get(vars), family="gaussian", data=all.new)
        beta[number]=as.numeric(coef(summary(model))[2,1])
        se[number]=as.numeric(coef(summary(model))[2,2])
        pval[number]=as.numeric(coef(summary(model))[2,4])
        var[number]=vars

        number=number+1
        }

results <- data.frame(var,beta,se,pval)
results.pc1 <- results[,c(1,4)]

rownames(results.pc1) <- results.pc1$var
results.pc1$neglog10p <- -(log10(results.pc1$pval))
res2 <- results.pc1[,c(3),drop=FALSE]
res3 <- res2[order(res2$neglog10p),,drop=FALSE]

png(file="expr_pc1_assoc2.png",units='px',width=3000,height=3000,res=300)
par(mar=c(4,12,2,2))
barplot(t(res3),horiz=T,xlim=c(0,4),xlab="-log10(p)",las=1,cex.names=0.75)
abline(v=1.3)
abline(v=2.87) #Bonferroni threshold changes because there is 1 less variable tested
dev.off()

#If additional variables are associated with PC1, repeat the above process to regress them out of all expression data, run PCA again on those residuals, and test remaining meta variables until no more are associated.

#If additional variables are NOT associated with PC1, move on to test PC2, PC3, etc. until you've tested all PCs that explain more than 5% of the variability in the data.  If you find a variable associated, repeat the above process to regress them out of all expression data, run PCA again on those residuals, and test remaining meta variables until no more are associated.

#assoc with PC2
pc2 <- all.new$PC2
start=21
end=ncol(all.new)
nvar=end-start+1

var=rep(NA,nvar)
beta=rep(NA,nvar)
se=rep(NA,nvar)
pval=rep(NA,nvar)

number=1

for (i in start:end){
        vars=colnames(all.new)[i]
        model <- glm(pc2 ~ get(vars), family="gaussian", data=all.new)
        beta[number]=as.numeric(coef(summary(model))[2,1])
        se[number]=as.numeric(coef(summary(model))[2,2])
        pval[number]=as.numeric(coef(summary(model))[2,4])
        var[number]=vars

        number=number+1
        }

results <- data.frame(var,beta,se,pval)
results.pc2 <- results[,c(1,4)]

rownames(results.pc2) <- results.pc2$var
results.pc2$neglog10p <- -(log10(results.pc2$pval))
res2 <- results.pc2[,c(3),drop=FALSE]
res3 <- res2[order(res2$neglog10p),,drop=FALSE]

png(file="expr_pc2_assoc2.png",units='px',width=3000,height=3000,res=300)
par(mar=c(4,12,2,2))
barplot(t(res3),horiz=T,xlim=c(0,3),xlab="-log10(p)",las=1,cex.names=0.75)
abline(v=1.3)
abline(v=2.87)
dev.off()

#assoc with PC3
pc3 <- all.new$PC3
start=21
end=ncol(all.new)
nvar=end-start+1

var=rep(NA,nvar)
beta=rep(NA,nvar)
se=rep(NA,nvar)
pval=rep(NA,nvar)

number=1

for (i in start:end){
        vars=colnames(all.new)[i]
        model <- glm(pc3 ~ get(vars), family="gaussian", data=all.new)
        beta[number]=as.numeric(coef(summary(model))[2,1])
        se[number]=as.numeric(coef(summary(model))[2,2])
        pval[number]=as.numeric(coef(summary(model))[2,4])
        var[number]=vars

        number=number+1
        }

results <- data.frame(var,beta,se,pval)
results.pc3 <- results[,c(1,4)]

rownames(results.pc3) <- results.pc3$var
results.pc3$neglog10p <- -(log10(results.pc3$pval))
res2 <- results.pc3[,c(3),drop=FALSE]
res3 <- res2[order(res2$neglog10p),,drop=FALSE]

png(file="expr_pc3_assoc2.png",units='px',width=3000,height=3000,res=300)
par(mar=c(4,12,2,2))
barplot(t(res3),horiz=T,xlim=c(0,3),xlab="-log10(p)",las=1,cex.names=0.75)
abline(v=1.3)
abline(v=2.87)
dev.off()

#assoc with PC4
pc4 <- all.new$PC4
start=21
end=ncol(all.new)
nvar=end-start+1

var=rep(NA,nvar)
beta=rep(NA,nvar)
se=rep(NA,nvar)
pval=rep(NA,nvar)

number=1

for (i in start:end){
        vars=colnames(all.new)[i]
        model <- glm(pc4 ~ get(vars), family="gaussian", data=all.new)
        beta[number]=as.numeric(coef(summary(model))[2,1])
        se[number]=as.numeric(coef(summary(model))[2,2])
        pval[number]=as.numeric(coef(summary(model))[2,4])
        var[number]=vars

        number=number+1
        }

results <- data.frame(var,beta,se,pval)
results.pc4 <- results[,c(1,4)]

rownames(results.pc4) <- results.pc4$var
results.pc4$neglog10p <- -(log10(results.pc4$pval))
res2 <- results.pc4[,c(3),drop=FALSE]
res3 <- res2[order(res2$neglog10p),,drop=FALSE]

png(file="expr_pc4_assoc2.png",units='px',width=3000,height=3000,res=300)
par(mar=c(4,12,2,2))
barplot(t(res3),horiz=T,xlim=c(0,3),xlab="-log10(p)",las=1,cex.names=0.75)
abline(v=1.3)
abline(v=2.87)
dev.off()

#If no meta variables are associated with any PCs explaining more than 5% of the variability in the expression data, you're done!  
# All of the variables that were regressed out are those that are important to include in your differential expression models as covariates.
