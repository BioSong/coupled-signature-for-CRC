####data preprocessing 
##probe files 2 gene profiles
library(Biobase) #install.packages("Biobase")
library(affy)
setwd(path)
dir1="affy//" #
len=length(disease.file)
#### RMA -------------------------------------------------------------------------------------------------
for(i in 1:len){
print(i)
OutputPath=paste("affy//",disease.file[i],"//",sep="")
InputPath=paste(OutputPath,"all//",sep="")
setwd(InputPath)                 
ExpData <- justRMA(normalized=F)      #对于特定一套数据标准化   
setwd(OutputPath)
write.exprs(ExpData,file=disease.file[i]) #write.exprs将expressionSet写出 这里将探针表达谱输出
}
rm(ExpData)

#### probid to geneid function-----------------------------------------------------------------------------
pid2gid = function(exp, probeid_geneid) 
{
	exp=as.matrix(exp)
    index = which(rowSums(is.na(exp)) == 0)
    subexp = exp[index, ]   
    probeid_geneid = unique(probeid_geneid)
    index1 = which(rowSums(is.na(probeid_geneid)) == 0)
    index2 = which(!is.na(as.numeric(as.character(probeid_geneid[, 2]))))
    id_index = intersect(index1, index2)
    
    pid_gid = probeid_geneid[id_index, ]
    pid = intersect(subexp[, 1], pid_gid[, 1])
    index2 = match(pid, subexp[, 1])
    subexp1 = subexp[index2, -1]
    mode(subexp1) = "numeric"
    index3 = match(pid, pid_gid[, 1])
    ### mean
    subgeneid = factor(pid_gid[index3, 2])
    geneexp = apply(subexp1, 2, function(x) tapply(x, subgeneid, mean))
    return(geneexp)
} 
##-------------------
dir1=dir("pro_exp//")
for (file1 in dir1){
	dir2=paste("platform//","processed//",file1,".txt",sep="")
	probeid_geneid=read.table(dir2,sep="\t",header=T,as.is=T,fill=T,quote="\"")
	dir3=paste("pro_exp//",file1,"//",sep="")
	dir4="gene_exp//"
	file3=dir(dir3)
	for(file2 in file3){
		proexp=read.table(paste(dir3,file2,sep=""),header=F,as.is=T,quote="\"",fill=T)
		proexp=as.matrix(proexp)
		proexp=proexp[which(proexp[,1]=="ID_REF"):nrow(proexp),]
		colnames(proexp)=proexp[1,]
		proexp=proexp[-1,]
		geneexp=pid2gid(proexp,probeid_geneid)
		write.table(geneexp,paste(dir4,file2,sep=""),sep="\t",col.names=T,row.names=T,quote=F)
	}
}
#-------------train relapse risk signature
#TCGA-------------DEGs between metastasis and non-metastasis
stage=read.table("data for metastasis genes//stage.txt",sep="\t",header=T,as.is=T,quote="\"")
stage1=stage[stage[,2]==1,1]
stage4=stage[stage[,2]==4,1]
#############3
count=read.table("data for metastasis genes//v2count_cancer.txt",sep="\t",header=T,row.names=1,quote="\"",as.is=T)
sams2=substr(names(count),1,1)
sams2=gsub("\\.","-",sams2)
sam1=intersect(sams2,stage1)
sam4=intersect(sams2,stage4)
can1=count[,match(sam1,sams2)]
can4=count[,match(sam4,sams2)]
exp1=cbind(can4,can1)
geneid=as.numeric(rownames(count))
group=factor(c(rep("meta",ncol(can4)),rep("nometa",ncol(can1))))
data1=DGEList(counts=exp1,genes=geneid,group=group)
index1=rowSums(cpm(data1)>1)>=(ncol(exp1)/2)
data1=data1[index1,]
data1=calcNormFactors(data1)
data1=estimateCommonDisp(data1)
data1=estimateTagwiseDisp(data1)
et1=exactTest(data1)
fdr=p.adjust(et1$table[,3],method="BH")
et1$table[,1]=-et1$table[,1]
write.table(cbind(et1$genes,et1$table,fdr),"tcga_deg_I_IV.txt",sep="\t",col.names=T,row.names=F,quote=F)
rm(list=ls())
gc()
#--------------------------metastasis-related gene pairs
gid=read.table("gid_used.txt",sep="\t",header=F,as.is=T)
gid=gid[,1]
tresult=read.table("tcga_deg_I_IV.txt",sep="\t",header=T,as.is=T)
mgene=tresult[tresult[,5]<0.1,1]
mgene=intersect(mgene,gid)
mpair=t(combn(mgene,2))
#TCGA GSE39582
mresult=list()
file1=c("TCGA","GSE39582")
for(i in 1:2){
	exp1=read.table(paste("gene_exp//",file1[i],".txt",sep=""),sep="\t",header=T,row.names=1,quote="\"",as.is=T)
	gid1=as.numeric(rownames(exp1))
	sam1=colnames(exp1)
	exp1=as.matrix(exp1)
	cli1=read.table(paste("clinical//",file1[i],".txt",sep=""),sep="\t",header=T,quote="\"",as.is=T)
	sam11=cli1[cli1$stage==1,1]
	sam14=cli1[cli1$stage==4,1]
	exp11=exp1[,na.omit(match(sam11,sam1))]
	com11=exp11[match(mpair[,1],gid1),]-exp11[match(mpair[,2],gid1),]
	index11=rowSums(com11>0)
	exp14=exp1[,na.omit(match(sam14,sam1))]
	com14=exp14[match(mpair[,1],gid1),]-exp14[match(mpair[,2],gid1),]
	index14=rowSums(com14>0)
	mresult1=c()
	for(j in 1:nrow(mpair)){
		ftest=fisher.test(matrix(c(index14[j],ncol(com14)-index14[j],index11[j],ncol(com11)-index11[j]),2,byrow=T))
		minus=(index14[j]/ncol(com14))-(index11[j]/ncol(com11))
		mresult1=rbind(mresult1,c(mpair[j,],minus,ftest$estimate,ftest$p.value))
	}
	mresult1=cbind(mresult1,p.adjust(mresult1[,5],method="BH"))
	mresult[[i]]=mresult1
}
save(mresult,file="mresult_tcga&39582.RData")
rm(list=ls())
gc()
########################
load("mresult_tcga&39582.RData")
mresult1=mresult[[1]]
mpair11=mresult1[mresult1[,4]>2&mresult1[,5]<0.1,c(1,2,3)]
mpair12=mresult1[mresult1[,4]<0.5&mresult1[,5]<0.1,c(2,1,3)]
mpair1=rbind(mpair11,mpair12)
mresult3=mresult[[2]]
mpair31=mresult3[mresult3[,4]>2&mresult3[,5]<0.1,c(1,2,3)]
mpair32=mresult3[mresult3[,4]<0.5&mresult3[,5]<0.1,c(2,1,3)]
mpair3=rbind(mpair31,mpair32)
#intersect gene pairs
colnames(mpair1)=c("a","b","c")
colnames(mpair3)=c("a","b","c")
interp=merge(mpair1,mpair3,by=c("a","b"))
or1=rowMeans(interp[,3:4])
mpair13=cbind(interp[,1:2],or1)
####redundancy removal process
minus=abs(mpair13[,3])
pair31=mpair13[order(minus,decreasing=T),1:2]
pair32=pair31
for(i in 1:nrow(pair31)){
	for(j in 1:2){
		g1=pair31[i,j]
		index1=which(pair32[,1]==g1|pair32[,2]==g1)
		if(length(index1)<2){
			next}
		pair32=pair32[-index1[-1],]
	}
}
#------------------
##################################################################################################################
load("data_for_44gps.RData")
interpair1=pair32
exp39582=geneexp[[1]]
exp39582=as.matrix(exp39582)
cli39582=clinical[[1]]
gid1=as.numeric(rownames(exp39582))
coms=exp39582[match(interpair1[,1],gid1),]-exp39582[match(interpair1[,2],gid1),]
coms[coms>0]=1
coms[coms<0]=0
freq=rowSums(coms)/ncol(coms)
library(survival)
survi=Surv(cli39582$survival,cli39582$survival.event)
cresult=matrix(,nrow(interpair1),3)
for (i in 1:nrow(interpair1)){
	label=coms[i,]
	cox=coxph(survi~label)
	b=summary(cox)
	hr=b$conf.int[1]
	p=b$sctest[3]
	cresult[i,1:2]=c(hr,p)
}
cresult[,3]=p.adjust(cresult[,2],method="BH")
propair1=interpair1[which(cresult[,3]<0.1),]
write.table(propair1,"GPS_44.txt")
#########
cpair1=read.table("GPS_44.txt")
file1=c("GSE39582","GSE92921","GSE30378","GSE14333")
par(mfrow=c(2,3))
#####discovery cohort
i=1
data1=geneexp[[i]]
cli1=clinical[[i]]
gid1=as.numeric(rownames(data1))
data1=as.matrix(data1)
data1=data1[,cli1$stage==2|cli1$stage==3]
cli1=cli1[cli1$stage==2|cli1$stage==3,]
com1=data1[match(cpair1[,1],gid1),]-data1[match(cpair1[,2],gid1),]
label1=colSums(com1>0)
label=c()
label[label1>s2]=1
label[label1<=s2]=0
cli=rbind(cli,cli1[,c("sample","survival","survival.event","stage")])
com=cbind(com,com1)
survi=Surv(as.numeric(cli1$survival),as.numeric(cli1$survival.event))
fit<-survfit(survi~label)
cox=coxph(survi~label)
b=summary(cox)
hr=b$conf.int[1]
lower.hr=b$conf.int[3]
upper.hr=b$conf.int[4]
p=b$sctest[3]
plot(fit,col =c("blue","red"),main=paste(file1[i],"II-III noCTX",sep=" "),ylab = "survival rate",xlab = "survival months",lty=c(1,2),lwd=1.5,mark.time=T)
legend("bottomright",pch=15:18,lty=c(1,2),merge=F,legend=c(paste("low-risk=",length(which(label==0)),sep=""),paste("high-risk=",length(which(label==1)),sep="")),col=c("blue","red"))
text(x=40,y=0.1,paste("log-rank P=",signif(p, digits = 3),"\n","HR=",signif(hr,4),"(","95%CI,",signif(lower.hr,4),"-",signif(upper.hr,4),")"),bty="n",font=2)
##ROC
library(survivalROC)
score=colMeans(com1>0)
fit <- survivalROC(Stime = cli1$survival, status = cli1$survival.event, predict.time = 40, marker=score,method = "KM")
ROC.1<-fit
plot(ROC.1$FP,ROC.1$TP,type='o',pch=16,xlim=c(0,1),ylim=c(0,1),xlab='FP',ylab='TP',main=file1[i],col='red')
abline(0,1)
text(x=0.5,y=0.2,paste("AUC=",signif(fit$AUC,4),seq=""),bty="n",font=2)
#######validation cohorts  Separately
cli=c()
com=c()
labels=c()
for(i in 2:4){
	data1=geneexp[[i]]
	cli1=clinical[[i]]
	gid1=as.numeric(rownames(data1))
	data1=as.matrix(data1)
	data1=data1[,cli1$stage==2|cli1$stage==3]
	cli1=cli1[cli1$stage==2|cli1$stage==3,]
	com1=data1[match(cpair1[,1],gid1),]-data1[match(cpair1[,2],gid1),]
	label1=colSums(com1>0)
	label=c()
	label[label1>s2]=1
	label[label1<=s2]=0
	##combine
	labels=c(labels,label)
	colnames(cli1)[1]="sample"
	cli=rbind(cli,cli1[,c("sample","survival","survival.event","stage")])
	com=cbind(com,com1)
	################
	survi=Surv(as.numeric(cli1$survival),as.numeric(cli1$survival.event))
	fit<-survfit(survi~label)
	cox=coxph(survi~label)
	b=summary(cox)
	hr=b$conf.int[1]
	lower.hr=b$conf.int[3]
	upper.hr=b$conf.int[4]
	p=b$sctest[3]
	plot(fit,col =c("blue","red"),main=paste(file1[i],"II-III noCTX",sep=" "),ylab = "survival rate",xlab = "survival months",lty=c(1,2),lwd=1.5,mark.time=T)
	legend("bottomright",pch=15:18,lty=c(1,2),merge=F,legend=c(paste("low-risk=",length(which(label==0)),sep=""),paste("high-risk=",length(which(label==1)),sep="")),col=c("blue","red"))
	text(x=40,y=0.1,paste("log-rank P=",signif(p, digits = 3),"\n","HR=",signif(hr,4),"(","95%CI,",signif(lower.hr,4),"-",signif(upper.hr,4),")"),bty="n",font=2)
}
#####validation cohort combined
survi=Surv(as.numeric(cli$survival),as.numeric(cli$survival.event))
fit<-survfit(survi~labels)
cox=coxph(survi~labels)
b=summary(cox)
hr=b$conf.int[1]
lower.hr=b$conf.int[3]
upper.hr=b$conf.int[4]
p=b$sctest[3]
plot(fit,col =c("blue","red"),main="combined validation",ylab = "survival rate",xlab = "survival months",lty=c(1,2),lwd=1.5,mark.time=T)
legend("bottomright",pch=15:18,lty=c(1,2),merge=F,legend=c(paste("low-risk=",length(which(label==0)),sep=""),paste("high-risk=",length(which(label==1)),sep="")),col=c("blue","red"))
text(x=40,y=0.1,paste("log-rank P=",signif(p, digits = 3),"\n","HR=",signif(hr,4),"(","95%CI,",signif(lower.hr,4),"-",signif(upper.hr,4),")"),bty="n",font=2)
#####ROC
score=colMeans(com>0)
fit <- survivalROC(Stime = cli$survival, status = cli$survival.event, predict.time = 40, marker=score,method = "KM")
ROC.1<-fit
plot(ROC.1$FP,ROC.1$TP,type='o',pch=16,xlim=c(0,1),ylim=c(0,1),xlab='FP',ylab='TP',main="combined validation",col='red')
abline(0,1)
text(x=0.5,y=0.2,paste("AUC=",signif(fit$AUC,4),seq=""),bty="n",font=2)
####
#分别分类II期和III期样本
file1=c("GSE39582 II","GSE39582 III","GSE14333 II","GSE14333 III","GSE30378 II","GSE30378 III")
par(mfrow=c(2,3))
slabel=list()
k=1
for(i in 2:4){
	data1=geneexp[[i]]
	cli1=clinical[[i]]
	gid1=as.numeric(rownames(data1))
	data1=as.matrix(data1)
	for(j in 2:3){
		data11=data1[,cli1$stage==j]
		cli11=cli1[cli1$stage==j,]
		com1=data11[match(cpair1[,1],gid1),]-data11[match(cpair1[,2],gid1),]
		label1=colSums(com1<0)
		label=c()
		label[label1>s2]=1
		label[label1<=s2]=0
		slabel[[i]]=label
		survi=Surv(as.numeric(cli11$survival),as.numeric(cli11$survival.event))
		fit<-survfit(survi~label)
		cox=coxph(survi~label)
		b=summary(cox)
		hr=b$conf.int[1]
		lower.hr=b$conf.int[3]
		upper.hr=b$conf.int[4]
		p=b$sctest[3]
		plot(fit,col =c("blue","red"),ylab = "survival rate",xlab = "survival months",lty=c(1,2),lwd=1.5,mark.time=T)
		legend("bottomright",pch=15:18,lty=c(1,2),merge=F,legend=c(paste("low-risk=",length(which(label==0)),sep=""),paste("high-risk=",length(which(label==1)),sep="")),col=c("blue","red"))
		text(x=40,y=0.1,paste("log-rank P=",signif(p, digits = 3),"\n","HR=",signif(hr,4),"(","95%CI,",signif(lower.hr,4),"-",signif(upper.hr,4),")"),bty="n",font=2)
		#k=k+1
	}
}
#TCGA 高低风险组和分期的关系
texp=read.table("gene_exp//TCGA.txt",sep="\t",header=T,row.names=1,quote="\"",as.is=T)
gid=as.numeric(rownames(texp))
sam1=colnames(texp)
tcli=read.table("clinical//TCGA.txt",sep="\t",header=T,row.names=1,quote="\"",as.is=T)
sam2=rownames(tcli)
sam12=intersect(sam1,sam2)
exp1=texp[,match(sam12,sam1)]
cli1=tcli[match(sam12,sam2),]
com1=exp1[match(cpair1[,1],gid),]-exp1[match(cpair1[,2],gid),]
label1=colSums(com1>0)
label=c()
label[label1>23]=0
label[label1<=23]=1
tongji=matrix(,4,2)
for(i in 1:4){
	tongji[i,1]=sum(label==0&cli1$stage==i)
	tongji[i,2]=sum(label==1&cli1$stage==i)
}
#分类II-III期用药样本--44GPS
load("data_for_4GPS.RData")
file1=c("GSE39582","GSE14333","GSE87211")
par(mfrow=c(2,3))
plabel=list()
for(i in c(2,1,3)){
	data1=datas[[i]]
	cli1=clis[[i]]
	gid1=as.numeric(rownames(data1))
	data1=as.matrix(data1)
	com1=data1[match(cpair1[,1],gid1),]-data1[match(cpair1[,2],gid1),]
	label1=colSums(com1>0)
	label=c()
	label[label1>23]=1
	label[label1<=23]=0
	#plabel[[i]]=label
	survi=Surv(as.numeric(cli1$survival),as.numeric(cli1$survival.event))
	fit<-survfit(survi~label)
	cox=coxph(survi~label)
	b=summary(cox)
	hr=b$conf.int[1]
	lower.hr=b$conf.int[3]
	upper.hr=b$conf.int[4]
	p=b$sctest[3]
	plot(fit,col =c("blue","red"),main=file1[i],ylab = "survival rate",xlab = "survival months",lty=c(1,2),lwd=1,mark.time=T)
	legend("bottomright",pch=15:18,lty=c(1,2),merge=F,legend=c(paste("low-risk=",length(which(label==0)),sep=""),paste("high-risk=",length(which(label==1)),sep="")),col=c("blue","red"))
	text(x=40,y=0.1,paste("log-rank P=",signif(p, digits = 3),"\n","HR=",signif(hr,4),"(","95%CI,",signif(lower.hr,4),"-",signif(upper.hr,4),")"),bty="n",font=2)
}
rm(list=ls())
gc()
###########train 5-fu signature
load("cell_pair.RData")
load("high_risk_ACT.RData")
exp1=hexp39582
cli1=hcli39582
library(survival)
survi=Surv(cli1$survival,cli1$survival.event)
sindex=c()
coms=exp1[match(cellpair[,1],gid1),]-exp1[match(cellpair[,2],gid1),]
coms[coms>0]=1
coms[coms<0]=0
for (m in 1:nrow(cellpair)){
	label=coms[m,]
	cox=coxph(survi~label)
	b=summary(cox)
	hr=b$conf.int[1]
	p=b$sctest[3]
	cindex=b$concordance[1]
	sindex=c(sindex,cindex)
}
cellpair[,3]=sindex
cellpair=cellpair[order(sindex,decreasing=T),]
seedpair=cellpair[1,,drop=F]
cellpair1=cellpair[2:nrow(cellpair),]
seed_cindex=seedpair[,3]
seed_cindex1=0.5
while(seed_cindex1<seed_cindex){
	cindex=c()
	for (i in 1:nrow(sigpair1)){
		pair1=rbind(seedpair,sigpair1[i,,drop=F])
		coms=exp1[match(pair1[,1],gid1),]-exp1[match(pair1[,2],gid1),]
		label1=colSums(coms>0)
		label=c()
		label[label1<=(nrow(pair1)/2)]=0
		label[label1>(nrow(pair1)/2)]=1
		cox=coxph(survi~label)
		b=summary(cox)
		cindex1=b$concordance[1]
		cindex=c(cindex,cindex1)
	}
	if(max(cindex)>seed_cindex){
		seed_cindex1=seed_cindex
		seed_cindex=max(cindex)
		index1=which(cindex==max(cindex))[1]
		seedpair=rbind(seedpair,sigpair1[index1,,drop=F])
		sigpair1=sigpair1[-index1,]
	}else{break}
}
write.table(sigpair1,"GPS_4.txt")
rm(list=ls())
gc()
################
cpair1=read.table("GPS_4.txt")
hexps=list(hexp39582,hexp14333,hexp87211)
hclis=list(hcli39582,hcli14333,hcli87211)
file1=c("GSE39582","GSE14333","GSE87211")
par(mfrow=c(2,3))
library(survival)
######discovery cohort
data1=hexps[[i]]
cli11=hclis[[i]]
colnames(cli11)[1]="sample"
gid1=as.numeric(rownames(data1))
data1=as.matrix(data1)
com11=data1[match(cpair1[,1],gid1),]-data1[match(cpair1[,2],gid1),]
label1=colSums(com11>0)
label=c()
label[label1>1]=1
label[label1<=1]=0
survi=Surv(cli11$survival,cli11$survival.event)
cli=rbind(cli,cli11[,c("sample","survival","survival.event")])
labels=c(labels,label)
fit<-survfit(survi~label)
cox=coxph(survi~label)
b=summary(cox)
hr=b$conf.int[1]
lower.hr=b$conf.int[3]
upper.hr=b$conf.int[4]
p=b$sctest[3]
plot(fit,col =c("blue","red"),main=paste(file1[i],"H-meta CTX",sep=" "),ylab = "survival rate",xlab = "survival months",lty=c(1,2),lwd=1.5,mark.time=T)
legend("bottomright",pch=15:18,lty=c(1,2),merge=F,legend=c(paste("low-risk=",length(which(label==0)),sep=""),paste("high-risk=",length(which(label==1)),sep="")),col=c("blue","red"))
text(x=40,y=0.1,paste("log-rank P=",signif(p, digits = 3),"\n","HR=",signif(hr,4),"(","95%CI,",signif(lower.hr,4),"-",signif(upper.hr,4),")"),bty="n",font=2)
##ROC
library(survivalROC)
score=colMeans(com11>0)
fit <- survivalROC(Stime = cli11$survival, status = cli11$survival.event, predict.time = 40, marker=score,method = "KM")
ROC.1<-fit
plot(ROC.1$FP,ROC.1$TP,type='o',pch=16,xlim=c(0,1),ylim=c(0,1),xlab='FP',ylab='TP',main=file1[i],col='red')
abline(0,1)
text(x=0.5,y=0.2,paste("AUC=",signif(fit$AUC,4),seq=""),bty="n",font=2)
##################validation cohorts seperately
cli=c()
labels=c()
com=c()
for(i in 2:3){
	data1=hexps[[i]]
	cli11=hclis[[i]]
	colnames(cli11)[1]="sample"
	gid1=as.numeric(rownames(data1))
	data1=as.matrix(data1)
	com11=data1[match(cpair1[,1],gid1),]-data1[match(cpair1[,2],gid1),]
	com=cbind(com,com11)
	label1=colSums(com11>0)
	label=c()
	label[label1>1]=1
	label[label1<=1]=0
	survi=Surv(cli11$survival,cli11$survival.event)
	cli=rbind(cli,cli11[,c("sample","survival","survival.event")])
	labels=c(labels,label)
	fit<-survfit(survi~label)
	cox=coxph(survi~label)
	b=summary(cox)
	hr=b$conf.int[1]
	lower.hr=b$conf.int[3]
	upper.hr=b$conf.int[4]
	p=b$sctest[3]
	plot(fit,col =c("blue","red"),main=paste(file1[i],"H-meta CTX",sep=" "),ylab = "survival rate",xlab = "survival months",lty=c(1,2),lwd=1.5,mark.time=T)
	legend("bottomright",pch=15:18,lty=c(1,2),merge=F,legend=c(paste("low-risk=",length(which(label==0)),sep=""),paste("high-risk=",length(which(label==1)),sep="")),col=c("blue","red"))
	text(x=40,y=0.1,paste("log-rank P=",signif(p, digits = 3),"\n","HR=",signif(hr,4),"(","95%CI,",signif(lower.hr,4),"-",signif(upper.hr,4),")"),bty="n",font=2)
}
#########validation combined cohort 
survi=Surv(cli$survival,cli$survival.event)
fit<-survfit(survi~labels)
cox=coxph(survi~labels)
b=summary(cox)
hr=b$conf.int[1]
lower.hr=b$conf.int[3]
upper.hr=b$conf.int[4]
p=b$sctest[3]
plot(fit,col =c("blue","red"),main="combined validation",ylab = "survival rate",xlab = "survival months",lty=c(1,2),lwd=1.5,mark.time=T)
legend("bottomright",pch=15:18,lty=c(1,2),merge=F,legend=c(paste("low-risk=",length(which(label==0)),sep=""),paste("high-risk=",length(which(label==1)),sep="")),col=c("blue","red"))
text(x=40,y=0.1,paste("log-rank P=",signif(p, digits = 3),"\n","HR=",signif(hr,4),"(","95%CI,",signif(lower.hr,4),"-",signif(upper.hr,4),")"),bty="n",font=2)
#ROC
score=colMeans(com>0)
fit <- survivalROC(Stime = cli$survival, status = cli$survival.event, predict.time = 40, marker=score,method = "KM")
ROC.1<-fit
plot(ROC.1$FP,ROC.1$TP,type='o',pch=16,xlim=c(0,1),ylim=c(0,1),xlab='FP',ylab='TP',main="combined validation",col='red')
abline(0,1)
text(x=0.5,y=0.2,paste("AUC=",signif(fit$AUC,4),seq=""),bty="n",font=2)
####################data with 5-fu response status
cpair1=read.table("GPS_4.txt")
###---------GSE19860
data1=read.table("gene_exp//GSE19860.txt",sep="\t",header=T,row.names=1,quote="\"",as.is=T)
gid1=as.numeric(rownames(data1))
data1=as.matrix(data1)
cli1=read.table("clinical//GSE19860.txt",sep="\t",header=T,as.is=T)
data1=data1[,match(cli1[,1],colnames(data1))]
com1=data1[match(cpair2[,1],gid1),]-data1[match(cpair2[,2],gid1),]
label1=colSums(com1>0)
label=c()
label[label1>1]=1
label[label1<=1]=0
sum(cli1$state==0&label==0)
sum(cli1$state==1&label==0)
sum(cli1$state==0&label==1)
sum(cli1$state==1&label==1)
######-------GSE28702
data1=read.table("F://crc_location//Test1/data//gene_exp//GSE28702.txt",sep="\t",header=T,row.names=1,quote="\"",as.is=T)
gid1=as.numeric(rownames(data1))
data1=as.matrix(data1)
cli1=read.table("F://crc_location//Test1//data//clinical//GSE28702.txt",sep="\t",header=T,as.is=T)
data11=data1[,cli1$lesion==1]
sample1=colnames(data11)
cli11=cli1[cli1$lesion==1,]
com1=data11[match(cpair1[,1],gid1),]-data11[match(cpair1[,2],gid1),]
label1=colSums(com1>0)
label=c()
label[label1>=1]=1
label[label1<1]=0
#筛选分类一致的标签
data2=read.table("F://crc_location//Test1/data//geneexp2//GSE28702.txt",sep="\t",header=T,row.names=1,quote="\"",as.is=T)
data2=as.matrix(data2)
data21=data2[,cli1$lesion==1]
#原始标签
label2=cli11$state
ssam=cli11[label==0&label2==0,1]
sexp1=data11[,match(ssam,sample1)]
rsam=cli11[label==1&label2==1,1]
rexp1=data11[,match(rsam,sample1)]
tresult=matrix(,length(gid1),4)
for(i in 1:length(gid1)){
	ttest=t.test(rexp1[i,],sexp1[i,])
	tresult[i,1:3]=c(gid1[i],ttest$statistic,ttest$p.value)
}
tresult[,4]=p.adjust(tresult[,3],method="BH")
##按照显著性排序
#tresult1=tresult[order(tresult[,4]),]
tresult=tresult[order(abs(tresult[,2]),decreasing=T),]
tresult1=tresult[tresult[,4]<0.05,]
gene1=tresult1[1:50,1]
#--------------------------
exp1=data21[match(gene1,gid1),]
kc = kmeans(t(exp1), 2)$cluster 
Kmeans=c()
Kmeans[kc==1]="class1"
Kmeans[kc==2]="class2"
GPS=c()
GPS[label==0&label2==0]="A"
GPS[label==1&label2==0]="C"
GPS[label==0&label2==1]="B"
GPS[label==1&label2==1]="D"
library(pheatmap)
ann_col=data.frame(GPS3_couple=GPS,Kmeans=Kmeans)
rownames(ann_col)=cli11$sample
ann_color=list(GPS3_couple=c(A=topo.colors(10)[3],B=heat.colors(10)[5],C=heat.colors(10)[2],D=topo.colors(10)[4]),Kmeans=c(class1=cm.colors(10)[2],class2=cm.colors(10)[9]))
exp1=exp1[,order(Kmeans)]
Kmeans=Kmeans[order(Kmeans)]
pheatmap(exp1,cluster_rows=T,cluster_cols=F,color=colorRampPalette(c("red","white","cornflowerblue"))(100),border_color="grey",legend = T,annotation_col=ann_col,annotation_colors=ann_color)






















