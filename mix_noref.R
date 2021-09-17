#!/usr/bin/env Rscript
library(optimx, quietly=TRUE,warn.conflicts=FALSE)
library(plyr, quietly=TRUE,warn.conflicts=FALSE)
library(dplyr, quietly=TRUE,warn.conflicts=FALSE)
library(argparser, quietly=TRUE,warn.conflicts=FALSE)
suppressPackageStartupMessages(library(tidyverse, quietly=TRUE,warn.conflicts=FALSE))

## ========================= IMPORT PARAMETERS ======================================================

argv<- arg_parser("Parse arguments")
argv <- add_argument(argv, "-i", help="input file")
argv <- add_argument(argv, "-o", help="output file",default="input file_mixAD_v0.tsv")
argv <- add_argument(argv, "-s", help="set a seed (numeric value) for pseudorandom number generator") 
#argv <- add_argument(argv, "-v", help="verbose", flag=TRUE)
#argv <- add_argument(argv, "-e", help="maximum error rate", default="allele_freq_311hu_10nea_3den_sima.txt filtered_D5276_100reads.txt")

args <- parse_args(argv)

infile<-args$i
x_verbose<-args$v
outputfile<-args$o
if ( !is.na(args$s) ){ set.seed(args$s) }
if ( !is.na(args$o) ){ outputfile<-paste0(infile,"_mixAD_v0.tsv") }

## ========================= DEFINE FUNCTIONS =======================================================
generate_genomes=function(size_genome, coverage, prob_genomes,probs_K.types,perr=0.001,summarize_data=TRUE) {
cov_distr<-rpois(size_genome,coverage)
K.types<-sample(names(probs_K.types),size_genome,prob=probs_K.types,replace=TRUE)
K<-K.types2K(K.types)[,c(2:length(prob_genomes),1)]
obs_probs<-probs2obs_probs(prob_genomes,K) %>% apply(MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)) %>% t
obs_probs_4bases<-matrix(0,ncol=4,nrow=nrow(obs_probs))
obs_probs_4bases[,1:ncol(obs_probs)]<-obs_probs
obs_probs_4bases<-(obs_probs_4bases+perr)/(1+perr*4)
res<-apply(cbind(cov_distr,obs_probs_4bases),MARGIN=1,FUN=function(x) rmultinom(1,x[1],prob=x[-1])) %>% t
if (summarize_data){
    res<-table(apply(res,1,FUN=function(x) paste(x,collapse="_")))
    res<- names(res) %>% strsplit("_") %>% unlist %>% as.numeric %>% matrix(nrow=length(res),byrow=TRUE) %>% cbind(res)
    } else
    {
    res<-cbind(1:size_genome,res)
    colnames(res)<-c("pos","A","G","C","T")
    }
res
}

generate_possible_K.types<-function(n_components){
#deterministic version, valid only 2 and 3
if (n_components==1) return(c("11")) else if (n_components==2) return(c("11","12")) else if (n_components==3) return(c("111","112","121","122","123")) else if (n_components==4) return(c("1111","1112","1121","1211","1122","1212","1221","1222","1223","1233","1234"))
}

generate_possible_K<-function(n_components,bases=1:4){
    completions<-generate_possible_K.types(n_components)
    for (i in 1:length(completions)){
    x_completion<-as.numeric(unlist(strsplit(completions[i],"")))
    x_completion<-expand.grid(bases,bases,bases,bases)[,x_completion]
    i_completion<-apply(x_completion,MARGIN=1,FUN=function(x) paste0(x,collapse="."))
    x_completion<-x_completion[!duplicated(i_completion),]
    }
return(x_completion)
}


ngenomes2K<-function(ngenomes,bases=1:4){
## give arrangement arrangement of bases for each genomes (i.e. K) given a number of genomes.
l_temp<-lapply(1:ngenomes,function(x) bases)
l_temp<-expand.grid(l_temp) #col are genomes, rows are different arrangement of bases for each genomes, i.e. K
l_temp
}

K2K.types<-function(K,printAsCharacter=TRUE){
## categorize arrangement of bases for each genomes (i.e. K) in types (previously called configuration, or configuration of divergences)
res<-apply(K,MARGIN=1,FUN=function(x) { seenbases<-c();K.types<-c();
for (iy in 1:ncol(K))
    {
    if ( ! (x[[iy]] %in% seenbases) ){ seenbases<-c(seenbases,x[[iy]]) };
    K.types<-c(K.types,which(x[[iy]] == seenbases)[1])
    };
K.types
})
tres<-t(res)
if (printAsCharacter){ tres<-unname(apply(tres,MARGIN=1,FUN=function(x) paste0(x,collapse="_")))}
tres
}

K.types2K<-function(K.types,ng=ngenomes){
sapply(K.types,function(x) strsplit(x,"_")) %>% unlist %>% as.numeric %>% matrix(ncol=ng,byrow=T)
}

probs2obs_probs<-function(xmixture,K){
## gives probs of observing a given base for each arrangement in K
alphabet<-names(table(unlist(K)))
obs_probs<-matrix(0,ncol=length(alphabet),nrow=nrow(K))
colnames(obs_probs)<-alphabet
xmixture_mat<-matrix(xmixture,nrow=nrow(K),ncol=length(xmixture),byrow=T)
for ( i in 1:length(alphabet)){
    xmixture_mat_t<-xmixture_mat
    xmixture_mat_t[K!=i]<-0
    obs_probs[,i]<-apply(xmixture_mat_t,MARGIN=1,sum)
}
obs_probs
}

log_lik_f<-function(xparams,xdata=mydata,Kt=K,K.typest=K.types,ngenomest=ngenomes)
{
    #parse proportions of genomes
    penalty<-0
    xmixture<-xparams[1:ngenomest-1]
    if (sum(xmixture)>1){penalty<-penalty+(10^9)*(sum(xmixture)-1)^2; xmixture<-xmixture/sum(xmixture)}
    xmixture<-c(xmixture,1-sum(xmixture))
    #parse K.types, i.e. divergence classes
    probs_K.types<-xparams[ngenomest:(length(xparams)-1)]
    if (sum(probs_K.types)>1){penalty<-penalty+(10^9)*(sum(probs_K.types)-1)^2; probs_K.types<-probs_K.types/sum(probs_K.types)}
    probs_K.types<-c(1-sum(probs_K.types),probs_K.types)
    names(probs_K.types)[1]<-paste0(rep(1,ngenomest),collapse="_")
    probK<-rep(0,length(K.typest))
    for (i in 1:length(table(K.typest)))
        {
        xK.type<-names(table(K.typest))[i]
        probK[K.typest==xK.type]<-probs_K.types[names(probs_K.types)==xK.type]/table(K.typest)[i]
        }
    #parse error rates
    perr<-xparams[length(xparams)]
    #compute likelihood likelihood
    obs_probs<-probs2obs_probs(xmixture,Kt)
    obs_probs<-(obs_probs+perr)/(1+perr*4)
    likelihood_all<-sapply(1:nrow(obs_probs), function(zz) apply(xdata[,1:4],MARGIN=1,FUN=function(x) dmultinom(x,prob=obs_probs[zz,]))) %>% t 
    ## VERSION 0 -> OK for two genomes
    log_lik<-likelihood_all*probK
    log_lik<-apply(log_lik,MARGIN=2,FUN=sum) %>% log 
    log_lik<-c(log_lik*xdata[,5]) %>% sum
    log_lik-penalty
}

log_lik_f1<-function(x) 
{ 
xparams<-c(1,0.999,x);
names(xparams)<-c("x1","1_2","perr");
log_lik_f(xparams)
}

mydata<-read.table(infile,header=TRUE)
mydata.table<-mydata[,2:5] %>% apply(MARGIN=1,FUN=function(x) paste0(x,collapse="_")) %>% table
mydata<- sapply(names(mydata.table),function(x) strsplit(x,"_")) %>% unlist %>% as.numeric %>% matrix(ncol=4,byrow=T) %>% cbind(mydata.table) 
#%>% cbind(mydata.table) %>% print

##==========================================================================================================

    ## 1 GENOME
    print("processing model with 1 genome")
    K<-generate_possible_K(2)
    K.types<-K2K.types(K)

    ngenomes<-2
    size_genome<-16500
    coverage<-3
    prob_genomes<-c(1-10^(-16),10^(-16))
    probs_K.types<-c(0.99,0.01)

    names(probs_K.types)<-c("1_1","1_2")
    probK<-rep(0,length(K.types))
    for (i in 1:length(table(K.types))){
    xK.type<-names(table(K.types))[i]
    probK[K.types==xK.type]<-probs_K.types[names(probs_K.types)==xK.type]/table(K.types)[i]
    }

    xparams<-0.01
    res1<-optimx(xparams,log_lik_f1,lower=c(10^-6),upper=c(0.5),method="L-BFGS-B",control=list(maximize=TRUE))
    names(res1)[1]<-"perr"
    res1<-c(res1,p1=1,p2=0,p3=0) %>% unlist
    namesallK<-c(names(table(K2K.types(generate_possible_K(2))))[-1],names(table(K2K.types(generate_possible_K(3))))[-1]) %>% paste0("X",.)
    missing_configurations_in1K<-rep(0,length(namesallK))
    names(missing_configurations_in1K)<-namesallK
    res1<-c(res1,missing_configurations_in1K)
## 2 GENOMES
print("processing model with 2 genomes")
#these are fixed for a given number of potential genomes
K<-generate_possible_K(2)
K.types<-K2K.types(K)

ngenomes<-2
size_genome<-16500
coverage<-3
prob_genomes<-c(0.95,0.05)
probs_K.types<-c(0.95,0.05)

names(probs_K.types)<-c("1_1","1_2")
probK<-rep(0,length(K.types))
for (i in 1:length(table(K.types))){
xK.type<-names(table(K.types))[i]
probK[K.types==xK.type]<-probs_K.types[names(probs_K.types)==xK.type]/table(K.types)[i]
}
xparams<-c(prob_genomes[1],probs_K.types[2],0.01)
names(xparams)<-c("p1","1_2","perr")

xparams<-c(0.85,0.05,0.001)
names(xparams)<-c("p1","1_2","perr")

size_exploration<-100
xparams.df<-data.frame(x1=runif(size_exploration,min=0.5,max=0.99),p1_2=rbeta(size_exploration,0.2,1),perr=10^runif(size_exploration,-6,-1))
names(xparams.df)<-c("p1","1_2","perr")
res_exploration<-apply(xparams.df,MARGIN=1,FUN=function(x) log_lik_f(unlist(x)))
xparams<-unlist(xparams.df[which.max(res_exploration),])
suppressWarnings(res2<-optimx(xparams,log_lik_f,lower=c(0,0,10^(-6)),upper=c(1,1,10^(-1)),method="L-BFGS-B",control=list(maximize=TRUE)))
res2<-c(res2,p2=1-res2$p1,p3=0) %>% unlist
namesallK<-c(names(table(K2K.types(generate_possible_K(3))))[-1]) %>% paste0("X",.)
missing_configurations_in1K<-rep(0,length(namesallK))
names(missing_configurations_in1K)<-namesallK
res2<-c(res2,missing_configurations_in1K)

## 3 GENOMES
print("processing model with 3 genomes")
ngenomes<-3
size_genome<-16500
coverage<-10

prob_genomes<-c(0.7,0.2,0.1)
K<-ngenomes2K(ngenomes)
K.types<-K2K.types(K)
probs_K.types<-table(K.types)
probs_K.types[1:5]<-c(0.8,0.05,0.05,0.05,0.05)
probK<-rep(0,length(K.types))
for (i in 1:length(table(K.types))){
xK.type<-names(table(K.types))[i]
probK[K.types==xK.type]<-probs_K.types[names(probs_K.types)==xK.type]/table(K.types)[i]
}
xparams<-c(prob_genomes[1:(ngenomes-1)],probs_K.types[-1],0.0001)
names(xparams)[length(xparams)]<-"perr"

xparams<-c(prob_genomes[1:(ngenomes-1)],probs_K.types[-1],0.0001)
names(xparams)[length(xparams)]<-"perr"
names(xparams)[1:2]<-c("p1","p2")

lowerlimits<-rep(0,length(xparams));lowerlimits[length(lowerlimits)]<-10^(-6)
upperlimits<-rep(1,length(xparams));upperlimits[length(upperlimits)]<-10^(-1)

#res<-optimx(xparams,log_lik_f,lower=lowerlimits,upper=upperlimits,itnmax=250,control=list(maximize=TRUE))

size_exploration<-100
xparams1<-runif(size_exploration,0.0001,0.9999)
xparams2<-runif(size_exploration,0.0001,1-xparams1)
xparams3<-1-xparams1-xparams2
xparamsx<-apply(cbind(xparams1,xparams2,xparams3),MARGIN=1,FUN=function(x) sort(x,decreasing=TRUE)) %>% t %>% .[,1:2]
xparamsk111<-runif(size_exploration,0.6,0.999)
xparamsk112<-sapply(1:4, function(x) runif(size_exploration,0.00001,1-xparamsk111)) %>% apply(MARGIN=1,FUN=function(x) x/sum(x)) %>% t
xparamsk112<-sapply(1:size_exploration,function(x) xparamsk112[x,]*(1-xparamsk111[x])) %>% t
perr<-10^runif(size_exploration,-6,-1)
xparams.df<-cbind(xparamsx,xparamsk112,perr) %>% as.data.frame
names(xparams.df)<-names(xparams)

res_exploration<-apply(xparams.df,MARGIN=1,FUN=function(x) log_lik_f(unlist(x)))
xparams<-unlist(xparams.df[which.max(res_exploration),])
lowerlimits<-rep(0,length(xparams));lowerlimits[length(lowerlimits)]<-10^(-6)
upperlimits<-rep(1,length(xparams));upperlimits[length(upperlimits)]<-10^(-1)

suppressWarnings(res3<-optimx(xparams,log_lik_f,lower=lowerlimits,upper=upperlimits,method="L-BFGS-B",itnmax=250,control=list(maximize=TRUE)))
res3<-c(res3,p3=1-res3$p1-res3$p2) %>% unlist

namesallK<-c(names(table(K2K.types(generate_possible_K(2))))[-1]) %>% paste0("X",.)
missing_configurations_in1K<-rep(0,length(namesallK))
names(missing_configurations_in1K)<-namesallK
res3<-c(res3,missing_configurations_in1K)

##=======COMBINING RESULTS=======
print("processing final output")
namesallK<-c(names(table(K2K.types(generate_possible_K(2))))[-1],names(table(K2K.types(generate_possible_K(3))))[-1]) %>% paste0("X",.)

allresults<-rbind(res1[order(names(res1))],res2[order(names(res2))],res3[order(names(res3))]) %>% as_tibble %>% select(p1,p2,p3,perr,all_of(namesallK),value,xtime)
allresults$ngenomes<-1:3
allresults$degree.freedom<-c(1,3,7)
allresults$value
allresults$AIC<-2*allresults$degree.freedom-2*allresults$value
suppressWarnings(allresults$lik_ratio_test_2g<-1-pchisq(as.numeric(exp(allresults$value[2]-allresults$value)),df=allresults$degree.freedom[2]-allresults$degree.freedom));allresults$lik_ratio_test_2g[2:3]<-c(1,NA)
suppressWarnings(allresults$lik_ratio_test_3g<-1-pchisq(as.numeric(exp(allresults$value[3]-allresults$value)),df=allresults$degree.freedom[3]-allresults$degree.freedom));allresults$lik_ratio_test_3g[3]<-1
allresults<-allresults[order(allresults$AIC),]
allresults$relL_vsbest<-exp((allresults$AIC[1]-allresults$AIC)/2)
allresults$relL_vs2ndbest<-exp((allresults$AIC[2]-allresults$AIC)/2)
allresults <-allresults %>% relocate(ngenomes,.before=p1)
write.table(allresults,file=outputfile,sep="\t",row.names=FALSE,col.names=TRUE,quote=FALSE)
