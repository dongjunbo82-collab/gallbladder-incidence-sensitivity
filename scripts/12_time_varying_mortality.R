options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(survival);library(cmprsk)})
args<-commandArgs(TRUE);package_dir<-normalizePath(args[1],winslash="/");raw_dir<-normalizePath(args[2],winslash="/")
td<-file.path(package_dir,"Analysis_Code_and_Source_Data/derived_results")
src<-readLines(file.path(package_dir,"Analysis_Code_and_Source_Data/scripts/09_competing_risk_analysis.R"))
raw_file<-file.path(raw_dir,"C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt")
eval(parse(text=src[grep('^d <- read.csv',src):(grep('^step_estimate <-',src)-1)]))
z<-droplevels(subset(d,histology==8140 & year>=2004 & year<=2018 & !dco_autopsy & !unknown_cod & is.finite(time_months)))
stopifnot(nrow(z)==9770)
z$time_months[z$time_months==0]<-.5;z$time_h<-pmin(z$time_months,60);z$status_h<-ifelse(z$time_months>60,0,z$status);z$year_c<-z$year-mean(z$year);z$young<-as.numeric(z$age_group=="20-49")
# Exploratory age-by-log-time sensitivity, same measured adjustment as primary model.
cs<-coxph(Surv(time_h,status_h==1)~young+tt(young)+sex+race+year_c+stage+grade,data=z,tt=function(x,t,...)x*log1p(t/12))
mm<-model.matrix(~young+sex+race+year_c+stage+grade,z)[,-1,drop=FALSE]
fg<-crr(z$time_h,z$status_h,cov1=mm,cov2=matrix(z$young,ncol=1,dimnames=list(NULL,"young")),tf=function(t)matrix(log1p(t/12),ncol=1))
stopifnot(fg$converged)
out<-list()
for(model in c("Cause-specific Cox","Fine-Gray")){
 b<-if(model=="Fine-Gray")fg$coef else coef(cs);v<-if(model=="Fine-Gray")fg$var else vcov(cs)
 main<-1;tv<-if(model=="Fine-Gray")length(b) else which(names(b)=="tt(young)")
 for(t in c(12,36,60)){cvec<-rep(0,length(b));cvec[main]<-1;cvec[tv]<-log1p(t/12);beta<-sum(cvec*b);se<-sqrt(drop(t(cvec)%*%v%*%cvec))
 out[[length(out)+1]]<-data.frame(model=model,month=t,ratio=exp(beta),lower=exp(beta-1.96*se),upper=exp(beta+1.96*se),age_time_p=2*pnorm(-abs(b[tv]/sqrt(v[tv,tv]))))}
}
res<-do.call(rbind,out);write.csv(res,file.path(td,"time_varying_age_mortality.csv"),row.names=FALSE);print(res)
