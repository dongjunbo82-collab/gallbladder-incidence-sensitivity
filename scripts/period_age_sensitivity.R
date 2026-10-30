# Contextual age/period sensitivity using current D1 adult exports, not exact replication.
options(stringsAsFactors=FALSE,scipen=999)
suppressPackageStartupMessages(library(splines))
args<-commandArgs(TRUE);root<-normalizePath(args[1],winslash='/');out<-file.path(root,'results');dir.create(out,recursive=TRUE,showWarnings=FALSE)
readcells<-function(path){p<-read.csv(path,check.names=FALSE);data.frame(race=p[['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)']],sex=if('Sex'%in%names(p))p$Sex else 'Combined',year=as.integer(p[['Year_Annual_2000_2023']]),age_start=suppressWarnings(as.integer(sub('^([0-9]+).*','\\1',p[['Age recode with <1 year olds and 90+']]))),count=p$Count,population=p$Population)}
x<-subset(readcells(args[2]),sex%in%c('Female','Male')&age_start>=20&!is.na(age_start))
p2<-subset(readcells(args[3]),age_start>=20&!is.na(age_start))
total<-aggregate(cbind(count,population)~race+year+age_start,x,sum)
ck<-merge(total,p2,by=c('race','year','age_start'),suffixes=c('_p3','_p2'),all=TRUE)
stopifnot(nrow(ck)==2160,all(ck$count_p3==ck$count_p2),all(ck$population_p3==ck$population_p2))
p2<-subset(p2,race!='Non-Hispanic Unknown Race');x<-subset(x,race!='Non-Hispanic Unknown Race')
q<-subset(p2,year==2000)[c('race','age_start','population')];q$q<-q$population/ave(q$population,q$age_start,FUN=sum)
combined<-p2[c('race','sex','year','age_start','count','population')]
x<-rbind(combined,x);x<-merge(x,q[c('race','age_start','q')],by=c('race','age_start'))
x<-subset(x,year%in%2000:2020);x$group<-ifelse(x$age_start<60,'20-59','60+');x$age_mid<-x$age_start+2;x$yc<-x$year-2010;x$race<-factor(x$race)
std<-c(18257225,17722067,19511370,22179956,22479229,19805793,17224359,13307234,10654272,9409940,8725574,7414559,4900234,2678567,1580606);names(std)<-seq(20,90,5)
slope<-function(r,v,yr){w<-r*r/pmax(v,1e-30);xx<-yr-2010;y<-log(pmax(r,1e-30));sw<-colSums(w);sx<-colSums(w*xx);sy<-colSums(w*y);(colSums(w*xx*y)-sx*sy/sw)/(colSums(w*xx*xx)-sx*sx/sw)}
B<-10000L;set.seed(261031);ans<-list();annual<-list();diag<-list()
for(sx in c('Combined','Female','Male')){
 train<-subset(x,sex==sx);warns<-character()
 fit<-withCallingHandlers(glm(count~ns(age_mid,4)*race+ns(yc,3)*race+ns(age_mid,4):yc+offset(log(population)),poisson(),data=train,control=glm.control(maxit=100)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
 stopifnot(fit$converged,all(fitted(fit)>0),all(is.finite(coef(fit))))
 train$mu<-fitted(fit)
 diag[[sx]]<-data.frame(sex=sx,window='2000-2020',converged=fit$converged,pearson_dispersion=sum(residuals(fit,type='pearson')^2)/df.residual(fit),min_mean=min(train$mu),warnings=paste(unique(warns),collapse=';'),training_years=paste(sort(unique(train$year)),collapse=';'))
 for(gr in c('20-59','60+')){
  d<-subset(train,group==gr);d<-d[order(d$year,d$age_start,d$race),];yr<-sort(unique(d$year))
  aw<-std[as.character(d$age_start)]/sum(std[as.character(unique(d$age_start))]);share<-d$population/ave(d$population,interaction(d$year,d$age_start),FUN=sum)
  ca<-1e5*aw*share/d$population;cj<-1e5*aw*d$q/d$population
  get<-function(N)list(ra=rowsum(N*ca,d$year),rj=rowsum(N*cj,d$year),va=rowsum(N*ca^2,d$year),vj=rowsum(N*cj^2,d$year))
  calc<-function(v)100*(expm1(slope(v$rj,v$vj,yr))-expm1(slope(v$ra,v$va,yr)))
  obs<-get(matrix(d$count,ncol=1));gen<-get(matrix(d$mu,ncol=1));N<-matrix(rpois(nrow(d)*B,d$mu),nrow=nrow(d));sim<-get(N)
  delta<-calc(obs);deltafit<-calc(gen);qq<-quantile(calc(sim)-deltafit,c(.025,.975));key<-paste(sx,gr)
  ans[[key]]<-data.frame(sex=sx,group=gr,window='2000-2020',cases=sum(d$count),age_only_apc=100*expm1(slope(obs$ra,obs$va,yr)),joint_apc=100*expm1(slope(obs$rj,obs$vj,yr)),difference_pp=delta,delta_fit=deltafit,error_q025=qq[1],error_q975=qq[2],lower=delta-qq[2],upper=delta-qq[1],replicates=B)
  for(m in c('Age only','Joint'))annual[[paste(key,m)]]<-data.frame(sex=sx,group=gr,method=m,year=yr,rate=if(m=='Age only')obs$ra[,1] else obs$rj[,1],variance=if(m=='Age only')obs$va[,1] else obs$vj[,1])
 }
}
a<-do.call(rbind,ans);r<-do.call(rbind,annual)
for(sx in c('Combined','Female','Male'))for(m in c('Age only','Joint')){
 z<-subset(r,sex==sx&method==m);z$group<-factor(z$group,levels=c('20-59','60+'));z$yc<-z$year-2010
 f<-lm(log(rate)~yc*group,weights=rate^2/variance,data=z);V<-vcov(f);b<-coef(f);crit<-qt(.975,df.residual(f))
 for(gr in c('20-59','60+')){
  L<-setNames(rep(0,length(b)),names(b));L['yc']<-1;if(gr=='60+')L['yc:group60+']<-1
  beta<-sum(L*b);se<-sqrt(drop(t(L)%*%V%*%L));ii<-a$sex==sx&a$group==gr;prefix<-if(m=='Age only')'age_only' else 'joint'
  a[ii,paste0(prefix,'_lower')]<-100*expm1(beta-crit*se);a[ii,paste0(prefix,'_upper')]<-100*expm1(beta+crit*se)
 }
}
write.csv(a,file.path(out,'period_age_sensitivity.csv'),row.names=FALSE)
write.csv(r,file.path(out,'period_age_annual_rates.csv'),row.names=FALSE)
write.csv(do.call(rbind,diag),file.path(out,'period_age_generator_diagnostics.csv'),row.names=FALSE)
print(a);print(do.call(rbind,diag))
