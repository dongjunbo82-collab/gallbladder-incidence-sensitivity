# Negative binomial generator sensitivity and annual rate intervals.
args<-commandArgs(TRUE);analysis_root<-normalizePath(args[1],winslash='/')
src<-readLines(args[4],warn=FALSE) # Path to paired_trends.R
eval(parse(text=src[seq_len(grep('^B<-',src)[1]-1)]),envir=.GlobalEnv)
suppressPackageStartupMessages(library(MASS))
B<-10000L;set.seed(261026)
windows<-list('2000-2023'=2000:2023,'2000-2019'=2000:2019,'2000-2021'=2000:2021,'Exclude 2020-2021'=c(2000:2019,2022:2023))
slope<-function(r,v,years){w<-r*r/pmax(v,1e-30);y<-log(pmax(r,1e-30));xx<-years-2011.5;sw<-colSums(w);sx<-colSums(w*xx);sy<-colSums(w*y);(colSums(w*xx*y)-sx*sy/sw)/(colSums(w*xx*xx)-sx*sx/sw)}
nbout<-list();nbdiag<-list();intervals<-list()
full<-subset(allcells,definition=='D1_all_malignant');full$race<-factor(full$race);full$yc<-full$year-2011.5
for(win in names(windows)){
 xfit<-subset(full,year%in%windows[[win]])
 warns<-character()
 ff<-withCallingHandlers(glm.nb(count~ns(age_mid,4)*race+ns(yc,3)*race+ns(age_mid,4):yc+offset(log(population)),data=xfit,control=glm.control(maxit=200)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
 stopifnot(ff$converged,is.finite(ff$theta),ff$theta>0,all(is.finite(coef(ff))),all(fitted(ff)>0))
 xfit$mu<-fitted(ff)
 nbdiag[[win]]<-data.frame(window=win,theta=ff$theta,se_theta=ff$SE.theta,converged=ff$converged,pearson_dispersion=sum(residuals(ff,type='pearson')^2)/df.residual(ff),log_likelihood=as.numeric(logLik(ff)),warnings=paste(unique(warns),collapse=';'),training_years=paste(sort(unique(xfit$year)),collapse=';'))
 for(gr in c('20-49','50+')){
  x<-subset(xfit,group==gr);x<-x[order(x$year,x$age_start,x$race),];years<-sort(unique(x$year))
  aw<-std[as.character(x$age_start)]/sum(std[as.character(unique(x$age_start))]);share<-x$population/ave(x$population,interaction(x$year,x$age_start),FUN=sum)
  ca<-1e5*aw*share/x$population;cj<-1e5*aw*x$ethnic_weight/x$population
  get<-function(N){list(ra=rowsum(N*ca,x$year),rj=rowsum(N*cj,x$year),va=rowsum(N*ca^2,x$year),vj=rowsum(N*cj^2,x$year))}
  calc<-function(r)100*(expm1(slope(r$rj,r$vj,years))-expm1(slope(r$ra,r$va,years)))
  obs<-get(matrix(x$count,ncol=1));gen<-get(matrix(x$mu,ncol=1));d0<-calc(obs);dfit<-calc(gen)
  N<-matrix(rnbinom(nrow(x)*B,mu=x$mu,size=ff$theta),nrow=nrow(x));ds<-calc(get(N));qq<-quantile(ds-dfit,c(.025,.975))
  nbout[[paste(win,gr)]]<-data.frame(definition='D1_all_malignant',group=gr,window=win,method='Paired negative binomial',difference_pp=d0,delta_fit=dfit,median_difference=median(ds),lower=d0-qq[2],upper=d0-qq[1],replicates=B,theta=ff$theta)
 }
}
wr(do.call(rbind,nbout),'negative_binomial_sensitivity.csv');wr(do.call(rbind,nbdiag),'negative_binomial_diagnostics.csv')
# Fay-Feuer gamma intervals for each annual rate; identical estimator to the main analysis.
for(de in defs)for(gr in c('20-49','50+')){
 x<-subset(allcells,definition==de&group==gr)
 aw<-std[as.character(x$age_start)]/sum(std[as.character(unique(x$age_start))]);share<-x$population/ave(x$population,interaction(x$year,x$age_start),FUN=sum)
 for(m in c('Age only','Joint')){
  a<-1e5*aw*if(m=='Age only')share/x$population else x$ethnic_weight/x$population
  for(yy in sort(unique(x$year))){ii<-x$year==yy;y<-sum(a[ii]*x$count[ii]);v<-sum(a[ii]^2*x$count[ii]);u<-max(a[ii])
   intervals[[paste(de,gr,m,yy)]]<-data.frame(definition=de,group=gr,method=m,year=yy,rate=y,variance=v,max_coefficient=u,lower=if(y==0)0 else qgamma(.025,y*y/v,scale=v/y),upper=qgamma(.975,(y+u)^2/(v+u*u),scale=(v+u*u)/(y+u)))
  }
 }
}
wr(do.call(rbind,intervals),'annual_rate_gamma_intervals.csv')
writeLines(capture.output(sessionInfo()),file.path(out,'session_info_uncertainty.txt'))
print(do.call(rbind,nbdiag));print(subset(do.call(rbind,nbout),group=='20-49'))
