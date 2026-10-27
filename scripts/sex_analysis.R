options(stringsAsFactors=FALSE,scipen=999)
suppressPackageStartupMessages(library(splines))
args<-commandArgs(TRUE);root<-normalizePath(args[1],winslash='/');out<-file.path(root,'results');dir.create(out,recursive=TRUE,showWarnings=FALSE)
p3<-read.csv(args[2],check.names=FALSE);p2<-read.csv(args[3],check.names=FALSE)
racecol<-'Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)';agecol<-'Age recode with <1 year olds and 90+'
prep<-function(d){data.frame(race=d[[racecol]],sex=if('Sex'%in%names(d))d$Sex else 'Combined',year=as.integer(d[['Year_Annual_2000_2023']]),age_start=suppressWarnings(as.integer(sub('^([0-9]+).*','\\1',d[[agecol]]))),count=d$Count,population=d$Population)}
x<-subset(prep(p3),sex%in%c('Male','Female')&age_start>=20&!is.na(age_start));base<-subset(prep(p2),age_start>=20&!is.na(age_start))
tot<-aggregate(cbind(count,population)~race+year+age_start,x,sum);ck<-merge(tot,base,by=c('race','year','age_start'),suffixes=c('_p3','_p2'),all=TRUE)
stopifnot(nrow(ck)==2160,all(ck$count_p3==ck$count_p2),all(ck$population_p3==ck$population_p2))
unknown<-aggregate(count~sex,subset(x,race=='Non-Hispanic Unknown Race'),sum)
x<-subset(x,race!='Non-Hispanic Unknown Race');stopifnot(nrow(x)==3600,!anyDuplicated(x[c('sex','race','year','age_start')]),all(x$population>0),all(x$count>=0))
x$group<-ifelse(x$age_start<50,'20-49','50+');x$age_mid<-x$age_start+2;x$yc<-x$year-2011.5;x$race<-factor(x$race)
q<-aggregate(population~race+age_start,subset(x,year==2000),sum);q$common_q<-q$population/ave(q$population,q$age_start,FUN=sum)
x<-merge(x,q[c('race','age_start','common_q')],by=c('race','age_start'))
qs<-subset(x,year==2000)[c('sex','race','age_start','population')];qs$own_q<-qs$population/ave(qs$population,interaction(qs$sex,qs$age_start),FUN=sum)
x<-merge(x,qs[c('sex','race','age_start','own_q')],by=c('sex','race','age_start'))
std<-c(18257225,17722067,19511370,22179956,22479229,19805793,17224359,13307234,10654272,9409940,8725574,7414559,4900234,2678567,1580606);names(std)<-seq(20,90,5)
B<-10000L;set.seed(2610261);windows<-list('2000-2023'=2000:2023,'2000-2019'=2000:2019)
wr<-function(z,n)write.csv(z,file.path(out,n),row.names=FALSE)
slope<-function(r,v,years){w<-r*r/pmax(v,1e-30);y<-log(pmax(r,1e-30));xx<-years-2011.5;sw<-colSums(w);sx<-colSums(w*xx);sy<-colSums(w*y);(colSums(w*xx*y)-sx*sy/sw)/(colSums(w*xx*xx)-sx*sx/sw)}
annual<-list();paired<-list();fitdiag<-list();diffs<-list();fitdiffs<-list();observed<-list()
for(win in names(windows))for(sx in c('Female','Male')){
 train<-subset(x,sex==sx&year%in%windows[[win]]);warns<-character()
 ff<-withCallingHandlers(glm(count~ns(age_mid,4)*race+ns(yc,3)*race+ns(age_mid,4):yc+offset(log(population)),family=poisson(),data=train,control=glm.control(maxit=100)),warning=function(w){warns<<-c(warns,conditionMessage(w));invokeRestart('muffleWarning')})
 stopifnot(ff$converged,all(fitted(ff)>0),all(is.finite(coef(ff))))
 train$mu<-fitted(ff)
 fitdiag[[paste(win,sx)]]<-data.frame(window=win,sex=sx,converged=ff$converged,pearson_dispersion=sum(residuals(ff,type='pearson')^2)/df.residual(ff),min_mean=min(train$mu),warnings=paste(unique(warns),collapse=';'),training_years=paste(sort(unique(train$year)),collapse=';'))
 for(gr in c('20-49','50+')){
  d<-subset(train,group==gr);d<-d[order(d$year,d$age_start,d$race),];yrs<-sort(unique(d$year));aw<-std[as.character(d$age_start)]/sum(std[as.character(unique(d$age_start))])
  share<-d$population/ave(d$population,interaction(d$year,d$age_start),FUN=sum);ca<-1e5*aw*share/d$population
  N<-matrix(rpois(nrow(d)*B,d$mu),nrow=nrow(d))
  for(ref in c('Common 2000','Sex-specific 2000')){
   cj<-1e5*aw*if(ref=='Common 2000')d$common_q/d$population else d$own_q/d$population
   get<-function(nn){list(ra=rowsum(nn*ca,d$year),rj=rowsum(nn*cj,d$year),va=rowsum(nn*ca^2,d$year),vj=rowsum(nn*cj^2,d$year))}
   obs<-get(matrix(d$count,ncol=1));gen<-get(matrix(d$mu,ncol=1));sim<-get(N)
   ba<-100*expm1(slope(obs$ra,obs$va,yrs));bj<-100*expm1(slope(obs$rj,obs$vj,yrs));d0<-bj-ba
   gd<-100*(expm1(slope(gen$rj,gen$vj,yrs))-expm1(slope(gen$ra,gen$va,yrs)))
   ds<-100*(expm1(slope(sim$rj,sim$vj,yrs))-expm1(slope(sim$ra,sim$va,yrs)));qq<-quantile(ds-gd,c(.025,.975))
   key<-paste(win,gr,ref,sx,sep='|');diffs[[key]]<-ds;fitdiffs[[key]]<-gd;observed[[key]]<-d0
   paired[[key]]<-data.frame(window=win,sex=sx,group=gr,reference=ref,age_only_apc=ba,joint_apc=bj,difference_pp=d0,delta_fit=gd,error_q025=qq[1],error_q975=qq[2],lower=d0-qq[2],upper=d0-qq[1],zero_cells=sum(d$count==0),total_cells=nrow(d),cases=sum(d$count),replicates=B)
   for(method in c('Age only','Joint')){
    rr<-if(method=='Age only')obs$ra[,1] else obs$rj[,1];vv<-if(method=='Age only')obs$va[,1] else obs$vj[,1]
    annual[[paste(key,method)]]<-data.frame(window=win,sex=sx,group=gr,reference=ref,method=method,year=yrs,rate=rr,variance=vv)
   }
  }
 }
}
annual<-do.call(rbind,annual);paired<-do.call(rbind,paired)
# Match primary paper: within each sex/method/reference/window fit both age groups together.
for(win in names(windows))for(sx in c('Female','Male'))for(ref in c('Common 2000','Sex-specific 2000'))for(m in c('Age only','Joint')){
 z<-subset(annual,window==win&sex==sx&reference==ref&method==m);z$group<-factor(z$group,levels=c('20-49','50+'));z$yc<-z$year-2011.5
 lmfit<-lm(log(rate)~yc*group,weights=rate^2/variance,data=z);V<-vcov(lmfit);beta<-coef(lmfit);tc<-qt(.975,df.residual(lmfit))
 for(gr in c('20-49','50+')){
  L<-setNames(rep(0,length(beta)),names(beta));L['yc']<-1;if(gr=='50+')L['yc:group50+']<-1
  est<-sum(L*beta);se<-sqrt(drop(t(L)%*%V%*%L));ii<-paired$window==win&paired$sex==sx&paired$reference==ref&paired$group==gr
  prefix<-if(m=='Age only')'age_only' else 'joint'
  paired[ii,paste0(prefix,'_lower')]<-100*expm1(est-tc*se);paired[ii,paste0(prefix,'_upper')]<-100*expm1(est+tc*se)
  stopifnot(abs(paired[ii,paste0(prefix,'_apc')]-100*expm1(est))<1e-10)
 }
}
between<-list()
for(win in names(windows))for(gr in c('20-49','50+'))for(ref in c('Common 2000','Sex-specific 2000')){
 fk<-paste(win,gr,ref,'Female',sep='|');mk<-paste(win,gr,ref,'Male',sep='|');d0<-observed[[fk]]-observed[[mk]];gf<-fitdiffs[[fk]]-fitdiffs[[mk]];ds<-diffs[[fk]]-diffs[[mk]];qq<-quantile(ds-gf,c(.025,.975))
 between[[paste(win,gr,ref)]]<-data.frame(window=win,group=gr,reference=ref,contrast='Female delta minus male delta',estimate=d0,delta_fit=gf,error_q025=qq[1],error_q975=qq[2],lower=d0-qq[2],upper=d0-qq[1],replicates=B)
}
wr(paired,'sex_paired_trends.csv');wr(annual,'sex_annual_rates.csv');wr(do.call(rbind,between),'sex_between_differences.csv');wr(do.call(rbind,fitdiag),'sex_generator_diagnostics.csv');wr(x,'sex_analysis_cells_PRIVATE.csv');wr(unknown,'sex_unknown_race_counts.csv')
wr(data.frame(check=c('Eligible known-ethnicity cells','P2 count mismatch','P2 population mismatch','Known-ethnicity cases','Unknown-ethnicity cases'),value=c(nrow(x),sum(ck$count_p3!=ck$count_p2),sum(ck$population_p3!=ck$population_p2),sum(x$count),sum(unknown$count))),'sex_export_audit.csv')
writeLines(capture.output(sessionInfo()),file.path(out,'session_info_sex_analysis.txt'))
print(subset(paired,group=='20-49'));print(subset(do.call(rbind,between),group=='20-49'));print(do.call(rbind,fitdiag))
