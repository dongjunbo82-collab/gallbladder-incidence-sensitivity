options(stringsAsFactors=FALSE,scipen=999)
suppressPackageStartupMessages(library(splines))
args<-commandArgs(TRUE);root<-normalizePath(args[1],winslash='/');out<-file.path(root,'results');dir.create(out,recursive=TRUE,showWarnings=FALSE)
td<-normalizePath(args[3],winslash='/') # Private inputs produced by prepare_inputs.R
rd<-function(n)read.csv(file.path(td,n),check.names=FALSE)
wr<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE)
std<-c(18257225,17722067,19511370,22179956,22479229,19805793,17224359,13307234,10654272,9409940,8725574,7414559,4900234,2678567,1580606);names(std)<-seq(20,90,5)
cells<-rd('ethnicity_year_age_definition_cells.csv')
# Reconstruct the official SEER variant classification from its code list.
rawfile<-args[2] # Authorized local C2 export; supply its path as argument 2
raw<-read.csv(rawfile,check.names=FALSE)
age<-as.numeric(sub('^([0-9]+).*','\\1',raw[['Age recode with single ages and 90+']]))
z<-data.frame(race=raw[['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)']],year=as.integer(raw[['Year of diagnosis']]),age_start=pmin(90,5*floor(age/5)),hist=as.integer(raw[['Histologic Type ICD-O-3']]))
codes<-rd('robustness_histology_code_list.csv')$code
v<-subset(z,hist%in%codes&race!='Non-Hispanic Unknown Race')
a<-aggregate(rep(1,nrow(v)),v[c('race','year','age_start')],sum);names(a)[4]<-'count'
base<-subset(cells,grepl('^D1',definition))[c('race','year','age_start','population')]
g<-merge(base,a,all.x=TRUE);g$count[is.na(g$count)]<-0;g$definition<-'SEER_variants';g$age_mid<-g$age_start+2
allcells<-rbind(cells,g);allcells$group<-ifelse(allcells$age_start<50,'20-49','50+')
weights<-subset(rd('ethnicity_standard_weights.csv'),reference=='2000')[c('race','age_start','ethnic_weight')]
allcells<-merge(allcells,weights,by=c('race','age_start'))
wr(aggregate(count~definition+group,allcells,sum),'case_counts_known.csv')
wr(data.frame(definition=c('D1_all_malignant','SEER_variants','D4_8140'),all_cases=c(nrow(z),sum(z$hist%in%codes),sum(z$hist==8140)),young_cases=c(sum(z$age_start<50),sum(z$hist%in%codes&z$age_start<50),sum(z$hist==8140&z$age_start<50)),unknown_race=c(sum(z$race=='Non-Hispanic Unknown Race'),sum(z$hist%in%codes&z$race=='Non-Hispanic Unknown Race'),sum(z$hist==8140&z$race=='Non-Hispanic Unknown Race'))),'case_counts_all.csv')
defs<-c('D1_all_malignant','SEER_variants','D4_strict_8140_3')
B<-10000L;set.seed(20261025)
rates<-list();trends<-list();paired<-list();diag<-list();smoothdiag<-list();plotfit<-list();fitrates<-list();refsens<-list()
windows<-list('2000-2023'=2000:2023,'2000-2019'=2000:2019,'2000-2021'=2000:2021,'Exclude 2020-2021'=c(2000:2019,2022:2023))
# Independent group intercepts and slopes in the original WLS model imply
# exactly the same point slopes as separate within-group weighted regressions.
slope<-function(r,v,years){w<-r*r/pmax(v,1e-30);y<-log(pmax(r,1e-30));x<-years-2011.5;sw<-colSums(w);sx<-colSums(w*x);sy<-colSums(w*y);(colSums(w*x*y)-sx*sy/sw)/(colSums(w*x*x)-sx*sx/sw)}
for(def in defs){
 full<-subset(allcells,definition==def);full$yc<-full$year-2011.5;full$race<-factor(full$race)
 fits<-list()
 for(win in names(windows)){
 train<-full[full$year%in%windows[[win]],]
 fitgen<-glm(count~ns(age_mid,df=4)*race+ns(yc,df=3)*race+ns(age_mid,df=4):yc+offset(log(population)),family=poisson(),data=train,control=glm.control(maxit=100))
 stopifnot(fitgen$converged,all(fitted(fitgen)>0),all(is.finite(coef(fitgen))))
 train$mu<-fitted(fitgen);fits[[win]]<-train
 smoothdiag[[paste(def,win)]]<-data.frame(definition=def,window=win,converged=fitgen$converged,pearson_dispersion=sum(residuals(fitgen,type='pearson')^2)/df.residual(fitgen),deviance=deviance(fitgen),df=df.residual(fitgen),min_mean=min(train$mu),training_years=paste(sort(unique(train$year)),collapse=';'))
 }
 full<-fits[['2000-2023']]
 for(group in c('20-49','50+')){
  x<-full[full$group==group,];x<-x[order(x$year,x$age_start,x$race),];yrs<-sort(unique(x$year));K<-nrow(x)/length(yrs)
  aw<-std[as.character(x$age_start)]/sum(std[as.character(unique(x$age_start))])
  # Age-only mixture uses observed ethnic shares within each age/year.
  share<-x$population/ave(x$population,interaction(x$year,x$age_start),FUN=sum)
  ca<-1e5*aw*share/x$population;cj<-1e5*aw*x$ethnic_weight/x$population
  sumyear<-function(mat,coef)rowsum(mat*coef,x$year,reorder=TRUE)
  getrates<-function(N){list(ra=sumyear(N,ca),rj=sumyear(N,cj),va=sumyear(N,ca^2),vj=sumyear(N,cj^2))}
  obs<-getrates(matrix(x$count,ncol=1));gen<-getrates(matrix(x$mu,ncol=1))
  fitrates[[paste(def,group)]]<-data.frame(definition=def,group=group,year=yrs,observed_age_only=obs$ra[,1],fitted_age_only=gen$ra[,1],observed_joint=obs$rj[,1],fitted_joint=gen$rj[,1],fitted_variance_age_only=gen$va[,1],fitted_variance_joint=gen$vj[,1])
  for(method in c('Age only','Joint')){
   r<-if(method=='Age only')obs$ra[,1] else obs$rj[,1];v0<-if(method=='Age only')obs$va[,1] else obs$vj[,1]
   rates[[paste(def,group,method)]]<-data.frame(definition=def,group=group,method=method,year=yrs,rate=r,se=sqrt(v0))
   zfit<-data.frame(year=yrs,yc=yrs-2011.5,r=r,se=sqrt(v0));f<-lm(log(r)~yc,weights=(r/se)^2,zfit);f3<-lm(log(r)~ns(yc,df=3),weights=(r/se)^2,zfit)
   u<-residuals(f)*sqrt(weights(f));acf1<-cor(u[-length(u)],u[-1]);pbox<-Box.test(u,lag=2,type='Ljung-Box')$p.value
   # Fixed-weight Bartlett HAC sandwich, lag 2; a sensitivity, not a new primary model.
   X<-model.matrix(f);ww<-weights(f);sc<-X*as.numeric(ww*residuals(f));meat<-crossprod(sc)
   for(lag in 1:2){tmp<-crossprod(sc[(lag+1):nrow(sc),,drop=FALSE],sc[1:(nrow(sc)-lag),,drop=FALSE]);meat<-meat+(1-lag/3)*(tmp+t(tmp))}
   bread<-solve(crossprod(X*sqrt(ww)));hac<-bread%*%meat%*%bread*nrow(X)/(nrow(X)-ncol(X));hs<-sqrt(hac[2,2]);hc<-qt(.975,22)
   diag[[paste(def,group,method)]]<-data.frame(definition=def,group=group,method=method,nonlinear_p=anova(f,f3)$`Pr(>F)`[2],residual_acf1=acf1,ljung_box_lag2_p=pbox,hac_apc=100*expm1(coef(f)[2]),hac_lower=100*expm1(coef(f)[2]-hc*hs),hac_upper=100*expm1(coef(f)[2]+hc*hs))
   nd<-data.frame(yc=seq(-11.5,11.5,length.out=150));pr<-predict(f3,nd,se.fit=TRUE)
   plotfit[[paste(def,group,method)]]<-data.frame(definition=def,group=group,method=method,year=nd$yc+2011.5,rate=exp(pr$fit),lower=exp(pr$fit-qt(.975,df.residual(f3))*pr$se.fit),upper=exp(pr$fit+qt(.975,df.residual(f3))*pr$se.fit))
  }
  for(win in names(windows)){
   x<-fits[[win]];x<-x[x$group==group,];x<-x[order(x$year,x$age_start,x$race),];yrs<-sort(unique(x$year))
   aw<-std[as.character(x$age_start)]/sum(std[as.character(unique(x$age_start))]);share<-x$population/ave(x$population,interaction(x$year,x$age_start),FUN=sum)
   ca<-1e5*aw*share/x$population;cj<-1e5*aw*x$ethnic_weight/x$population
   obs<-getrates(matrix(x$count,ncol=1));gen<-getrates(matrix(x$mu,ncol=1));idx<-rep(TRUE,length(yrs))
   simN<-list('Paired plug-in Poisson'=matrix(rpois(nrow(x)*B,x$count),nrow=nrow(x)),
             'Paired smooth Poisson'=matrix(rpois(nrow(x)*B,x$mu),nrow=nrow(x)),
             'Paired Jeffreys Gamma'=matrix(rgamma(nrow(x)*B,shape=x$count+.5,rate=1),nrow=nrow(x)))
   ba<-slope(obs$ra[idx,,drop=FALSE],obs$va[idx,,drop=FALSE],yrs[idx]);bj<-slope(obs$rj[idx,,drop=FALSE],obs$vj[idx,,drop=FALSE],yrs[idx])
   observed<-100*(expm1(bj)-expm1(ba))
   ga<-slope(gen$ra[idx,,drop=FALSE],gen$va[idx,,drop=FALSE],yrs[idx]);gj<-slope(gen$rj[idx,,drop=FALSE],gen$vj[idx,,drop=FALSE],yrs[idx]);gd<-100*(expm1(gj)-expm1(ga))
   for(m in names(simN)){
    sr<-getrates(simN[[m]]);sa<-100*expm1(slope(sr$ra[idx,,drop=FALSE],sr$va[idx,,drop=FALSE],yrs[idx]));sj<-100*expm1(slope(sr$rj[idx,,drop=FALSE],sr$vj[idx,,drop=FALSE],yrs[idx]));dif<-sj-sa;q<-quantile(dif,c(.025,.975))
    if(m=='Paired smooth Poisson'){lo<-observed-(q[2]-gd);hi<-observed-(q[1]-gd)}else{lo<-q[1];hi<-q[2]}
    paired[[paste(def,group,win,m)]]<-data.frame(definition=def,group=group,window=win,method=m,age_only_apc=100*expm1(ba),joint_apc=100*expm1(bj),difference_pp=observed,delta_fit=gd,median_difference=median(dif),lower=lo,upper=hi,correlation=cor(sa,sj),replicates=B,zero_cells=sum(x$count==0),total_cells=nrow(x))
   }
   if(def=='D1_all_malignant'&&group=='20-49'&&win=='2000-2023'){
    for(ref in c('2000','2020','Pooled 2000-2023')){
     rp<-if(ref=='Pooled 2000-2023')aggregate(population~race+age_start,x,sum) else subset(x,year==as.integer(ref))[c('race','age_start','population')]
     rp$q<-rp$population/ave(rp$population,rp$age_start,FUN=sum)
     qref<-rp$q[match(paste(x$race,x$age_start),paste(rp$race,rp$age_start))];cj<-1e5*aw*qref/x$population
     oo<-getrates(matrix(x$count,ncol=1));gg<-getrates(matrix(x$mu,ncol=1));ba0<-slope(oo$ra,oo$va,yrs);bj0<-slope(oo$rj,oo$vj,yrs);d0<-100*(expm1(bj0)-expm1(ba0));dfit<-100*(expm1(slope(gg$rj,gg$vj,yrs))-expm1(slope(gg$ra,gg$va,yrs)))
     sr<-getrates(simN[['Paired smooth Poisson']]);sa<-100*expm1(slope(sr$ra,sr$va,yrs));sj<-100*expm1(slope(sr$rj,sr$vj,yrs));qq<-quantile(sj-sa,c(.025,.975))
     # Match the manuscript's joint two-age-group WLS residual variance and df.
     af<-fits[['2000-2023']]
     rpall<-if(ref=='Pooled 2000-2023')aggregate(population~race+age_start,af,sum) else subset(af,year==as.integer(ref))[c('race','age_start','population')]
     rpall$q<-rpall$population/ave(rpall$population,rpall$age_start,FUN=sum)
     af$q<-rpall$q[match(paste(af$race,af$age_start),paste(rpall$race,rpall$age_start))]
     af$coef<-1e5*std[as.character(af$age_start)]/ifelse(af$group=='20-49',sum(std[names(std)%in%as.character(seq(20,45,5))]),sum(std[names(std)%in%as.character(seq(50,90,5))]))*af$q/af$population
     af$rate<-af$coef*af$count;af$variance<-af$coef^2*af$count
     zz<-aggregate(cbind(rate,variance)~year+group,af,sum);zz$group<-factor(zz$group,levels=c('20-49','50+'));zz$yc<-zz$year-2011.5
     rf<-lm(log(rate)~yc*group,weights=rate^2/variance,data=zz);bc<-confint(rf)['yc',]
     refsens[[ref]]<-data.frame(reference=ref,age_only_apc=100*expm1(ba0),joint_apc=100*expm1(bj0),joint_lower=100*expm1(bc[1]),joint_upper=100*expm1(bc[2]),difference_pp=d0,delta_fit=dfit,lower=d0-(qq[2]-dfit),upper=d0-(qq[1]-dfit),replicates=B)
    }
   }
  }
 }
}
wr(do.call(rbind,rates),'annual_rates.csv');wr(do.call(rbind,paired),'paired_apc_differences.csv');wr(do.call(rbind,diag),'time_diagnostics.csv');wr(do.call(rbind,smoothdiag),'bootstrap_count_model_diagnostics.csv');wr(do.call(rbind,plotfit),'nonlinear_time_curves.csv');wr(do.call(rbind,fitrates),'generator_fitted_annual_rates.csv');wr(do.call(rbind,refsens),'reference_weight_sensitivity.csv')
# Group-specific WLS intervals using the two-age-group interaction model.
main<-subset(rd('ethnicity_joint_standardization_trends.csv'),window=='full'&reference%in%c('age_only_known','2000')&definition%in%defs)
vv<-rd('robustness_glandular_trends.csv');vv$definition<-'SEER_variants';vv$reference<-ifelse(vv$method=='Age only','age_only_known','2000');vv$window<-'full'
main<-rbind(main,vv[names(main)]);wr(main,'original_wls_intervals.csv')
detail<-rd('robustness_decomposition_cells.csv');detail$definition[detail$definition=='Glandular_subset']<-'SEER_variants'
detail$group<-ifelse(detail$group=='younger','20-49','50+')
detail<-subset(detail,definition%in%defs)
set.seed(20261026) # Independent stream: period decomposition is unchanged by trend-loop order.
comp<-list();pooled<-list()
for(key in unique(paste(detail$definition,detail$group,detail$comparison))){
 x<-detail[paste(detail$definition,detail$group,detail$comparison)==key,];w<-x$age_weight;p0<-x$share0;p1<-x$share1
 for(m in c('Plug-in Poisson','Jeffreys Gamma')){
  draw<-function(n,p)if(m=='Plug-in Poisson')matrix(rpois(length(n)*B,n),nrow=length(n))/p*1e5 else matrix(rgamma(length(n)*B,n+.5,rate=p),nrow=length(n))*1e5
  s0<-draw(x$count0,x$population0);s1<-draw(x$count1,x$population1)
  simC<-w*(p1-p0)*(s0+s1)/2;simW<-w*(p1+p0)*(s1-s0)/2
  for(race in c('Total',unique(x$race))){ii<-if(race=='Total')rep(TRUE,nrow(x)) else x$race==race
   for(term in c('Composition','Within-stratum','Total change')){
    sim<-if(term=='Composition')colSums(simC[ii,,drop=FALSE]) else if(term=='Within-stratum')colSums(simW[ii,,drop=FALSE]) else colSums((simC+simW)[ii,,drop=FALSE])
    est<-if(term=='Composition')sum(x$composition[ii]) else if(term=='Within-stratum')sum(x$within_stratum[ii]) else sum(x$composition[ii]+x$within_stratum[ii]);q<-quantile(sim,c(.025,.975))
    comp[[paste(key,m,race,term)]]<-data.frame(definition=x$definition[1],group=x$group[1],comparison=x$comparison[1],race=race,component=term,method=m,estimate=est,median=median(sim),lower=q[1],upper=q[2])
   }
  }
 }
}
wr(do.call(rbind,comp),'ethnic_decomposition.csv')
# Population shares and ethnicity-specific absolute rates, including zero cells.
x<-subset(allcells,definition=='D1_all_malignant');x$period<-ifelse(x$year<=2004,'2000-2004',ifelse(x$year>=2019,'2019-2023','Other'))
pop<-aggregate(population~race+group+year,x,sum);pop$share<-pop$population/ave(pop$population,interaction(pop$year,pop$group),FUN=sum);wr(pop,'population_shares.csv')
annual<-aggregate(cbind(count,population)~race+group+year+age_start,x,sum)
ar<-do.call(rbind,lapply(split(annual,interaction(annual$race,annual$group,annual$year,drop=TRUE)),function(v){w<-std[as.character(v$age_start)]/sum(std[as.character(v$age_start)]);a<-1e5*w/v$population;y<-sum(a*v$count);v0<-sum(a*a*v$count);u<-max(a);data.frame(race=v$race[1],group=v$group[1],year=v$year[1],count=sum(v$count),rate=y,lower=if(y==0)0 else qgamma(.025,y*y/v0,scale=v0/y),upper=qgamma(.975,(y+u)^2/(v0+u*u),scale=(v0+u*u)/(y+u)))}));wr(ar,'ethnic_absolute_rates.csv')
pp<-aggregate(cbind(count,population)~race+group+period+age_start,subset(x,period!='Other'),sum)
pr<-do.call(rbind,lapply(split(pp,interaction(pp$race,pp$group,pp$period,drop=TRUE)),function(v){w<-std[as.character(v$age_start)]/sum(std[as.character(v$age_start)]);a<-1e5*w/v$population;y<-sum(a*v$count);v0<-sum(a*a*v$count);u<-max(a);data.frame(race=v$race[1],group=v$group[1],period=v$period[1],count=sum(v$count),population=sum(v$population),rate=y,lower=if(y==0)0 else qgamma(.025,y*y/v0,scale=v0/y),upper=qgamma(.975,(y+u)^2/(v0+u*u),scale=(v0+u*u)/(y+u)))}));pr$share<-pr$population/ave(pr$population,interaction(pr$group,pr$period),FUN=sum);wr(pr,'pooled_ethnic_rates_shares.csv')
writeLines(capture.output(sessionInfo()),file.path(out,'session_info.txt'))
cat('PAIRED D1\n');print(subset(do.call(rbind,paired),definition=='D1_all_malignant'&group=='20-49'&window=='2000-2023'))
cat('TIME DIAGNOSTICS\n');print(subset(do.call(rbind,diag),definition=='D1_all_malignant'&group=='20-49'))
cat('BOOTSTRAP MODEL\n');print(do.call(rbind,smoothdiag))
