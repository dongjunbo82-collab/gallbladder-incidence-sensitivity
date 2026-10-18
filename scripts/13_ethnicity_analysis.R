options(stringsAsFactors=FALSE)
suppressPackageStartupMessages(library(splines))
args<-commandArgs(TRUE)
pkg<-normalizePath(args[1],winslash='/');raw<-normalizePath(args[2],winslash='/')
td<-file.path(pkg,'Analysis_Code_and_Source_Data/derived_results')
wr<-function(x,n)write.csv(x,file.path(td,n),row.names=FALSE)
read<-function(n)read.csv(file.path(raw,n),check.names=FALSE,na.strings=c('NA','~',''),stringsAsFactors=FALSE)
racevar<-'Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)'
p<-read('P2_GB_AllMalignant_YearAge5_RaceEthnicity_SEER17_2000-2023.txt')
p$age_start<-suppressWarnings(as.numeric(sub('^([0-9]+).*','\\1',p[['Age recode with <1 year olds and 90+']])))
p$year<-as.integer(p[['Year_Annual_2000_2023']]);p$race<-p[[racevar]]
excluded<-is.na(p$age_start)|p$age_start<20
stopifnot(all(p$Count[excluded]==0),all(p$Population[excluded]==0))
p<-p[!excluded,];stopifnot(nrow(p)==2160,!anyDuplicated(p[c('race','year','age_start')]))
c<-read('C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt')
c$age<-suppressWarnings(as.numeric(sub('^([0-9]+).*','\\1',c[['Age recode with single ages and 90+']])))
c$year<-as.integer(c[['Year of diagnosis']]);c$race<-c[[racevar]]
c$hist<-as.integer(c[['Histologic Type ICD-O-3']]);c$age_start<-pmin(90,floor(c$age/5)*5)
stopifnot(all(c$age>=20),all(c$year%in%2000:2023),nrow(c)==24509)
counts<-aggregate(rep(1,nrow(c)),c[c('race','year','age_start')],sum);names(counts)[4]<-'c2_count'
qa<-merge(p,counts,by=c('race','year','age_start'),all.x=TRUE);qa$c2_count[is.na(qa$c2_count)]<-0
stopifnot(all(qa$Count==qa$c2_count))
base<-read.csv(file.path(td,'incidence_year_age_definition_cells.csv'))
allp<-aggregate(cbind(Count,Population)~year+age_start,p,sum)
q<-merge(allp,subset(base,definition=='D1_all_malignant'),by=c('year','age_start'))
stopifnot(nrow(q)==360,all(q$Count==q$count),all(q$Population==q$population))
known<-p$race!='Non-Hispanic Unknown Race';stopifnot(all(p$Population[known]>0),all(p$Population[!known]==0))
defs<-unique(base$definition)
select_def<-function(d)if(grepl('^D1',d)) rep(TRUE,nrow(c)) else if(grepl('^D2',d)) c$hist>=8010 & c$hist<=8576 else if(grepl('^D3',d)) c$hist>=8140 & c$hist<=8389 else c$hist==8140
cells<-list();unknown<-list()
for(d in defs){
 z<-c[select_def(d),];unknown[[d]]<-data.frame(definition=d,total=nrow(z),unknown=sum(z$race=='Non-Hispanic Unknown Race'))
 z<-z[z$race!='Non-Hispanic Unknown Race',]
 a<-aggregate(rep(1,nrow(z)),z[c('race','year','age_start')],sum);names(a)[4]<-'count'
 b<-merge(p[known,c('race','year','age_start','Population')],a,all.x=TRUE)
 b$count[is.na(b$count)]<-0;b$definition<-d;names(b)[4]<-'population';b$age_mid<-b$age_start+2
 cells[[d]]<-b
}
cells<-do.call(rbind,cells);rownames(cells)<-NULL
stopifnot(nrow(cells)==7200)
wr(cells,'ethnicity_year_age_definition_cells.csv');wr(do.call(rbind,unknown),'ethnicity_unknown_case_audit.csv')
wr(data.frame(check=c('P2 C2 matching cells','P2 P1 matching cells','Known ethnicity cells','Excluded zero placeholders'),value=c(2160,360,1800,sum(excluded)),passed=TRUE),'ethnicity_export_audit.csv')
std<-c(18257225,17722067,19511370,22179956,22479229,19805793,17224359,13307234,10654272,9409940,8725574,7414559,4900234,2678567,1580606);names(std)<-seq(20,90,5)
pop<-p[known,c('race','year','age_start','Population')];pop$group<-ifelse(pop$age_start<50,'20-49','50+')
composition<-aggregate(Population~race+year+group,pop,sum)
tot<-aggregate(Population~year+group,pop,sum);names(tot)[3]<-'total'
composition<-merge(composition,tot);composition$share<-composition$Population/composition$total
wr(composition,'ethnicity_population_composition.csv')
weights<-list()
for(ref in c('2000','2020','pooled')){
 z<-if(ref=='pooled')pop else subset(pop,year==as.integer(ref))
 a<-aggregate(Population~race+age_start,z,sum);tot<-aggregate(Population~age_start,a,sum);names(tot)[2]<-'total'
 a<-merge(a,tot);a$ethnic_weight<-a$Population/a$total;a$reference<-ref;weights[[ref]]<-a
}
wr(do.call(rbind,weights),'ethnicity_standard_weights.csv')
rates<-list()
for(d in defs){
 z<-subset(cells,definition==d);z$group<-ifelse(z$age_start<50,'younger','older')
 # Age-only comparator uses the same known-ethnicity cases and population as joint rates.
 a<-aggregate(cbind(count,population)~year+age_start+group,z,sum)
 for(ref in c('age_only_known','2000','2020','pooled')){
  zz<-if(ref=='age_only_known')a else merge(z,weights[[ref]][c('race','age_start','ethnic_weight')],by=c('race','age_start'))
  if(ref=='age_only_known')zz$ethnic_weight<-1
  rr<-do.call(rbind,lapply(split(zz,interaction(zz$year,zz$group,drop=TRUE)),function(v){
   w<-std[as.character(v$age_start)]/sum(std[as.character(sort(unique(v$age_start)))])*v$ethnic_weight
   stopifnot(abs(sum(w)-1)<1e-10)
   data.frame(definition=d,reference=ref,year=v$year[1],group=v$group[1],rate=1e5*sum(w*v$count/v$population),se=1e5*sqrt(sum(w^2*v$count/v$population^2)))
  }))
  rates[[paste(d,ref)]]<-rr
 }
}
rates<-do.call(rbind,rates);wr(rates,'ethnicity_joint_standardized_rates.csv')
wins<-list(full=2000:2023,through2019=2000:2019,through2021=2000:2021,exclude2020_2021=c(2000:2019,2022:2023))
fitw<-function(z){
 stopifnot(all(z$rate>0),all(z$se>0));z$group<-factor(z$group,levels=c('older','younger'));z$yc<-z$year-2011.5
 f<-lm(log(rate)~yc*group,weights=(rate/se)^2,data=z);b<-coef(f);v<-vcov(f);q<-qt(.975,df.residual(f))
 do.call(rbind,lapply(c('older','younger','rate_ratio'),function(term){
  L<-setNames(rep(0,length(b)),names(b));if(term!='rate_ratio')L['yc']<-1;if(term!='older')L['yc:groupyounger']<-1
  beta<-sum(L*b);s<-sqrt(drop(t(L)%*%v%*%L))
  data.frame(estimand=term,estimate=100*expm1(beta),lower=100*expm1(beta-q*s),upper=100*expm1(beta+q*s),p=2*pt(-abs(beta/s),df.residual(f)))
 }))
}
res<-list()
for(d in defs)for(ref in unique(rates$reference))for(win in names(wins)){
 z<-subset(rates,definition==d & reference==ref & year%in%wins[[win]])
 res[[paste(d,ref,win)]]<-cbind(definition=d,reference=ref,window=win,fitw(z))
}
res<-do.call(rbind,res);wr(res,'ethnicity_joint_standardization_trends.csv')
# Race-specific grouped slopes: count-based model retains zero-count years.
strat<-list();cont<-list();tests<-list()
contrast<-function(f,L){b<-sum(L*coef(f));s<-sqrt(drop(t(L)%*%vcov(f)%*%L));data.frame(estimate=100*expm1(b),lower=100*expm1(b-1.96*s),upper=100*expm1(b+1.96*s),p=2*pnorm(-abs(b/s)))}
for(d in defs){
 z<-subset(cells,definition==d);z$yc<-z$year-2011.5;z$young<-as.integer(z$age_start<50)
 for(r in unique(z$race)){
  a<-subset(z,race==r)
  f<-glm(count~factor(age_start)+yc+yc:young+offset(log(population)),quasipoisson(),a)
  for(g in c('older','younger','rate_ratio')){
   L<-setNames(rep(0,length(coef(f))),names(coef(f)));if(g!='rate_ratio')L['yc']<-1;if(g!='older')L['yc:young']<-1
   strat[[paste(d,r,g)]]<-cbind(definition=d,race=r,estimand=g,cases=sum(a$count[a$young==ifelse(g=='older',0,1)]),contrast(f,L))
  }
 }
 for(k in 4:6)for(win in names(wins)){
  a<-subset(z,age_start<90 & year%in%wins[[win]]);a$race<-factor(a$race)
  full<-glm(count~race*ns(age_mid,df=k)+race*yc+ns(age_mid,df=k):yc+offset(log(population)),quasipoisson(),a)
  red<-glm(count~race*ns(age_mid,df=k)+race*yc+offset(log(population)),quasipoisson(),a)
  pg<-anova(red,full,test='F')$`Pr(>F)`[2]
  # Shared age-by-year shape with race-specific baseline age profiles and linear slopes.
  # Standardize fitted age-specific incidence to fixed 2000 ethnic weights.
  for(age in c(30,40,50,60,70,80)){
   refw<-weights[['2000']];refw<-subset(refw,age_start==age);rlev<-levels(a$race)
   nd<-expand.grid(race=rlev,yc=c(0,1));nd$race<-factor(nd$race,levels=rlev);nd$age_mid<-age;nd$population<-1
   X<-model.matrix(delete.response(terms(full)),nd);mu<-exp(drop(X%*%coef(full)));w<-refw$ethnic_weight[match(as.character(nd$race),refw$race)]
   i0<-nd$yc==0;i1<-nd$yc==1;R0<-sum(w[i0]*mu[i0]);R1<-sum(w[i1]*mu[i1])
   grad<-colSums(X[i1,,drop=FALSE]*(w[i1]*mu[i1]))/R1-colSums(X[i0,,drop=FALSE]*(w[i0]*mu[i0]))/R0
   b<-log(R1/R0);s<-sqrt(drop(t(grad)%*%vcov(full)%*%grad))
   cont[[paste(d,k,win,age)]]<-data.frame(definition=d,df=k,window=win,age=age,estimate=100*expm1(b),lower=100*expm1(b-1.96*s),upper=100*expm1(b+1.96*s),global_p=pg)
  }
  if(k==4 && win=='full'){
   hetero<-glm(count~race*ns(age_mid,df=k)*yc+offset(log(population)),quasipoisson(),a)
   tests[[d]]<-data.frame(definition=d,age_year_p=pg,race_age_year_p=anova(full,hetero,test='F')$`Pr(>F)`[2],dispersion=summary(full)$dispersion)
   # Fully race-specific age-time shapes: check the shared-shape restriction.
   for(age in c(30,40,50,60,70,80)){
    refw<-subset(weights[['2000']],age_start==age);rlev<-levels(a$race)
    nd<-expand.grid(race=rlev,yc=c(0,1));nd$race<-factor(nd$race,levels=rlev);nd$age_mid<-age;nd$population<-1
    X<-model.matrix(delete.response(terms(hetero)),nd);mu<-exp(drop(X%*%coef(hetero)));w<-refw$ethnic_weight[match(as.character(nd$race),refw$race)]
    i0<-nd$yc==0;i1<-!i0;R0<-sum(w[i0]*mu[i0]);R1<-sum(w[i1]*mu[i1]);grad<-colSums(X[i1,,drop=FALSE]*(w[i1]*mu[i1]))/R1-colSums(X[i0,,drop=FALSE]*(w[i0]*mu[i0]))/R0
    b<-log(R1/R0);s<-sqrt(drop(t(grad)%*%vcov(hetero)%*%grad))
    cont[[paste(d,k,'full_race_specific',age)]]<-data.frame(definition=d,df=k,window='full_race_specific',age=age,estimate=100*expm1(b),lower=100*expm1(b-1.96*s),upper=100*expm1(b+1.96*s),global_p=anova(red,hetero,test='F')$`Pr(>F)`[2])
   }
  }
 }
}
wr(do.call(rbind,strat),'ethnicity_stratified_group_trends.csv');wr(do.call(rbind,cont),'ethnicity_adjusted_continuous_sensitivity.csv');wr(do.call(rbind,tests),'ethnicity_interaction_tests.csv')
writeLines(capture.output(sessionInfo()),file.path(td,'ethnicity_session_info.txt'))
print(subset(res,grepl('^D4',definition)&window=='full'))
print(subset(do.call(rbind,strat),grepl('^D4',definition)))
print(subset(do.call(rbind,cont),grepl('^D4',definition)&df==4&window=='full'))
print(do.call(rbind,tests));print(do.call(rbind,unknown))
