# Build private analysis cells directly from authorized SEER exports.
options(stringsAsFactors=FALSE, scipen=999)
args<-commandArgs(TRUE)
stopifnot(length(args)==3L)
out<-args[1];dir.create(out,recursive=TRUE,showWarnings=FALSE)
wr<-function(x,n)write.csv(x,file.path(out,n),row.names=FALSE)
raw<-read.csv(args[2],check.names=FALSE)
p<-read.csv(args[3],check.names=FALSE)
racevar<-'Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)'
required<-c(racevar,'Year of diagnosis','Age recode with single ages and 90+','Histologic Type ICD-O-3')
stopifnot(all(required%in%names(raw)))
z<-data.frame(race=raw[[racevar]],year=as.integer(raw[['Year of diagnosis']]),age=suppressWarnings(as.numeric(sub('^([0-9]+).*','\\1',raw[['Age recode with single ages and 90+']]))),hist=as.integer(raw[['Histologic Type ICD-O-3']]))
stopifnot(!anyNA(z),all(z$age>=20),all(z$year%in%2000:2023))
z$age_start<-pmin(90,5*floor(z$age/5))
required<-c(racevar,'Year_Annual_2000_2023','Age recode with <1 year olds and 90+','Count','Population')
stopifnot(all(required%in%names(p)))
p<-data.frame(race=p[[racevar]],year=as.integer(p[['Year_Annual_2000_2023']]),age_start=suppressWarnings(as.numeric(sub('^([0-9]+).*','\\1',p[['Age recode with <1 year olds and 90+']]))),count=p$Count,population=p$Population)
drop<-is.na(p$age_start)|p$age_start<20
stopifnot(all(p$count[drop]==0),all(p$population[drop]==0))
p<-p[!drop,]
stopifnot(nrow(p)==2160,!anyNA(p),!anyDuplicated(p[c('race','year','age_start')]),all(p$count>=0),all(p$population>=0))
a<-aggregate(rep(1,nrow(z)),z[c('race','year','age_start')],sum);names(a)[4]<-'case_count'
qa<-merge(p,a,by=c('race','year','age_start'),all=TRUE);qa$case_count[is.na(qa$case_count)]<-0
stopifnot(nrow(qa)==2160,all(qa$count==qa$case_count))
known<-p$race!='Non-Hispanic Unknown Race'
stopifnot(all(p$population[known]>0),all(p$population[!known]==0))
std<-c(18257225,17722067,19511370,22179956,22479229,19805793,17224359,13307234,10654272,9409940,8725574,7414559,4900234,2678567,1580606);names(std)<-seq(20,90,5)
# SEER Rare Cancer Classification 15.1 / recode 60, restricted to behavior /3.
codes<-sort(unique(c(8020:8022,8050,8140:8141,8143:8144,8147,8160:8162,8190,8200:8201,8210:8211,8221,8230:8231,8255,8260:8263,8290,8310,8315,8320,8323,8333,8380:8384,8430,8440:8441,8450,8470,8480:8482,8490,8500,8503:8504,8510,8512,8514,8525,8542,8550:8551,8560,8562,8570:8576)))
wr(data.frame(code=codes),'robustness_histology_code_list.csv')
defs<-c('D1_all_malignant','D2_carcinoma_8010_8576','D3_range_8140_8389','D4_strict_8140_3','SEER_variants')
select<-function(d)switch(d,D1_all_malignant=rep(TRUE,nrow(z)),D2_carcinoma_8010_8576=z$hist>=8010&z$hist<=8576,D3_range_8140_8389=z$hist>=8140&z$hist<=8389,D4_strict_8140_3=z$hist==8140,SEER_variants=z$hist%in%codes)
cells<-list()
for(d in defs){
 v<-z[select(d)&z$race!='Non-Hispanic Unknown Race',]
 a<-aggregate(rep(1,nrow(v)),v[c('race','year','age_start')],sum);names(a)[4]<-'count'
 g<-merge(p[known,c('race','year','age_start','population')],a,all.x=TRUE);g$count[is.na(g$count)]<-0
 g$definition<-d;g$age_mid<-g$age_start+2;cells[[d]]<-g
}
allcells<-do.call(rbind,cells);rownames(allcells)<-NULL
wr(subset(allcells,definition!='SEER_variants'),'ethnicity_year_age_definition_cells.csv')
q<-subset(p,known&year==2000);q$ethnic_weight<-q$population/ave(q$population,q$age_start,FUN=sum);q$reference<-'2000'
wr(q[c('race','age_start','ethnic_weight','reference')],'ethnicity_standard_weights.csv')
fitw<-function(x){
 x$group<-factor(x$group,levels=c('older','younger'));x$yc<-x$year-2011.5
 f<-lm(log(rate)~yc*group,weights=(rate/se)^2,data=x);b<-coef(f);V<-vcov(f);crit<-qt(.975,df.residual(f))
 do.call(rbind,lapply(c('older','younger','rate_ratio'),function(term){
 L<-setNames(rep(0,length(b)),names(b));if(term!='rate_ratio')L['yc']<-1;if(term!='older')L['yc:groupyounger']<-1
 beta<-sum(L*b);s<-sqrt(drop(t(L)%*%V%*%L));data.frame(estimand=term,estimate=100*expm1(beta),lower=100*expm1(beta-crit*s),upper=100*expm1(beta+crit*s),p=2*pt(-abs(beta/s),df.residual(f)))
 }))
}
trends<-list();variants<-list();detail<-list();strat<-list()
comparisons<-list(pooled_latest=list(a=2000:2004,b=2019:2023),pooled_prepandemic=list(a=2000:2004,b=2015:2019),single_year=list(a=2000,b=2023))
for(d in defs){
 g<-cells[[d]];g$group<-ifelse(g$age_start<50,'younger','older');g<-merge(g,q[c('race','age_start','ethnic_weight')])
 for(m in c('Age only','Joint')){
 rr<-do.call(rbind,lapply(split(g,interaction(g$year,g$group,drop=TRUE)),function(x){
 if(m=='Age only')x<-aggregate(cbind(count,population)~age_start+year+group,x,sum)
 w<-std[as.character(x$age_start)]/sum(std[as.character(unique(x$age_start))]);if(m=='Joint')w<-w*x$ethnic_weight
 data.frame(year=x$year[1],group=x$group[1],rate=1e5*sum(w*x$count/x$population),se=1e5*sqrt(sum(w^2*x$count/x$population^2)))
 }))
 f<-fitw(rr)
 if(d=='SEER_variants')variants[[m]]<-cbind(method=m,f,total_cases=sum(select(d)),younger_cases=sum(select(d)&z$age<50),unknown_race=sum(select(d)&z$race=='Non-Hispanic Unknown Race')) else trends[[paste(d,m)]]<-cbind(definition=d,reference=if(m=='Age only')'age_only_known' else '2000',window='full',f)
 }
 for(gr in c('younger','older'))for(cmp in names(comparisons)){
 x<-g[g$group==gr,];pool<-function(yrs){a<-aggregate(cbind(count,population)~race+age_start,subset(x,year%in%yrs),sum);a[order(a$age_start,a$race),]}
 a<-pool(comparisons[[cmp]]$a);b<-pool(comparisons[[cmp]]$b);stopifnot(identical(a[c('race','age_start')],b[c('race','age_start')]))
 w<-std[as.character(a$age_start)]/sum(std[as.character(unique(a$age_start))]);p0<-a$population/ave(a$population,a$age_start,FUN=sum);p1<-b$population/ave(b$population,b$age_start,FUN=sum);r0<-1e5*a$count/a$population;r1<-1e5*b$count/b$population
 detail[[paste(d,gr,cmp)]]<-data.frame(definition=if(d=='SEER_variants')'Glandular_subset' else d,group=gr,comparison=cmp,race=a$race,age_start=a$age_start,count0=a$count,count1=b$count,population0=a$population,population1=b$population,age_weight=w,share0=p0,share1=p1,rate0=r0,rate1=r1,composition=w*(p1-p0)*(r0+r1)/2,within_stratum=w*(p1+p0)*(r1-r0)/2)
 }
 if(d=='D1_all_malignant')for(r in unique(g$race)){
 a<-subset(g,race==r);a$yc<-a$year-2011.5;a$young<-as.integer(a$age_start<50)
 f<-glm(count~factor(age_start)+yc+yc:young+offset(log(population)),quasipoisson(),a)
 for(gr in c('older','younger','rate_ratio')){
 L<-setNames(rep(0,length(coef(f))),names(coef(f)));if(gr!='rate_ratio')L['yc']<-1;if(gr!='older')L['yc:young']<-1
 beta<-sum(L*coef(f));s<-sqrt(drop(t(L)%*%vcov(f)%*%L))
 strat[[paste(r,gr)]]<-data.frame(definition=d,race=r,estimand=gr,cases=sum(a$count[a$young==ifelse(gr=='older',0,1)]),estimate=100*expm1(beta),lower=100*expm1(beta-1.96*s),upper=100*expm1(beta+1.96*s),p=2*pnorm(-abs(beta/s)))
 }
 }
}
wr(do.call(rbind,trends),'ethnicity_joint_standardization_trends.csv')
wr(do.call(rbind,variants),'robustness_glandular_trends.csv')
wr(do.call(rbind,detail),'robustness_decomposition_cells.csv')
wr(do.call(rbind,strat),'ethnicity_stratified_group_trends.csv')
audit<-aggregate(rep(1,nrow(z)),list(code=z$hist),sum);names(audit)[2]<-'cases';audit$D2<-audit$code>=8010&audit$code<=8576;audit$D3<-audit$code>=8140&audit$code<=8389;audit$D4<-audit$code==8140;audit$SEER_variants<-audit$code%in%codes
wr(audit,'morphology_audit.csv')
wr(data.frame(check=c('P2 C2 matching cells','Known ethnicity cells','Unknown-race cases'),value=c(nrow(qa),sum(known),sum(z$race=='Non-Hispanic Unknown Race'))),'input_audit.csv')
cat('Authorized exports validated; private inputs prepared.\n')
