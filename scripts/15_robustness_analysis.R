options(stringsAsFactors=FALSE,scipen=999)
args<-commandArgs(TRUE);pkg<-normalizePath(args[1],winslash='/');raw<-normalizePath(args[2],winslash='/');td<-file.path(pkg,'Analysis_Code_and_Source_Data/derived_results')
rd<-function(n)read.csv(file.path(td,n),check.names=FALSE);wr<-function(x,n)write.csv(x,file.path(td,n),row.names=FALSE)
std<-c(18257225,17722067,19511370,22179956,22479229,19805793,17224359,13307234,10654272,9409940,8725574,7414559,4900234,2678567,1580606);names(std)<-seq(20,90,5)
# Exact SEER Rare Cancer Classification 15.1: adenocarcinoma with variants,
# C23.9. Apply malignant /3 from the case cohort (source permits 0-3).
# Includes designated variants such as adenosquamous and neuroendocrine
# differentiation; not a histologically pure adenocarcinoma definition.
codes<-sort(unique(c(8020:8022,8050,8140:8141,8143:8144,8147,8160:8162,8190,8200:8201,8210:8211,8221,8230:8231,8255,8260:8263,8290,8310,8315,8320,8323,8333,8380:8384,8430,8440:8441,8450,8470,8480:8482,8490,8500,8503:8504,8510,8512,8514,8525,8542,8550:8551,8560,8562,8570:8576)))
c<-read.csv(file.path(raw,'C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt'),check.names=FALSE)
c$year<-as.integer(c[['Year of diagnosis']]);c$age<-as.numeric(sub('^([0-9]+).*','\\1',c[['Age recode with single ages and 90+']]))
c$age_start<-pmin(90,5*floor(c$age/5));c$hist<-as.integer(c[['Histologic Type ICD-O-3']]);c$race<-c[['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)']]
sel<-c$hist%in%codes;z<-c[sel&c$race!='Non-Hispanic Unknown Race',]
audit<-rd('reframe_morphology_audit.csv');audit$glandular_subset<-audit$code%in%codes;wr(audit,'robustness_histology_membership.csv');wr(data.frame(code=codes),'robustness_histology_code_list.csv')
cells<-rd('ethnicity_year_age_definition_cells.csv');pop<-subset(cells,grepl('^D1',definition))[c('race','year','age_start','population')]
a<-aggregate(rep(1,nrow(z)),z[c('race','year','age_start')],sum);names(a)[4]<-'count';g<-merge(pop,a,all.x=TRUE);g$count[is.na(g$count)]<-0;g$definition<-'Glandular_subset';g$group<-ifelse(g$age_start<50,'younger','older')
weights<-subset(rd('ethnicity_standard_weights.csv'),reference=='2000')[c('race','age_start','ethnic_weight')]
g<-merge(g,weights);rates<-list()
for(m in c('Age only','Joint'))for(key in unique(paste(g$year,g$group))){x<-g[paste(g$year,g$group)==key,];w<-std[as.character(x$age_start)]/sum(std[as.character(unique(x$age_start))]);if(m=='Joint')w<-w*x$ethnic_weight else {x<-aggregate(cbind(count,population)~age_start+year+group,x,sum);w<-std[as.character(x$age_start)]/sum(std[as.character(x$age_start)])};rates[[paste(m,key)]]<-data.frame(method=m,year=x$year[1],group=x$group[1],rate=1e5*sum(w*x$count/x$population),se=1e5*sqrt(sum(w^2*x$count/x$population^2)))}
rates<-do.call(rbind,rates);wr(rates,'robustness_glandular_rates.csv');trends<-list()
for(m in unique(rates$method)){x<-subset(rates,method==m);x$yc<-x$year-2011.5;x$group<-factor(x$group,levels=c('older','younger'));f<-lm(log(rate)~yc*group,weights=(rate/se)^2,x);b<-coef(f);v<-vcov(f)
 for(term in c('older','younger','rate_ratio')){L<-setNames(rep(0,length(b)),names(b));if(term!='rate_ratio')L['yc']<-1;if(term!='older')L['yc:groupyounger']<-1;beta<-sum(L*b);se<-sqrt(drop(t(L)%*%v%*%L));trends[[paste(m,term)]]<-data.frame(method=m,estimand=term,estimate=100*expm1(beta),lower=100*expm1(beta-qt(.975,df.residual(f))*se),upper=100*expm1(beta+qt(.975,df.residual(f))*se),p=2*pt(-abs(beta/se),df.residual(f)),total_cases=sum(sel),younger_cases=sum(sel&c$age<50),unknown_race=sum(sel&c$race=='Non-Hispanic Unknown Race'))}}
wr(do.call(rbind,trends),'robustness_glandular_trends.csv')
detail<-rd('reframe_decomposition_cells.csv');cmp<-list(pooled_latest=list(a=2000:2004,b=2019:2023),pooled_prepandemic=list(a=2000:2004,b=2015:2019),single_year=list(a=2000,b=2023))
for(group in c('younger','older'))for(comparison in names(cmp)){x<-g[g$group==group,];pool<-function(yrs){a<-aggregate(cbind(count,population)~race+age_start,subset(x,year%in%yrs),sum);a[order(a$age_start,a$race),]};a<-pool(cmp[[comparison]]$a);b<-pool(cmp[[comparison]]$b);w<-std[as.character(a$age_start)]/sum(std[as.character(unique(a$age_start))]);p0<-a$population/ave(a$population,a$age_start,FUN=sum);p1<-b$population/ave(b$population,b$age_start,FUN=sum);r0<-1e5*a$count/a$population;r1<-1e5*b$count/b$population
 detail<-rbind(detail,data.frame(definition='Glandular_subset',group=group,comparison=comparison,race=a$race,age_start=a$age_start,count0=a$count,count1=b$count,population0=a$population,population1=b$population,age_weight=w,share0=p0,share1=p1,rate0=r0,rate1=r1,composition=w*(p1-p0)*(r0+r1)/2,within_stratum=w*(p1+p0)*(r1-r0)/2))}
set.seed(20261018);B<-50000L;out<-list();zero<-list()
for(key in unique(paste(detail$definition,detail$group,detail$comparison))){x<-detail[paste(detail$definition,detail$group,detail$comparison)==key,];w<-x$age_weight;p0<-x$share0;p1<-x$share1;r0<-x$rate0;r1<-x$rate1
 est<-c(total=sum(w*(p1*r1-p0*r0)),composition=sum(x$composition),within_stratum=sum(x$within_stratum))
 for(t in 0:1){n<-x[[paste0('count',t)]];share<-x[[paste0('share',t)]];zero[[paste(key,t)]]<-data.frame(definition=x$definition[1],group=x$group[1],comparison=x$comparison[1],period=t,cells=length(n),zero_cells=sum(n==0),zero_age_ethnic_weight=sum(w[n==0]*share[n==0]),zero_abs_composition_weight=sum(abs(w[n==0]*(p1[n==0]-p0[n==0])/2)),zero_within_weight=sum(w[n==0]*(p1[n==0]+p0[n==0])/2))}
 for(method in c('Plugin Poisson','Jeffreys Gamma posterior')){
  draw<-function(n,p)if(method=='Plugin Poisson')matrix(rpois(length(n)*B,n),nrow=length(n))/p*1e5 else matrix(rgamma(length(n)*B,shape=n+.5,rate=p),nrow=length(n))*1e5
  s0<-draw(x$count0,x$population0);s1<-draw(x$count1,x$population1)
  sim<-rbind(total=colSums(w*(p1*s1-p0*s0)),composition=colSums(w*(p1-p0)*(s0+s1)/2),within_stratum=colSums(w*(p1+p0)*(s1-s0)/2))
  for(term in names(est)){ci<-quantile(sim[term,],c(.025,.975));out[[paste(key,method,term)]]<-data.frame(definition=x$definition[1],group=x$group[1],comparison=x$comparison[1],component=term,method=method,observed_estimate=est[term],simulation_median=median(sim[term,]),lower=ci[1],upper=ci[2],replicates=B)}
 }}
wr(detail,'robustness_decomposition_cells.csv');wr(do.call(rbind,zero),'robustness_zero_cells.csv');wr(do.call(rbind,out),'robustness_decomposition_intervals.csv')
# Fay-Feuer (1997) conservative gamma interval, max weight over ALL age bands.
grid<-rd('reframe_stage_year_age_cells.csv');pg<-aggregate(cbind(count,population)~age_start+group+stage+period,grid,sum)
stage<-do.call(rbind,lapply(split(pg,interaction(pg$group,pg$stage,pg$period,drop=TRUE)),function(x){w<-std[as.character(x$age_start)]/sum(std[as.character(x$age_start)]);a<-1e5*w/x$population;y<-sum(a*x$count);v<-sum(a^2*x$count);u<-max(a);lo<-if(y==0)0 else qgamma(.025,shape=y^2/v,scale=v/y);hi<-qgamma(.975,shape=(y+u)^2/(v+u^2),scale=(v+u^2)/(y+u));data.frame(group=x$group[1],stage=x$stage[1],period=x$period[1],cases=sum(x$count),rate=y,se=sqrt(v),lower=lo,upper=hi,max_weight=u)}))
wr(stage,'robustness_stage_gamma_intervals.csv');writeLines(capture.output(sessionInfo()),file.path(td,'robustness_session_info.txt'))
print(subset(do.call(rbind,out),grepl('^D4|^Gland',definition)&group=='younger'));print(do.call(rbind,trends));print(subset(do.call(rbind,zero),grepl('^D4',definition)&group=='younger'));print(subset(stage,group=='younger'&stage=='Unknown'))
