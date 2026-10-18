options(stringsAsFactors=FALSE,scipen=999)
args<-commandArgs(TRUE);pkg<-normalizePath(args[1],winslash='/');raw<-normalizePath(args[2],winslash='/')
td<-file.path(pkg,'Analysis_Code_and_Source_Data/derived_results')
rd<-function(n)read.csv(file.path(td,n),check.names=FALSE)
wr<-function(x,n)write.csv(x,file.path(td,n),row.names=FALSE)
c<-read.csv(file.path(raw,'C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt'),check.names=FALSE)
c$year<-as.integer(c[['Year of diagnosis']]);c$age<-as.numeric(sub('^([0-9]+).*','\\1',c[['Age recode with single ages and 90+']]))
c$age_start<-pmin(90,5*floor(c$age/5));c$hist<-as.integer(c[['Histologic Type ICD-O-3']]);c$young<-c$age<50
c$race<-c[['Race and origin recode (NHW, NHB, NHAIAN, NHAPI, Hispanic)']]
cells<-rd('ethnicity_year_age_definition_cells.csv');std<-c(18257225,17722067,19511370,22179956,22479229,19805793,17224359,13307234,10654272,9409940,8725574,7414559,4900234,2678567,1580606);names(std)<-seq(20,90,5)
defs<-unique(cells$definition)
in_def<-function(d)if(grepl('^D1',d))rep(TRUE,nrow(c))else if(grepl('^D2',d))c$hist>=8010&c$hist<=8576 else if(grepl('^D3',d))c$hist>=8140&c$hist<=8389 else c$hist==8140
# Exact observed morphology labels, rather than clinical assumptions from ranges.
audit<-aggregate(rep(1,nrow(c)),list(code=c$hist,label=c[['ICD-O-3 Hist/behav, malignant']]),sum);names(audit)[3]<-'cases'
audit$D1<-TRUE;audit$D2<-audit$code>=8010&audit$code<=8576;audit$D3<-audit$code>=8140&audit$code<=8389;audit$D4<-audit$code==8140
audit$neuroendocrine_sensitivity_exclusion<-audit$code%in%c(8150:8159,8240:8249)
audit<-audit[order(-audit$cases),];wr(audit,'reframe_morphology_audit.csv')
summary<-do.call(rbind,lapply(defs,function(d){q<-c[in_def(d),];data.frame(definition=d,cases=nrow(q),younger=sum(q$young),unknown_race=sum(q$race=='Non-Hispanic Unknown Race'))}))
wr(summary,'reframe_definition_counts.csv')
wr(subset(rd('ethnicity_joint_standardization_trends.csv'),window=='full'&reference%in%c('age_only_known','2000')),'reframe_main_trends.csv')
# Compare continuous-age models in the identical known-ethnicity population.
continuous<-list();suppressPackageStartupMessages(library(splines))
for(d in defs){
 a<-aggregate(cbind(count,population)~year+age_start,subset(cells,definition==d&age_start<90),sum);a$age_mid<-a$age_start+2;a$yc<-a$year-2011.5
 f<-glm(count~ns(age_mid,df=4)*yc+offset(log(population)),quasipoisson(),a)
 f0<-glm(count~ns(age_mid,df=4)+yc+offset(log(population)),quasipoisson(),a);pg<-anova(f0,f,test='F')$`Pr(>F)`[2];stopifnot(f$converged)
 for(age in c(30,40,50,60,70,80)){
  nd<-data.frame(age_mid=c(age,age),yc=c(0,1),population=1);X<-model.matrix(delete.response(terms(f)),nd);L<-X[2,]-X[1,];beta<-sum(L*coef(f));se<-sqrt(drop(t(L)%*%vcov(f)%*%L))
  continuous[[paste(d,age)]]<-data.frame(definition=d,age=age,model='Age only',estimate=100*expm1(beta),lower=100*expm1(beta-1.96*se),upper=100*expm1(beta+1.96*se),global_p=pg)
 }
}
adjusted<-subset(rd('ethnicity_adjusted_continuous_sensitivity.csv'),df==4&window=='full');adjusted$model<-'Race adjusted'
continuous<-rbind(do.call(rbind,continuous),adjusted[c('definition','age','model','estimate','lower','upper','global_p')]);wr(continuous,'reframe_continuous_comparison.csv')
# Symmetric two-component decomposition, conditional on fixed population counts.
# Each comparison uses non-overlapping periods; 5000 plug-in Poisson replicates.
set.seed(20261017);B<-5000L;decomp<-list();detail<-list()
comparisons<-list(pooled_latest=list(a=2000:2004,b=2019:2023),pooled_prepandemic=list(a=2000:2004,b=2015:2019),single_year=list(a=2000,b=2023))
for(d in defs)for(g in c('younger','older'))for(cmp in names(comparisons)){
 z<-subset(cells,definition==d);z<-z[if(g=='younger')z$age_start<50 else z$age_start>=50,]
 pool<-function(yrs){a<-aggregate(cbind(count,population)~race+age_start,subset(z,year%in%yrs),sum);a[order(a$age_start,a$race),]}
 a<-pool(comparisons[[cmp]]$a);b<-pool(comparisons[[cmp]]$b)
 stopifnot(identical(a[c('race','age_start')],b[c('race','age_start')]))
 w<-std[as.character(a$age_start)]/sum(std[as.character(unique(a$age_start))])
 p0<-a$population/ave(a$population,a$age_start,FUN=sum);p1<-b$population/ave(b$population,b$age_start,FUN=sum)
 r0<-1e5*a$count/a$population;r1<-1e5*b$count/b$population
 R0<-sum(w*p0*r0);R1<-sum(w*p1*r1);C<-sum(w*(p1-p0)*(r0+r1)/2);W<-sum(w*(p1+p0)*(r1-r0)/2)
 stopifnot(abs((R1-R0)-C-W)<1e-12)
 s0<-matrix(rpois(nrow(a)*B,a$count),nrow=nrow(a))/a$population*1e5
 s1<-matrix(rpois(nrow(b)*B,b$count),nrow=nrow(b))/b$population*1e5
 sim<-rbind(total=colSums(w*(p1*s1-p0*s0)),composition=colSums(w*(p1-p0)*(s0+s1)/2),within_stratum=colSums(w*(p1+p0)*(s1-s0)/2))
 estimates<-c(total=R1-R0,composition=C,within_stratum=W)
 for(term in names(estimates)){
  ci<-quantile(sim[term,],c(.025,.975));decomp[[paste(d,g,cmp,term)]]<-data.frame(definition=d,group=g,comparison=cmp,component=term,baseline_rate=R0,final_rate=R1,estimate=estimates[[term]],lower=ci[1],upper=ci[2],baseline_cases=sum(a$count),final_cases=sum(b$count),replicates=B)
 }
 detail[[paste(d,g,cmp)]]<-data.frame(definition=d,group=g,comparison=cmp,race=a$race,age_start=a$age_start,count0=a$count,count1=b$count,population0=a$population,population1=b$population,age_weight=w,share0=p0,share1=p1,rate0=r0,rate1=r1,composition=w*(p1-p0)*(r0+r1)/2,within_stratum=w*(p1+p0)*(r1-r0)/2)
}
decomp<-do.call(rbind,decomp);wr(decomp,'reframe_population_decomposition.csv');wr(do.call(rbind,detail),'reframe_decomposition_cells.csv')
# Sex strata: annual directly age-adjusted rates already exported for ages20-49.
s<-read.csv(file.path(raw,'R2_GBAC8140_EO20-49_AnnualRates_bySex_SEER17_2000-2023.csv'),check.names=FALSE)
q<-s[s$Sex%in%c('Male','Female')&grepl('^[0-9]{4} Rate$',s[['Year of diagnosis']]),]
sx<-data.frame(sex=q$Sex,year=as.integer(substr(q[['Year of diagnosis']],1,4)),rate=q[['Age-Adjusted Rate/Trend']],se=q[['Standard Error/P-Value']])
stopifnot(nrow(sx)==48,all(sx$rate>0),all(sx$se>0));sx$sex<-factor(sx$sex,levels=c('Female','Male'));sx$yc<-sx$year-2011.5
fit<-lm(log(rate)~yc*sex,weights=(rate/se)^2,data=sx)
sxres<-list();v<-vcov(fit);b<-coef(fit);tcrit<-qt(.975,df.residual(fit))
for(term in c('Female','Male','male_female_rate_ratio')){
 L<-setNames(rep(0,length(b)),names(b));if(term!='male_female_rate_ratio')L['yc']<-1;if(term!='Female')L['yc:sexMale']<-1
 beta<-sum(L*b);se<-sqrt(drop(t(L)%*%v%*%L));sxres[[term]]<-data.frame(estimand=term,estimate=100*expm1(beta),lower=100*expm1(beta-tcrit*se),upper=100*expm1(beta+tcrit*se),p=2*pt(-abs(beta/se),df.residual(fit)))
}
wr(sx,'reframe_sex_annual_rates.csv');wr(do.call(rbind,sxres),'reframe_sex_trends.csv')
# Native SEER*Stat sex-specific summaries are retained to show rounding sensitivity.
wr(s[s$Sex%in%c('Male','Female')&grepl('APC',s[['Year of diagnosis']]),],'reframe_sex_native_apc.csv')
# Stage-specific population incidence, both age groups, harmonized combined stage.
base<-subset(rd('incidence_year_age_definition_cells.csv'),grepl('^D1',definition));base<-base[c('year','age_start','population')]
z<-c[c$hist==8140&c$year>=2004,];stage<-z[['Combined Summary Stage with Expanded Regional Codes (2004+)']]
z$stage<-ifelse(grepl('^Localized',stage),'Localized',ifelse(grepl('^Regional',stage),'Regional',ifelse(grepl('^Distant',stage),'Distant','Unknown')))
agg<-aggregate(rep(1,nrow(z)),z[c('year','age_start','stage')],sum);names(agg)[4]<-'count'
grid<-merge(subset(base,year>=2004),data.frame(stage=c('Localized','Regional','Distant','Unknown')))
grid<-merge(grid,agg,all.x=TRUE);grid$count[is.na(grid$count)]<-0;grid$group<-ifelse(grid$age_start<50,'younger','older')
grid$period<-paste0(2004+5*floor((grid$year-2004)/5),'-',2008+5*floor((grid$year-2004)/5))
pg<-aggregate(cbind(count,population)~age_start+group+stage+period,grid,sum)
stagerates<-do.call(rbind,lapply(split(pg,interaction(pg$group,pg$stage,pg$period,drop=TRUE)),function(a){
 w<-std[as.character(a$age_start)]/sum(std[as.character(a$age_start)]);rate<-1e5*sum(w*a$count/a$population);se<-1e5*sqrt(sum(w^2*a$count/a$population^2))
 data.frame(group=a$group[1],stage=a$stage[1],period=a$period[1],cases=sum(a$count),population=sum(a$population),rate=rate,se=se,lower=max(0,rate-1.96*se),upper=rate+1.96*se)
}))
wr(stagerates,'reframe_stage_specific_rates.csv');wr(grid,'reframe_stage_year_age_cells.csv')
# Validate against the independent existing stage session, avoiding interpretation of proportions as rates.
old<-read.csv(file.path(raw,'R4b_GBAC8140_EO20-49_StageSpecificRates_SEER17_2004-2023.csv'),check.names=FALSE)
check<-merge(subset(stagerates,group=='younger'),old,by.x=c('stage','period'),by.y=c('STAGE_LRD','PERIOD_STAGE_5Y'))
stopifnot(nrow(check)==16,all(check$cases==check$Count),all(abs(check$rate-check[['Age-Adjusted Rate']])<.00051))
wr(check,'reframe_stage_export_validation.csv')
# D3 range after removal of explicitly neuroendocrine morphology families.
newsel<-c$hist>=8140&c$hist<=8389&!c$hist%in%c(8150:8159,8240:8249)
z<-c[newsel&c$race!='Non-Hispanic Unknown Race',];a<-aggregate(rep(1,nrow(z)),z[c('race','year','age_start')],sum);names(a)[4]<-'count'
pop<-subset(cells,grepl('^D1',definition))[c('race','year','age_start','population')]
clean<-merge(pop,a,all.x=TRUE);clean$count[is.na(clean$count)]<-0
weights<-rd('ethnicity_standard_weights.csv');weights<-subset(weights,reference=='2000')
clean<-merge(clean,weights[c('race','age_start','ethnic_weight')]);clean$group<-ifelse(clean$age_start<50,'younger','older')
rr<-do.call(rbind,lapply(split(clean,interaction(clean$year,clean$group,drop=TRUE)),function(a){w<-std[as.character(a$age_start)]/sum(std[as.character(unique(a$age_start))])*a$ethnic_weight;data.frame(year=a$year[1],group=a$group[1],rate=1e5*sum(w*a$count/a$population),se=1e5*sqrt(sum(w^2*a$count/a$population^2)))}))
rr$yc<-rr$year-2011.5;rr$group<-factor(rr$group,levels=c('older','younger'));f<-lm(log(rate)~yc*group,weights=(rate/se)^2,rr)
b<-coef(f);v<-vcov(f);out<-list()
for(term in c('older','younger','rate_ratio')){L<-setNames(rep(0,length(b)),names(b));if(term!='rate_ratio')L['yc']<-1;if(term!='older')L['yc:groupyounger']<-1;beta<-sum(L*b);se<-sqrt(drop(t(L)%*%v%*%L));out[[term]]<-data.frame(estimand=term,estimate=100*expm1(beta),lower=100*expm1(beta-qt(.975,df.residual(f))*se),upper=100*expm1(beta+qt(.975,df.residual(f))*se),p=2*pt(-abs(beta/se),df.residual(f)),total_cases=sum(newsel),younger_cases=sum(newsel&c$young))}
wr(do.call(rbind,out),'reframe_D3_neuroendocrine_exclusion_sensitivity.csv')
writeLines(capture.output(sessionInfo()),file.path(td,'reframe_session_info.txt'))
print(subset(decomp,grepl('^D4',definition)&comparison=='pooled_latest'));print(do.call(rbind,sxres));print(subset(stagerates,group=='younger'));print(subset(audit,D3&neuroendocrine_sensitivity_exclusion));print(do.call(rbind,out))
