options(stringsAsFactors=FALSE)
suppressPackageStartupMessages({library(splines);library(survival);library(cmprsk)})
args <- commandArgs(TRUE)
package_dir <- normalizePath(args[1], winslash="/")
raw_dir <- normalizePath(args[2], winslash="/")
td <- file.path(package_dir,"Analysis_Code_and_Source_Data","derived_results")
wr <- function(x,n) write.csv(x,file.path(td,n),row.names=FALSE)
cells <- read.csv(file.path(td,"incidence_year_age_definition_cells.csv"))
defs <- unique(cells$definition)
# Census P25-1130, SEER 20-age-group standard; normalize within each display group.
std <- c(18257225,17722067,19511370,22179956,22479229,19805793,17224359,13307234,10654272,9409940,8725574,7414559,4900234,2678567,1580606)
names(std) <- seq(20,90,5)
windows <- list(`2000-2023`=2000:2023,`2000-2019`=2000:2019,`2000-2021`=2000:2021,`2000-2023_excluding_2020-2021`=c(2000:2019,2022:2023))
rates <- do.call(rbind,lapply(c(40,45,50,55),function(cutoff) {
 z<-cells;z$group<-ifelse(z$age_start<cutoff,"younger","older")
 do.call(rbind,lapply(split(z,interaction(z$definition,z$year,z$group,drop=TRUE)),function(a){
  w<-std[as.character(a$age_start)];w<-w/sum(w)
  data.frame(definition=a$definition[1],year=a$year[1],cutoff=cutoff,group=a$group[1],
   rate=1e5*sum(w*a$count/a$population),se=1e5*sqrt(sum(w^2*a$count/a$population^2)),count=sum(a$count),population=sum(a$population))
 }))
}))
wr(rates,"reconstructed_age_standardized_rates.csv")
fit_wls<-function(z,window){
 stopifnot(all(z$rate>0),all(z$se>0))
 z$group<-factor(z$group,levels=c("older","younger"));z$year_c<-z$year-2011.5
 f<-lm(log(rate)~year_c*group,weights=(rate/se)^2,data=z)
 b<-coef(f);v<-vcov(f);term<-"year_c:groupyounger";q<-qt(.975,df.residual(f))
 do.call(rbind,lapply(c("older","younger"),function(g){
  cvec<-setNames(rep(0,length(b)),names(b));cvec["year_c"]<-1;if(g=="younger")cvec[term]<-1
  beta<-sum(cvec*b);s<-sqrt(drop(t(cvec)%*%v%*%cvec));si<-sqrt(v[term,term])
  data.frame(definition=z$definition[1],cutoff=z$cutoff[1],window=window,group=g,
   apc=100*expm1(beta),lower_95=100*expm1(beta-q*s),upper_95=100*expm1(beta+q*s),p_value=2*pt(-abs(beta/s),df.residual(f)),
   differential_apc=100*expm1(b[term]),differential_lower_95=100*expm1(b[term]-q*si),differential_upper_95=100*expm1(b[term]+q*si),
   p_interaction=2*pt(-abs(b[term]/si),df.residual(f)),dispersion=summary(f)$sigma^2,
   method="Inverse-variance WLS of directly age-standardized rates",estimand="Annual percent change in younger/older rate ratio")
 }))
}
cutoffs<-do.call(rbind,lapply(split(rates,interaction(rates$definition,rates$cutoff,drop=TRUE)),function(z)
 do.call(rbind,lapply(names(windows),function(w)fit_wls(z[z$year%in%windows[[w]],],w)))))
wr(cutoffs,"age_cutoff_calendar_window_sensitivity.csv")
ints<-read.csv(file.path(td,"age_group_trend_interactions_by_definition_registry.csv"))
for(d in defs){a<-subset(cutoffs,definition==d & cutoff==50 & window=="2000-2023")[1,];i<-which(ints$definition==d & ints$registry=="SEER17")
 ints[i,c("differential_apc","lower_95","upper_95","p_interaction")]<-a[c("differential_apc","differential_lower_95","differential_upper_95","p_interaction")]}
wr(ints,"age_group_trend_interactions_by_definition_registry.csv")
trend<-read.csv(file.path(td,"loglinear_apc_by_definition_registry.csv"))
for(d in defs)for(g in c("younger","older")){
 z<-subset(rates,definition==d & cutoff==50 & group==g);f<-lm(log(rate)~year,weights=(rate/se)^2,data=z);b<-coef(f)[2];s<-sqrt(vcov(f)[2,2]);q<-qt(.975,df.residual(f))
 i<-which(trend$registry=="SEER17" & trend$definition==d & trend$group==ifelse(g=="younger","20-49","50+"))
 trend[i,c("apc","lower_95","upper_95","p_value")]<-c(100*expm1(b),100*expm1(b-q*s),100*expm1(b+q*s),2*pt(-abs(b/s),df.residual(f)))
}
wr(trend,"loglinear_apc_by_definition_registry.csv")
off<-read.csv(file.path(td,"official_age_group_rates_seer17_seer8.csv"));a<-subset(rates,cutoff==50);a$group<-ifelse(a$group=="younger","20-49","50+")
qa<-merge(a,subset(off,registry=="SEER17"),by=c("definition","year","group"),suffixes=c("_reconstructed","_export"))
stopifnot(all(qa$count_reconstructed==qa$count_export),max(abs(qa$rate_reconstructed-qa$rate_export))<0.000501)
wr(qa,"age_standardization_reconstruction_audit.csv")
cont<-do.call(rbind,lapply(defs,function(d)do.call(rbind,lapply(4:6,function(k)do.call(rbind,lapply(names(windows),function(w){
 z<-subset(cells,definition==d & age_start<90 & year%in%windows[[w]]);z$year_c<-z$year-2011.5
 full<-glm(count~ns(age_mid,df=k)*year_c+offset(log(population)),quasipoisson(),z)
 red<-glm(count~ns(age_mid,df=k)+year_c+offset(log(population)),quasipoisson(),z)
 pg<-anova(red,full,test="F")$`Pr(>F)`[2]
 do.call(rbind,lapply(c(30,40,50,60,70,80),function(age){
  nd<-data.frame(age_mid=age,year_c=c(0,1),population=1e5)
  X<-model.matrix(delete.response(terms(full)),nd);v<-X[2,]-X[1,];b<-sum(v*coef(full));s<-sqrt(drop(t(v)%*%vcov(full)%*%v))
  data.frame(definition=d,df=k,window=w,age=age,apc=100*expm1(b),lower_95=100*expm1(b-1.96*s),upper_95=100*expm1(b+1.96*s),global_p=pg,dispersion=summary(full)$dispersion)
 }))
}))))))
wr(cont,"continuous_age_df_window_sensitivity.csv")
# Reuse the exact original variable derivations, but not its analysis/output side effects.
src<-readLines(file.path(package_dir,"Analysis_Code_and_Source_Data/scripts/09_competing_risk_analysis.R"))
start<-grep('^d <- read.csv',src);end<-grep('^step_estimate <-',src)-1
raw_file<-file.path(raw_dir,"C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt")
eval(parse(text=src[start:end]))
flow<-list();coh<-list()
for(def in c("D4","D3"))for(endyr in c(2018,2020)){
 z<-d[if(def=="D4") d$histology==8140 else d$histology>=8140 & d$histology<=8389,]
 lab<-paste0(def,"_2004-",endyr);n0<-nrow(z)
 rows<-data.frame(cohort=lab,step="Definition cohort 2000-2023",excluded=0,remaining=n0)
 flags<-list(`Outside diagnosis window`=!(z$year>=2004 & z$year<=endyr),`Death certificate or autopsy only`=z$dco_autopsy,`Unknown cause of death`=z$unknown_cod,`Missing survival months`=!is.finite(z$time_months))
 keep<-rep(TRUE,n0)
 for(nm in names(flags)){ex<-keep & flags[[nm]];stopifnot(!anyNA(ex));keep[ex]<-FALSE;rows<-rbind(rows,data.frame(cohort=lab,step=nm,excluded=sum(ex),remaining=sum(keep)))}
 flow[[lab]]<-rows;coh[[lab]]<-droplevels(z[keep,])
}
wr(do.call(rbind,flow),"outcome_cohort_sequential_flow.csv")
z<-coh[["D4_2004-2018"]];stopifnot(nrow(z)==9770,!anyNA(z$race))
z$time_months[z$time_months==0]<-.5;z$time_h<-pmin(z$time_months,60);z$status_h<-ifelse(z$time_months>60,0,z$status);z$year_c<-z$year-mean(z$year)
forms<-c("age_group","age_group + sex + race + year_c","age_group + sex + race + year_c + stage","age_group + sex + race + year_c + stage + grade")
out<-list();ph<-list()
for(i in seq_along(forms)){
 rhs<-forms[i];f<-coxph(as.formula(paste("Surv(time_h,status_h==1)~",rhs)),z,x=TRUE)
 mm<-model.matrix(as.formula(paste("~",rhs)),z)[,-1,drop=FALSE]
 fg<-crr(z$time_h,z$status_h,cov1=mm);stopifnot(fg$converged)
 b<-coef(f)[1];s<-sqrt(vcov(f)[1,1]);b2<-fg$coef[1];s2<-sqrt(fg$var[1,1])
 out[[i]]<-data.frame(model=i,covariates=rhs,n=nrow(z),hr=exp(b),hr_lower=exp(b-1.96*s),hr_upper=exp(b+1.96*s),shr=exp(b2),shr_lower=exp(b2-1.96*s2),shr_upper=exp(b2+1.96*s2))
 pt<-cox.zph(f)$table;ph[[i]]<-data.frame(model=i,term=rownames(pt),pt,row.names=NULL)
}
wr(do.call(rbind,out),"sequential_adjustment_mortality.csv");wr(do.call(rbind,ph),"cause_specific_PH_diagnostics.csv")
race<-read.csv(file.path(raw_dir,"R3_GBAC8140_EO20-49_5YRates_byRaceEthnicity_SEER17_2000-2023.csv"),check.names=FALSE)
wr(race,"existing_young_race_period_rates.csv")
cat("Cutoffs D4 full period:\n");print(subset(cutoffs,definition=="D4_strict_8140_3" & window=="2000-2023" & group=="younger"))
print(do.call(rbind,flow));print(do.call(rbind,out));print(ph[[4]])
cat("D4 continuous sensitivity:\n");print(subset(cont,definition=="D4_strict_8140_3" & age%in%c(40,50)))
writeLines(capture.output(sessionInfo()),file.path(td,"revision_session_info.txt"))
