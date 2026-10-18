options(stringsAsFactors=FALSE,scipen=999)
suppressPackageStartupMessages({library(ggplot2);library(patchwork);library(ragg)})
args<-commandArgs(TRUE);pkg<-normalizePath(args[1],winslash='/');qa<-args[2]
td<-file.path(pkg,'Analysis_Code_and_Source_Data/derived_results');fd<-file.path(pkg,'Figures');sd<-file.path(pkg,'Analysis_Code_and_Source_Data/figure_source_data')
dir.create(qa,recursive=TRUE,showWarnings=FALSE)
rd<-function(n)read.csv(file.path(td,n));src<-function(x,n)write.csv(x,file.path(sd,n),row.names=FALSE)
# Contract: quantitative grids, same known-ethnicity population in paired models.
# Figure 1: standardization affects each age group's absolute incidence trend.
# Figure 2: continuous age differences remain but age-specific estimates attenuate.
# Figure 3: composition and within-stratum changes can oppose each other.
# 183-mm width, 600-dpi TIFF only; PNG previews are internal QA products.
theme_set(theme_classic(base_size=8,base_family='Arial')+theme(axis.text=element_text(colour='black'),strip.background=element_blank(),strip.text=element_text(face='bold'),legend.position='bottom',legend.title=element_blank(),plot.tag=element_text(face='bold',size=10),plot.margin=margin(5,6,5,5)))
savefig<-function(g,name,height=125){agg_tiff(file.path(fd,paste0(name,'.tiff')),width=183,height=height,units='mm',res=600,compression='lzw',background='white');print(g);dev.off();agg_png(file.path(qa,paste0(name,'.png')),width=1400,height=round(1400*height/183),res=180,background='white');print(g);dev.off()}
cols<-c('Age only'='#666666','Age + race/ethnicity'='#246B93','Race adjusted'='#246B93')
x<-subset(rd('ethnicity_joint_standardized_rates.csv'),grepl('^D4',definition)&reference%in%c('age_only_known','2000'))
x$method<-ifelse(x$reference=='2000','Age + race/ethnicity','Age only');x$age_group<-factor(x$group,levels=c('younger','older'),labels=c('Ages 20-49','Ages 50+'));x$lower<-pmax(0,x$rate-1.96*x$se);x$upper<-x$rate+1.96*x$se
g<-ggplot(x,aes(year,rate,colour=method,fill=method,linetype=method))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.08,colour=NA)+geom_line(linewidth=.6)+geom_point(size=.7)+facet_wrap(~age_group,scales='free_y',ncol=1)+scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+labs(x='Year of diagnosis',y='Incidence per 100,000',title='Adenocarcinoma NOS incidence by age group')+theme(plot.title=element_text(size=10,face='bold'))
savefig(g,'Figure_1_Age_dependent_incidence',150);src(x,'Figure1_annual_standardized_rates.csv')
x<-subset(rd('reframe_continuous_comparison.csv'),grepl('^D4',definition));x$model<-factor(x$model,levels=c('Age only','Race adjusted'));pd<-position_dodge(width=2)
g<-ggplot(x,aes(age,estimate,colour=model,shape=model))+geom_hline(yintercept=0,linetype=2,colour='#777777')+geom_errorbar(aes(ymin=lower,ymax=upper),position=pd,width=.8,linewidth=.45)+geom_point(position=pd,size=2)+scale_colour_manual(values=cols)+scale_x_continuous(breaks=c(30,40,50,60,70,80))+labs(x='Age (years)',y='Model-derived one-year change (%)',title='Continuous-age estimates before and after race adjustment',subtitle='Identical known-ethnicity population; 95% confidence intervals')+theme(plot.title=element_text(size=10,face='bold'))
savefig(g,'Figure_2_Robustness_and_sensitivity',110);src(x,'Figure2_continuous_age_comparison.csv')
x<-subset(rd('reframe_population_decomposition.csv'),comparison=='pooled_latest');x$definition<-factor(substr(x$definition,1,2),levels=c('D4','D3','D2','D1'));x$age_group<-factor(x$group,levels=c('younger','older'),labels=c('Ages 20-49','Ages 50+'));x$component<-factor(x$component,levels=c('total','composition','within_stratum'),labels=c('Total change','Population composition','Within-stratum rates'))
g<-ggplot(x,aes(estimate,definition,colour=component,shape=component))+geom_vline(xintercept=0,linetype=2,colour='#777777')+geom_errorbar(aes(xmin=lower,xmax=upper),orientation='y',position=position_dodge(width=.55),width=.14,linewidth=.4)+geom_point(position=position_dodge(width=.55),size=1.8)+facet_wrap(~age_group,scales='free_x',ncol=2)+scale_colour_manual(values=c('Total change'='#3A3A3A','Population composition'='#CB7E34','Within-stratum rates'='#246B93'))+labs(x='Absolute incidence change per 100,000',y='Definition',title='Descriptive decomposition of incidence differences',subtitle='2000-2004 to 2019-2023; fixed age standard; 95% bootstrap intervals')+theme(plot.title=element_text(size=10,face='bold'),legend.text=element_text(size=7))
savefig(g,'Figure_3_Population_decomposition',110);src(x,'Figure3_population_decomposition.csv')
# Histology composition is descriptive of code-range membership, not purity.
x<-rd('definition_composition_by_age_period.csv');x$category<-x$definition_category
x$category<-gsub('D4 strict 8140/3','D4 adenocarcinoma NOS',x$category,fixed=TRUE)
x$category<-gsub('D3 other adenocarcinoma','Other D3-range codes',x$category,fixed=TRUE)
x$category<-gsub('D2 other carcinoma','Other D2-range codes',x$category,fixed=TRUE)
g<-ggplot(x,aes(period,percent_within_age_period/100,fill=category))+geom_col(width=.75)+facet_wrap(~age_group,ncol=2)+scale_y_continuous(labels=scales::label_percent())+scale_fill_manual(values=c('#D97706','#2A9D8F','#4C78A8','#AAAAAA'))+labs(x='Diagnosis period',y='Proportion of all malignant tumors',fill=NULL)+theme(axis.text.x=element_text(angle=30,hjust=1),legend.text=element_text(size=6.5))
savefig(g,'Figure_S1_Definition_composition',105);src(x,'FigureS1_definition_composition.csv')
x<-subset(rd('crude_cancer_death_cif.csv'),definition=='D4_strict_8140_3' & window=='2004-2018_5year')
a<-ggplot(x,aes(months,cif,colour=age_group,group=age_group))+geom_errorbar(aes(ymin=lower_95,ymax=upper_95),width=1.8)+geom_line()+geom_point()+scale_x_continuous(breaks=c(12,36,60),labels=c('1','3','5'))+scale_y_continuous(labels=scales::label_percent())+labs(x='Years after diagnosis',y='Cumulative incidence (95% CI)',title='Cancer-specific mortality')
src(x,'FigureS2A_mortality.csv')
x<-rd('adjusted_competing_risk_models.csv');x$label<-factor(c('D4, 5-year','D4, 3-year','D3, 5-year','D3, 3-year'),levels=c('D3, 3-year','D3, 5-year','D4, 3-year','D4, 5-year'))
b<-ggplot(x,aes(fine_gray_shr,label))+geom_vline(xintercept=1,linetype=2)+geom_errorbar(aes(xmin=fine_gray_lower_95,xmax=fine_gray_upper_95),orientation='y',width=.15)+geom_point()+labs(x='Adjusted subdistribution hazard ratio (95% CI)',y=NULL,title='Constant-coefficient summaries')
savefig(a+b+plot_annotation(tag_levels='A'),'Figure_S2_Cancer_specific_mortality',105);src(x,'FigureS2B_mortality.csv')
