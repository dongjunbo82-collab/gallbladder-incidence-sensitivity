# Code-only entry point; all inputs and generated results remain outside this repository.
args<-commandArgs(TRUE)
if(length(args)<2L||any(args%in%c('-h','--help'))){
 cat('Usage: Rscript run_analysis.R OUTPUT_DIRECTORY RAW_DIRECTORY [--figures]\n')
 quit(status=if(any(args%in%c('-h','--help')))0L else 1L)
}
f<-grep('^--file=',commandArgs(FALSE),value=TRUE)
project<-dirname(normalizePath(sub('^--file=','',f),winslash='/'))
raw<-normalizePath(args[2],winslash='/',mustWork=TRUE)
out<-file.path(normalizePath(dirname(args[1]),winslash='/',mustWork=TRUE),basename(args[1]))
inside<-function(a,b)identical(tolower(a),tolower(b))||startsWith(tolower(a),paste0(tolower(b),'/'))
if(inside(out,raw)||inside(out,project)||inside(raw,out)||inside(project,out))stop('Output must be separate from input and code directories.')
if(dir.exists(out)&&length(list.files(out,all.files=TRUE,no..=TRUE)))stop('Choose a new or empty output directory.')
pkgs<-c('splines','MASS');if('--figures'%in%args)pkgs<-c(pkgs,'ggplot2','patchwork')
missing<-pkgs[!vapply(pkgs,requireNamespace,quietly=TRUE,FUN.VALUE=logical(1))]
if(length(missing))stop('Required packages unavailable: ',paste(missing,collapse=', '),'. No installation attempted.')
c2<-file.path(raw,'C2_GB_AllMalignant_CaseListing_SEER17_2000-2023_labeled.txt')
p2<-file.path(raw,'P2_GB_AllMalignant_YearAge5_RaceEthnicity_SEER17_2000-2023.txt')
p3<-file.path(raw,'P3_GB_AllMalignant_YearAge5_Sex_RaceEthnicity_SEER17_2000-2023.txt')
stopifnot(all(file.exists(c(c2,p2,p3))))
dir.create(out,recursive=TRUE);inputs<-file.path(out,'private_inputs');logs<-file.path(out,'logs');dir.create(logs)
rscript<-file.path(R.home('bin'),if(.Platform$OS.type=='windows')'Rscript.exe' else 'Rscript')
run<-function(name,a){
 cat('Running ',name,'\n',sep='');flush.console()
 code<-system2(rscript,shQuote(c(file.path(project,'scripts',name),a)),stdout=file.path(logs,paste0(name,'.log')),stderr=file.path(logs,paste0(name,'.log')))
 if(code!=0L)stop('Analysis failed; inspect local log for ',name)
}
run('prepare_inputs.R',c(inputs,c2,p2))
run('paired_trends.R',c(out,c2,inputs))
run('uncertainty_checks.R',c(out,c2,inputs,file.path(project,'scripts/paired_trends.R')))
run('sex_analysis.R',c(out,p3,p2))
if('--figures'%in%args)run('figures.R',out)
cat('Analysis completed. Keep generated cells, results and logs private.\n')
