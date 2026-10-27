# Entirely synthetic tests: no study data or external files are read.
stopifnot(abs(sum(c(.2,.3,.5))-1)<1e-12)
w<-c(.2,.3,.5);p0<-c(.2,.4,.6);p1<-c(.3,.5,.7);r0<-c(1,3,4);r1<-c(2,2,5)
C<-sum(w*(p1-p0)*(r0+r1)/2);W<-sum(w*(p1+p0)*(r1-r0)/2)
stopifnot(abs(C+W-sum(w*(p1*r1-p0*r0)))<1e-12)
n<-c(0,2,4);a<-c(.1,.2,.3);y<-sum(a*n);v<-sum(a*a*n);u<-max(a)
lo<-qgamma(.025,y*y/v,scale=v/y);hi<-qgamma(.975,(y+u)^2/(v+u*u),scale=(v+u*u)/(y+u))
stopifnot(lo>=0,hi>lo,is.finite(hi))
slope<-function(r,v,years){w<-r*r/v;x<-years-mean(years);y<-log(r);(sum(w*x*y)-sum(w*x)*sum(w*y)/sum(w))/(sum(w*x*x)-sum(w*x)^2/sum(w))}
yr<-2000:2023;r<-exp(.02*(yr-2000));v<-r/100
stopifnot(abs(slope(r,v,yr)-coef(lm(log(r)~yr,weights=r*r/v))[2])<1e-12)
obs<--.75;fit<--.70;draw<-seq(-1,-.4,length.out=100)
q<-quantile(draw-fit,c(.025,.975));ci<-obs-rev(q)
stopifnot(ci[1]<ci[2],abs(mean(ci)-obs)<1e-12)
cat('Synthetic mathematical checks passed.\n')
