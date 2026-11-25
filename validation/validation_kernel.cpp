#include <cmath>
#include <algorithm>
#include <vector>
#include <limits>

static double evtime(double s,double target,double eq,double rate){
 if(rate<=0 || std::abs(s-eq)<1e-300)return INFINITY;
 double r=(target-eq)/(s-eq);
 return r>0 && r<1 ? -std::log(r)/rate : INFINITY;
}
static double source(double B,double c,double al,double be,double dt,const double *p,int &holds,int &hits){
 const double m=p[0],chi=p[6],lo=p[7],mid=p[8],cap=p[9];
 double g=al*c,rate=be/m,eq=be>0?g/be:0;
 if(p[12]<.5)return be>0?eq+(B-eq)*std::exp(-rate*dt):B+g*dt/m;
 double b=B,rem=dt;
 for(int it=0;it<8 && rem>1e-15;it++){
  if(b>=cap)return cap;
  if(b==lo && g-be*lo<=0){holds++;return b;}
  if(b<lo){
   double v=chi*g/m,event=v>0?(lo-b)/v:INFINITY,step=std::min(rem,event);
   b+=v*step;if(event<=rem)b=lo;rem-=step;
  }else if(b<mid || (b==mid && g-be*mid<=0)){
   double q=g-be*b,target=q>=0?mid:lo;
   double event=be>0?evtime(b,target,eq,rate):(g>0?(target-b)*m/g:INFINITY);
   double step=std::min(rem,event);
   b=be>0?eq+(b-eq)*std::exp(-rate*step):b+g*step/m;
   if(event<=rem)b=target;rem-=step;
  }else{
   double eq2=eq*cap,target=g*cap-be*b*b>=0?cap:mid;
   double event=be>0?evtime(b*b,target*target,eq2,2*rate):(g>0?(target*target-b*b)*m/(2*g*cap):INFINITY);
   double step=std::min(rem,event);
   double b2=be>0?eq2+(b*b-eq2)*std::exp(-2*rate*step):b*b+2*g*cap*step/m;
   b=std::sqrt(std::max(b2,0.));
   if(event<=rem){b=target;if(target==cap)hits++;}rem-=step;
  }
 }
 return rem>1e-12?NAN:b;
}

extern "C" int filtration_run(const double *p,int n,const double *samples,int ns,int axis,double cin,double *out,double *profile_a,double *profile_b){
 // out columns: t*,PV,u*,dp*,meanB,inventory,input,output,maxmasserror,
 // minimum_all_saturations,minimum_S1,maxc,keff/k0,hydraulicconductance,steps,holds,caphits
 const int nc=17;const double m=p[0],eps=p[1],wr=1/(1-eps),dx=1./n,cfl=p[13],dtmax=p[14];
 std::vector<double>a(n,0),b(n,0),phi(n),c(n),k(n),mu(n),an(n),bn(n),cs(n);
 double t=0,pv=0,in=0,ou=0,maxerr=0,minsat=0,minS1=1,maxc=0;
 int steps=0,holds=0,hits=0,j=0;
 while(j<ns){
  double resistance=0,invk=0,minphi=1,meanB=0,inventory=0;
  for(int i=0;i<n;i++){
   phi[i]=1-wr*b[i];if(phi[i]<=0)return -2;
   c[i]=a[i]/phi[i];k[i]=std::pow(1-std::sqrt(std::max(wr*b[i],0.)),3);mu[i]=1+2.5*c[i];
   if(k[i]<=0)return -3;resistance+=mu[i]/k[i]/n;invk+=1/k[i]/n;
   minphi=std::min(minphi,phi[i]);meanB+=b[i]/n;inventory+=(a[i]+b[i])/n;
  }
  double u=p[11]<.5?p[10]/resistance:p[10],dp=u*resistance;
  double current=axis==0?t:pv,target=samples[j];
  if(current>=target-1e-12*std::max(1.,target)){
   double *r=out+j*nc;
   r[0]=t;r[1]=pv;r[2]=u;r[3]=dp;r[4]=meanB;r[5]=inventory;r[6]=in;r[7]=ou;r[8]=maxerr;
   r[9]=minsat;r[10]=minS1;r[11]=maxc;r[12]=1/invk;r[13]=1/resistance;r[14]=steps;r[15]=holds;r[16]=hits;
   for(int i=0;i<n;i++){profile_a[j*n+i]=a[i];profile_b[j*n+i]=b[i];}j++;continue;
  }
  if(steps>20000000)return -4;
  if(u<=0 && axis==1)return -5;
  double boost=p[12]<.5?1:std::max(1.,std::max(p[6],p[9]/p[8]));
  double alpha=p[2]+p[3]*u;
  double dt=std::min(cfl*m*minphi*dx/std::max(u,1e-300),cfl*m*minphi/(alpha*boost+1e-300));
  dt=std::min(dt,dtmax);dt=std::min(dt,axis==0?target-t:(target-pv)*m/u);
  bool accepted=false;int step_holds=0,step_hits=0;
  for(int retry=0;retry<30 && !accepted;retry++){
   double up=cin;
   for(int i=0;i<n;i++){an[i]=a[i]-dt*u*(c[i]-up)/(m*dx);up=c[i];cs[i]=an[i]/phi[i];}
   // As in the original routine, dynamic coefficients use the post-transport hydraulic state.
   double us=u;
   if(p[3]!=0 || p[5]!=0){
    double res2=0;for(int i=0;i<n;i++)res2+=(1+2.5*cs[i])/k[i]/n;
    us=p[11]<.5?p[10]/res2:p[10];
   }
   step_holds=0;step_hits=0;accepted=true;
   for(int i=0;i<n;i++){
    double al=p[2]+p[3]*us,be=p[4]+p[5]*us*(1+2.5*cs[i])/k[i];
    bn[i]=source(b[i],cs[i],al,be,dt,p,step_holds,step_hits);
    double mobile=an[i]-(bn[i]-b[i]);
    if(!std::isfinite(bn[i]) || mobile<-1e-13 || bn[i]<-1e-13 || 1-wr*bn[i]-mobile<-1e-12){accepted=false;break;}
   }
   if(!accepted)dt*=.5;
  }
  if(!accepted)return -6;
  in+=dt*u*cin/m;ou+=dt*u*c[n-1]/m;
  t+=dt;pv+=dt*u/m;steps++;holds+=step_holds;hits+=step_hits;
  double inventory_now=0;
  for(int i=0;i<n;i++){
   a[i]=an[i]-(bn[i]-b[i]);b[i]=bn[i];inventory_now+=(a[i]+b[i])/n;
   double ph=1-wr*b[i],carrier=ph-a[i];
   minS1=std::min(minS1,carrier);minsat=std::min(minsat,std::min(a[i],std::min(b[i],carrier)));maxc=std::max(maxc,a[i]/ph);
  }
  maxerr=std::max(maxerr,std::abs(inventory_now-in+ou));
 }
 return 0;
}
