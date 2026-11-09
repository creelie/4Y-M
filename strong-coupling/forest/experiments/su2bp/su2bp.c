/* Block-Poincare constant for 3D SU(2) Wilson action (beta/2) sum tr U_p, periodic L^3, quaternion links.
   Same estimator as u1bp.c: varpi = E Var(f|cube^c) / (n^3 Var(f)/V), cube = links with base in [x0,x0+n)^3,
   nested KP heat-bath + overrelaxation resampling of non-adjacent cubes (spacing n+4, so no w<=3 loop touches two).
   Observables (w = 1,2,3; W = (1/2) tr of w x w loop):  3(w-1)+0: sum_all-planes W, k=0
   3(w-1)+1: sum_{12} W e^{i k x1}, 3(w-1)+2: sum_{12} W e^{i k x3}, k = 2pi/L.
   usage: su2bp L beta n Nconf sep R s seed */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>
typedef struct{double a[4];} Q;
static int L; static double beta; static Q *U;
static uint64_t rs[2];
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline double R1(void){uint64_t s0=rs[0],s1=rs[1],r=s0+s1;s1^=s0;rs[0]=rotl(s0,55)^s1^(s1<<14);rs[1]=rotl(s1,36);return ((r>>11)+0.5)*(1.0/9007199254740992.0);}
static inline Q mul(Q x,Q y){Q z;double*a=x.a,*b=y.a;
  z.a[0]=a[0]*b[0]-a[1]*b[1]-a[2]*b[2]-a[3]*b[3];
  z.a[1]=a[0]*b[1]+b[0]*a[1]-(a[2]*b[3]-a[3]*b[2]);
  z.a[2]=a[0]*b[2]+b[0]*a[2]-(a[3]*b[1]-a[1]*b[3]);
  z.a[3]=a[0]*b[3]+b[0]*a[3]-(a[1]*b[2]-a[2]*b[1]); return z;}
static inline Q dag(Q x){x.a[1]=-x.a[1];x.a[2]=-x.a[2];x.a[3]=-x.a[3];return x;}
static inline int md(int a){a%=L;return a<0?a+L:a;}
static inline Q*lk(int x,int y,int z,int mu){return &U[((md(x)*L+md(y))*L+md(z))*3+mu];}
static const int E[3][3]={{1,0,0},{0,1,0},{0,0,1}};
static Q rnd_su2(void){Q q;double n;do{for(int i=0;i<4;i++)q.a[i]=2*R1()-1;n=q.a[0]*q.a[0]+q.a[1]*q.a[1]+q.a[2]*q.a[2]+q.a[3]*q.a[3];}while(n>1||n<1e-6);
  n=sqrt(n);for(int i=0;i<4;i++)q.a[i]/=n;return q;}
static void update(int x,int y,int z,int mu,int over){
  Q A={{0,0,0,0}}; const int*a=E[mu];
  for(int nu=0;nu<3;nu++){if(nu==mu)continue;const int*b=E[nu];
    Q t=mul(mul(*lk(x+a[0],y+a[1],z+a[2],nu),dag(*lk(x+b[0],y+b[1],z+b[2],mu))),dag(*lk(x,y,z,nu)));
    Q u=mul(mul(dag(*lk(x+a[0]-b[0],y+a[1]-b[1],z+a[2]-b[2],nu)),dag(*lk(x-b[0],y-b[1],z-b[2],mu))),*lk(x-b[0],y-b[1],z-b[2],nu));
    for(int i=0;i<4;i++)A.a[i]+=t.a[i]+u.a[i];}
  double k=sqrt(A.a[0]*A.a[0]+A.a[1]*A.a[1]+A.a[2]*A.a[2]+A.a[3]*A.a[3]); Q*l=lk(x,y,z,mu);
  if(k<1e-12){*l=rnd_su2();return;}
  Q V=A; for(int i=0;i<4;i++)V.a[i]/=k;          /* weight exp(beta k (1/2) tr(U V)) */
  if(over){ *l=mul(mul(dag(V),dag(*l)),dag(V)); }
  else { double al=beta*k,w0,x1;
    for(;;){double r1=R1(),r2=R1(),r3=R1(),r4=R1(),c=cos(2*M_PI*r2);
      x1=-(log(r1)+c*c*log(r3))/al; if(r4*r4<=1-0.5*x1)break;}
    w0=1-x1; double rr=sqrt(fmax(0,1-w0*w0)),ct=2*R1()-1,ph=2*M_PI*R1(),st=sqrt(1-ct*ct);
    Q W={{w0,rr*st*cos(ph),rr*st*sin(ph),rr*ct}}; *l=mul(W,dag(V)); }
  double n=sqrt(l->a[0]*l->a[0]+l->a[1]*l->a[1]+l->a[2]*l->a[2]+l->a[3]*l->a[3]); for(int i=0;i<4;i++)l->a[i]/=n;}
static void sweep_box(int x0,int y0,int z0,int n,int over){
  for(int x=0;x<n;x++)for(int y=0;y<n;y++)for(int z=0;z<n;z++)for(int mu=0;mu<3;mu++)update(x0+x,y0+y,z0+z,mu,over);}
static void sweep_all(void){for(int c=0;c<2;c++)sweep_box(0,0,0,L,c);}
static double loop(int x,int y,int z,int mu,int nu,int w){
  Q P={{1,0,0,0}}; const int*a=E[mu],*b=E[nu]; int p[3]={x,y,z};
  for(int i=0;i<w;i++){P=mul(P,*lk(p[0],p[1],p[2],mu));for(int j=0;j<3;j++)p[j]+=a[j];}
  for(int i=0;i<w;i++){P=mul(P,*lk(p[0],p[1],p[2],nu));for(int j=0;j<3;j++)p[j]+=b[j];}
  for(int i=0;i<w;i++){for(int j=0;j<3;j++)p[j]-=a[j];P=mul(P,dag(*lk(p[0],p[1],p[2],mu)));}
  for(int i=0;i<w;i++){for(int j=0;j<3;j++)p[j]-=b[j];P=mul(P,dag(*lk(p[0],p[1],p[2],nu)));}
  return P.a[0];}
#define NOBS 9
static void obs_box(int x0,int y0,int z0,int nb,double complex*o){
  double k=2*M_PI/L; for(int i=0;i<NOBS;i++)o[i]=0;
  for(int x=0;x<nb;x++)for(int y=0;y<nb;y++)for(int z=0;z<nb;z++){int X=x0+x,Y=y0+y,Z=z0+z;
    for(int w=1;w<=3;w++){double w12=loop(X,Y,Z,0,1,w);
      o[3*(w-1)]+=w12+loop(X,Y,Z,0,2,w)+loop(X,Y,Z,1,2,w);
      o[3*(w-1)+1]+=w12*cexp(I*k*md(X)); o[3*(w-1)+2]+=w12*cexp(I*k*md(Z));}}}
int main(int argc,char**argv){
  if(argc<9){fprintf(stderr,"usage\n");return 1;}
  L=atoi(argv[1]);beta=atof(argv[2]);int n=atoi(argv[3]),Nc=atoi(argv[4]),sep=atoi(argv[5]),R=atoi(argv[6]),s=atoi(argv[7]);
  rs[0]=0x9E3779B97F4A7C15ULL^strtoull(argv[8],0,10);rs[1]=0xD1B54A32D192ED03ULL;for(int i=0;i<20;i++)R1();
  int NL=L*L*L*3; U=malloc(sizeof(Q)*NL); for(int i=0;i<NL;i++)U[i]=rnd_su2();
  for(int i=0;i<300;i++)sweep_all();
  int full=(n>=L), sp=n+4, nc=full?1:L/sp, nb=full?L:n+3, off=full?0:-3, ncube=nc*nc*nc;
  double complex (*g)[NOBS]=malloc(sizeof(*g)*ncube*R); Q*save=malloc(sizeof(Q)*NL);
  printf("# L=%d beta=%g n=%d Nconf=%d sep=%d R=%d s=%d ncube=%d\n",L,beta,n,Nc,sep,R,s,ncube);
  for(int c=0;c<Nc;c++){
    for(int i=0;i<sep;i++)sweep_all();
    double complex F[NOBS]; obs_box(0,0,0,L,F);
    for(int i=0;i<NL;i++)save[i]=U[i];
    int sh[3]={(int)(R1()*L),(int)(R1()*L),(int)(R1()*L)};
    for(int r=0;r<R;r++){
      if(r>0)for(int t=0;t<s;t++)for(int ov=0;ov<2;ov++)
        for(int q=0;q<ncube;q++){int a=q/(nc*nc),b=(q/nc)%nc,cc=q%nc;sweep_box(sh[0]+a*sp,sh[1]+b*sp,sh[2]+cc*sp,full?L:n,ov);}
      for(int q=0;q<ncube;q++){int a=q/(nc*nc),b=(q/nc)%nc,cc=q%nc;obs_box(sh[0]+a*sp+off,sh[1]+b*sp+off,sh[2]+cc*sp+off,nb,g[q*R+r]);}}
    double num[NOBS]={0};
    for(int q=0;q<ncube;q++)for(int o=0;o<NOBS;o++){double complex m=0;for(int r=0;r<R;r++)m+=g[q*R+r][o];m/=R;
      double v=0;for(int r=0;r<R;r++){double complex d=g[q*R+r][o]-m;v+=creal(d*conj(d));}num[o]+=v/(R-1)/ncube;}
    for(int i=0;i<NL;i++)U[i]=save[i];
    printf("C");for(int o=0;o<NOBS;o++)printf(" %.10g",num[o]);
    for(int o=0;o<NOBS;o++)printf(" %.10g %.10g",creal(F[o]),cimag(F[o]));printf("\n");fflush(stdout);}
  return 0;}
