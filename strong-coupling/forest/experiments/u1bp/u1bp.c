/* Block-Poincare constant for 3D compact U(1) (Wilson action beta*sum cos theta_p), periodic L^3.
   varpi_n(f) = E_n(f)/Var f with E_n(f) = n^-3 sum_{cube positions} E Var(f | outside cube), cube = links with base in
   [x0,x0+n)^3 (all 3 directions). By translation invariance varpi = E Var(f|cube^c) / (n^3 * Var(f)/V).
   E Var(f|cube^c) is estimated by R inner heat-bath+overrelaxation samples of the cube links (outside frozen), unbiased
   sample variance; non-adjacent cubes (spacing n+2) are resampled simultaneously (conditionally independent).
   Observables (complex, phases at plaquette base point):
     0: sum sin th_12            k=0       3: sum cos th_p (all)  k=0
     1: sum sin th_12 e^{i k x1} k=2pi/L   4: sum cos th_12       k=0
     2: sum sin th_12 e^{i k x3} k=2pi/L
   Mass: slice sums along x3 of sum_{all p} cos th_p and of sin th_12; per-config output, analysis offline.
   usage: u1bp L beta n Nconf sep R s seed */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include <stdint.h>
static int L; static double beta; static double *th;
static uint64_t rs[2];
static inline uint64_t rotl(uint64_t x,int k){return (x<<k)|(x>>(64-k));}
static inline double U(void){uint64_t s0=rs[0],s1=rs[1],r=s0+s1;s1^=s0;rs[0]=rotl(s0,55)^s1^(s1<<14);rs[1]=rotl(s1,36);return ((r>>11)+0.5)*(1.0/9007199254740992.0);}
static inline int md(int a){a%=L; return a<0?a+L:a;}
static inline int S(int x,int y,int z){return (md(x)*L+md(y))*L+md(z);}
static inline double *ln(int x,int y,int z,int mu){return &th[S(x,y,z)*3+mu];}
static const int E[3][3]={{1,0,0},{0,1,0},{0,0,1}};
static double plaq(int x,int y,int z,int mu,int nu){
  const int*a=E[mu],*b=E[nu];
  return *ln(x,y,z,mu)+*ln(x+a[0],y+a[1],z+a[2],nu)-*ln(x+b[0],y+b[1],z+b[2],mu)-*ln(x,y,z,nu);}
static double vonmises(double k){
  if(k<1e-8) return M_PI*(2*U()-1);
  double tau=1+sqrt(1+4*k*k),rho=(tau-sqrt(2*tau))/(2*k),r=(1+rho*rho)/(2*rho),f;
  for(;;){double z=cos(M_PI*U()),u2=U();f=(1+r*z)/(r+z);double c=k*(r-f);
    if(c*(2-c)-u2>0||log(c/u2)+1-c>=0)break;}
  return (U()<0.5?-1:1)*acos(f>1?1:(f<-1?-1:f));}
static void update(int x,int y,int z,int mu,int over){
  double complex s=0; const int*a=E[mu];
  for(int nu=0;nu<3;nu++){ if(nu==mu)continue; const int*b=E[nu];
    s+=cexp(I*(*ln(x+a[0],y+a[1],z+a[2],nu)-*ln(x+b[0],y+b[1],z+b[2],mu)-*ln(x,y,z,nu)));
    s+=cexp(I*(-*ln(x+a[0]-b[0],y+a[1]-b[1],z+a[2]-b[2],nu)-*ln(x-b[0],y-b[1],z-b[2],mu)+*ln(x-b[0],y-b[1],z-b[2],nu)));}
  double m=-carg(s),k=beta*cabs(s),*t=ln(x,y,z,mu);
  if(over){ *t=2*m-*t; } else { *t=m+vonmises(k); }
  *t=remainder(*t,2*M_PI);}
static void sweep_box(int x0,int y0,int z0,int n,int over){
  for(int x=0;x<n;x++)for(int y=0;y<n;y++)for(int z=0;z<n;z++)for(int mu=0;mu<3;mu++)update(x0+x,y0+y,z0+z,mu,over);}
static void sweep_all(void){ for(int c=0;c<2;c++) sweep_box(0,0,0,L,c); }
#define NOBS 5
static void obs_box(int x0,int y0,int z0,int nb,double complex *o){
  double k=2*M_PI/L; for(int i=0;i<NOBS;i++)o[i]=0;
  for(int x=0;x<nb;x++)for(int y=0;y<nb;y++)for(int z=0;z<nb;z++){
    int X=x0+x,Y=y0+y,Z=z0+z; double p12=plaq(X,Y,Z,0,1), s12=sin(p12);
    o[0]+=s12; o[1]+=s12*cexp(I*k*md(X)); o[2]+=s12*cexp(I*k*md(Z));
    o[3]+=cos(p12)+cos(plaq(X,Y,Z,0,2))+cos(plaq(X,Y,Z,1,2)); o[4]+=cos(p12);}}
int main(int argc,char**argv){
  if(argc<9){fprintf(stderr,"usage\n");return 1;}
  L=atoi(argv[1]);beta=atof(argv[2]);int n=atoi(argv[3]),Nc=atoi(argv[4]),sep=atoi(argv[5]),R=atoi(argv[6]),s=atoi(argv[7]);
  rs[0]=0x9E3779B97F4A7C15ULL^strtoull(argv[8],0,10);rs[1]=0xD1B54A32D192ED03ULL;for(int i=0;i<20;i++)U();
  th=malloc(sizeof(double)*L*L*L*3); for(int i=0;i<L*L*L*3;i++)th[i]=M_PI*(2*U()-1);
  for(int i=0;i<400;i++)sweep_all();
  int sp=n+2, nc=(n>=L)?1:L/sp, nb=(n>=L)?L:n+1, off=(n>=L)?0:-1;
  int ncube=nc*nc*nc; double complex (*g)[NOBS]=malloc(sizeof(*g)*ncube*R);
  double *save=malloc(sizeof(double)*L*L*L*3);
  printf("# L=%d beta=%g n=%d Nconf=%d sep=%d R=%d s=%d ncube=%d\n",L,beta,n,Nc,sep,R,s,ncube);
  for(int c=0;c<Nc;c++){
    for(int i=0;i<sep;i++)sweep_all();
    double complex F[NOBS]; obs_box(0,0,0,L,F);
    double sl0[256],sl1[256];
    for(int z=0;z<L;z++){sl0[z]=sl1[z]=0;for(int x=0;x<L;x++)for(int y=0;y<L;y++){double p=plaq(x,y,z,0,1);
      sl1[z]+=sin(p); sl0[z]+=cos(p)+cos(plaq(x,y,z,0,2))+cos(plaq(x,y,z,1,2));}}
    for(int i=0;i<L*L*L*3;i++)save[i]=th[i];
    int sh[3]={(int)(U()*L),(int)(U()*L),(int)(U()*L)};   /* random global offset of the cube sublattice */
    for(int r=0;r<R;r++){
      if(r>0) for(int t=0;t<s;t++) for(int ov=0;ov<2;ov++)
        for(int q=0;q<ncube;q++){int a=q/(nc*nc),b=(q/nc)%nc,cc=q%nc; sweep_box(sh[0]+a*sp,sh[1]+b*sp,sh[2]+cc*sp,(n>=L)?L:n,ov);}
      for(int q=0;q<ncube;q++){int a=q/(nc*nc),b=(q/nc)%nc,cc=q%nc; obs_box(sh[0]+a*sp+off,sh[1]+b*sp+off,sh[2]+cc*sp+off,nb,g[q*R+r]);}}
    double num[NOBS]={0};
    for(int q=0;q<ncube;q++)for(int o=0;o<NOBS;o++){double complex m=0;for(int r=0;r<R;r++)m+=g[q*R+r][o];m/=R;
      double v=0;for(int r=0;r<R;r++){double complex d=g[q*R+r][o]-m;v+=creal(d*conj(d));} num[o]+=v/(R-1)/ncube;}
    for(int i=0;i<L*L*L*3;i++)th[i]=save[i];          /* outer chain continues from the unperturbed config */
    printf("C");for(int o=0;o<NOBS;o++)printf(" %.10g",num[o]);
    for(int o=0;o<NOBS;o++)printf(" %.10g %.10g",creal(F[o]),cimag(F[o]));
    for(int z=0;z<L;z++)printf(" %.8g",sl0[z]);for(int z=0;z<L;z++)printf(" %.8g",sl1[z]);printf("\n");fflush(stdout);}
  return 0;}
