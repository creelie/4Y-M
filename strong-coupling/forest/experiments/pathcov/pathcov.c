/* pathcov.c -- posterior covariances of path scores (Hexagon v7, Thm 164 / Prop 168).
 *
 * Fine torus Lambda_N, N = b L, d = 4, SU(2) Wilson exp(beta * sum_p Re U_p)  (Re = 1/2 ReTr).
 * Straight-path map: coarse link (y,mu) = product of the b fine links (b y + i e_mu, mu), i = 0..b-1.
 * Interiors I (1 <= i <= b-2) are integrated out exactly (Lemma bppathint):
 *   pi_V(W) ~ exp(beta sum_{p in P0} Re U_p) * prod_{I} Z(beta|S_I|) * prod_e p_e(V_e | W)
 * with P0 = plaquettes containing no interior link, Z(B) = 2 I_1(B)/B and
 *   p_e(V|W) = F_e(Re(R^-1 V T)),  T = R'^-1 W_m ... W_1,  W_i = S_i/|S_i|,
 *   F_e(c) = sum_n (n+1) prod_i a_n(B_i) U_n(c),  a_n(B) = I_{n+1}(B)/I_1(B).
 * Score a_e^X = D_X log p_e(V_e e^{tX}), X = q/sqrt2, q in {i,j,k}:  a_e = -(F'/F) * vec(T G0 T^-1)/sqrt2.
 *
 * Usage: pathcov b L beta Vmode(id|haar) ntherm nmeas every seed outprefix
 *        pathcov selftest
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>

typedef struct { double w, x, y, z; } q_t;
static inline q_t qm(q_t a, q_t b) {
  q_t r = { a.w*b.w - a.x*b.x - a.y*b.y - a.z*b.z,
            a.w*b.x + a.x*b.w + a.y*b.z - a.z*b.y,
            a.w*b.y - a.x*b.z + a.y*b.w + a.z*b.x,
            a.w*b.z + a.x*b.y - a.y*b.x + a.z*b.w };
  return r;
}
static inline q_t qc(q_t a) { q_t r = { a.w, -a.x, -a.y, -a.z }; return r; }
static inline q_t qadd(q_t a, q_t b) { q_t r = { a.w+b.w, a.x+b.x, a.y+b.y, a.z+b.z }; return r; }
static inline double qn(q_t a) { return sqrt(a.w*a.w + a.x*a.x + a.y*a.y + a.z*a.z); }
static inline q_t qs(q_t a, double s) { q_t r = { a.w*s, a.x*s, a.y*s, a.z*s }; return r; }
static const q_t Q1 = {1,0,0,0};

/* xoshiro256** */
static uint64_t rs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t rnext(void) {
  uint64_t r = rotl(rs[1] * 5, 7) * 9, t = rs[1] << 17;
  rs[2] ^= rs[0]; rs[3] ^= rs[1]; rs[1] ^= rs[2]; rs[0] ^= rs[3]; rs[2] ^= t; rs[3] = rotl(rs[3], 45);
  return r;
}
static inline double ru(void) { return (rnext() >> 11) * (1.0 / 9007199254740992.0); }
static void rseed(uint64_t s) { for (int i = 0; i < 4; i++) { s += 0x9E3779B97F4A7C15ULL; uint64_t z = s;
  z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL; z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL; rs[i] = z ^ (z >> 31); } }
static double rn(void) { double u = ru(), v = ru(); if (u < 1e-300) u = 1e-300; return sqrt(-2*log(u))*cos(2*M_PI*v); }
static q_t rhaar(void) { q_t q = { rn(), rn(), rn(), rn() }; return qs(q, 1.0/qn(q)); }
/* sample g with density ~ exp(B Re g) on S^3 */
static q_t rvmf(double B) {
  if (B < 1e-10) return rhaar();
  double w;
  for (;;) {
    double u = ru();
    w = 1.0 + log(u + (1.0 - u) * exp(-2.0*B)) / B;
    if (w < -1) w = -1;
    if (ru() * ru() < 1.0 - w*w) { /* accept with prob sqrt(1-w^2): compare u^2 < 1-w^2 */ break; }
  }
  double s = sqrt(fmax(0.0, 1.0 - w*w));
  double x = rn(), y = rn(), z = rn(), r = sqrt(x*x + y*y + z*z);
  q_t q = { w, s*x/r, s*y/r, s*z/r };
  return q;
}

/* ---------- Bessel pieces ---------- */
static int NMAX;
static double logZfun(double B) { /* log(2 I_1(B)/B) */
  if (B < 1e-12) return 0.0;
  double h = B/2, t = 1.0, s = 1.0; /* I_1(B)/(B/2) = sum_k h^{2k}/(k!(k+1)!) */
  for (int k = 1; k < 200; k++) { t *= h*h / (k*(k+1.0)); s += t; if (t < 1e-17*s) break; }
  return log(s);
}
static void acoef(double B, double *a) { /* a[n] = I_{n+1}(B)/I_1(B), n = 0..NMAX */
  a[0] = 1.0;
  if (B < 1e-12) { for (int n = 1; n <= NMAX; n++) a[n] = 0.0; return; }
  int K = NMAX + 60 + (int)(2*B);
  static double *rr = 0; static int rrK = 0;
  if (K + 2 > rrK) { rrK = K + 2 + 256; rr = realloc(rr, sizeof(double) * rrK); }
  double r = 0.0;
  for (int k = K; k >= 1; k--) { r = 1.0 / (2.0*(k+1)/B + r); rr[k] = r; } /* rr[k] = I_{k+1}/I_k */
  for (int n = 1; n <= NMAX; n++) { a[n] = a[n-1] * rr[n]; if (a[n] < 1e-300) a[n] = 0.0; }
}
/* F(c) and F'(c) with weights P[n] */
static void Fval(const double *P, double c, double *F, double *dF) {
  double u0 = 1, u1 = 2*c, d0 = 0, d1 = 2, f = P[0], df = 0;
  if (NMAX >= 1) { f += 2*P[1]*u1; df += 2*P[1]*d1; }
  for (int n = 2; n <= NMAX; n++) {
    double u2 = 2*c*u1 - u0, d2 = 2*u1 + 2*c*d1 - d0;
    f += (n+1)*P[n]*u2; df += (n+1)*P[n]*d2;
    u0 = u1; u1 = u2; d0 = d1; d1 = d2;
    if (P[n] == 0.0) break;
  }
  *F = f; *dF = df;
}

/* ---------- lattice ---------- */
static int b, L, N, D = 4, m;
static long V4, NL;
static double beta;
static q_t *U;
static int *up, *dn;           /* up[s*4+mu] */
static signed char *ltype;     /* 0 free, 1 W-nonfree, 2 interior */
static int *intid;             /* link -> interior id or -1 */
static long *intlink;          /* interior id -> link */
static int nint;
static q_t *Sint; static double *Bint, *aint, *logZint; /* per interior */
static long npath;
static long *pR, *pRp; static int *pint; /* path -> R link, R' link, interiors (m each) */
static q_t *Vc; static double *logF;
#define CAP 8
static int *affI, *naffI;     /* per link */
static long *affP; static int *naffP;
static unsigned short *p0mask;

static inline long site(const int *x) { return x[0] + (long)N*(x[1] + (long)N*(x[2] + (long)N*x[3])); }
static void coords(long s, int *x) { for (int i = 0; i < 4; i++) { x[i] = s % N; s /= N; } }

/* staple of plaquette k (0..11) for link (s,mu): k = 2*(nu index) + dir */
static q_t staple1(long s, int mu, int nu, int fwd) {
  if (fwd) {
    long smu = up[s*4+mu], snu = up[s*4+nu];
    return qm(qm(U[smu*4+nu], qc(U[snu*4+mu])), qc(U[s*4+nu]));
  } else {
    long sm = dn[s*4+nu], smm = up[sm*4+mu];
    return qm(qm(qc(U[smm*4+nu]), qc(U[sm*4+mu])), U[sm*4+nu]);
  }
}
/* the three other links of plaquette */
static void plaqlinks(long s, int mu, int nu, int fwd, long *o) {
  if (fwd) { long smu = up[s*4+mu], snu = up[s*4+nu]; o[0] = smu*4+nu; o[1] = snu*4+mu; o[2] = s*4+nu; }
  else { long sm = dn[s*4+nu], smm = up[sm*4+mu]; o[0] = smm*4+nu; o[1] = sm*4+mu; o[2] = sm*4+nu; }
}
static q_t fullstaple(long l) {
  long s = l / 4; int mu = l % 4; q_t S = {0,0,0,0};
  for (int nu = 0; nu < 4; nu++) if (nu != mu) { S = qadd(S, staple1(s, mu, nu, 1)); S = qadd(S, staple1(s, mu, nu, 0)); }
  return S;
}
static q_t p0staple(long l) {
  long s = l / 4; int mu = l % 4; q_t S = {0,0,0,0}; int k = 0;
  for (int nu = 0; nu < 4; nu++) if (nu != mu) for (int f = 1; f >= 0; f--, k++)
    if (p0mask[l] >> k & 1) S = qadd(S, staple1(s, mu, nu, f));
  return S;
}

static void update_interior(int I) {
  q_t S = fullstaple(intlink[I]);
  Sint[I] = S; double n = qn(S); Bint[I] = beta * n;
  acoef(Bint[I], aint + (long)I*(NMAX+1)); logZint[I] = logZfun(Bint[I]);
}
static q_t Tpath(long e) { /* T = R'^-1 W_m ... W_1 */
  q_t T = qc(U[pRp[e]]);
  for (int i = m - 1; i >= 0; i--) {
    int I = pint[e*m+i]; double n = qn(Sint[I]);
    q_t W = n > 1e-14 ? qs(Sint[I], 1.0/n) : Q1;
    T = qm(T, W);
  }
  return T;
}
static double *Pbuf;
static void pathP(long e) {
  for (int n = 0; n <= NMAX; n++) { double p = 1; for (int i = 0; i < m; i++) p *= aint[(long)pint[e*m+i]*(NMAX+1)+n]; Pbuf[n] = p; }
}
static int badF = 0;
static double pathlogF(long e) {
  pathP(e);
  q_t T = Tpath(e), G = qm(qm(qc(U[pR[e]]), Vc[e]), T);
  double c = fmax(-1.0, fmin(1.0, G.w)), F, dF; Fval(Pbuf, c, &F, &dF);
  if (!(F > 0)) { badF++; return -1e300; }
  return log(F);
}
static void pathscore(long e, float *out) {
  pathP(e);
  q_t T = Tpath(e), G0 = qm(qm(qc(U[pR[e]]), Vc[e]), T);
  double c = fmax(-1.0, fmin(1.0, G0.w)), F, dF; Fval(Pbuf, c, &F, &dF);
  q_t H = qm(qm(T, G0), qc(T));
  double f = dF / F / sqrt(2.0);
  out[0] = (float)(-f*H.x); out[1] = (float)(-f*H.y); out[2] = (float)(-f*H.z);
}

static void build(void) {
  N = b * L; V4 = (long)N*N*N*N; NL = V4 * 4; m = b - 2;
  U = malloc(sizeof(q_t) * NL); up = malloc(sizeof(int) * NL); dn = malloc(sizeof(int) * NL);
  ltype = calloc(NL, 1); intid = malloc(sizeof(int) * NL); p0mask = malloc(sizeof(short) * NL);
  for (long s = 0; s < V4; s++) { int x[4]; coords(s, x);
    for (int mu = 0; mu < 4; mu++) { int y[4]; memcpy(y, x, sizeof y); y[mu] = (x[mu]+1)%N; up[s*4+mu] = site(y);
      y[mu] = (x[mu]+N-1)%N; dn[s*4+mu] = site(y); } }
  for (long l = 0; l < NL; l++) U[l] = Q1;
  npath = (long)L*L*L*L*4; nint = 0;
  pR = malloc(sizeof(long)*npath); pRp = malloc(sizeof(long)*npath); pint = malloc(sizeof(int)*npath*m);
  Vc = malloc(sizeof(q_t)*npath); logF = malloc(sizeof(double)*npath);
  long *linkpath = malloc(sizeof(long)*NL); for (long l = 0; l < NL; l++) { intid[l] = -1; linkpath[l] = -1; }
  intlink = malloc(sizeof(long)*npath*m);
  for (long y = 0; y < (long)L*L*L*L; y++) {
    int c[4]; long t = y; for (int i = 0; i < 4; i++) { c[i] = t % L; t /= L; }
    for (int mu = 0; mu < 4; mu++) { long e = y*4+mu; int x[4];
      for (int i = 0; i < 4; i++) x[i] = b*c[i];
      for (int i = 0; i < b; i++) { x[mu] = b*c[mu] + i; long l = site(x)*4+mu; linkpath[l] = e;
        if (i == 0) pR[e] = l; else if (i == b-1) pRp[e] = l;
        else { intid[l] = nint; intlink[nint] = l; pint[e*m + i-1] = nint; ltype[l] = 2; nint++; } } }
  }
  affI = malloc(sizeof(int)*NL*CAP); naffI = calloc(NL, sizeof(int));
  affP = malloc(sizeof(long)*NL*CAP); naffP = calloc(NL, sizeof(int));
  for (long l = 0; l < NL; l++) {
    p0mask[l] = 0; if (ltype[l] == 2) continue;
    long s = l/4; int mu = l%4, k = 0;
    for (int nu = 0; nu < 4; nu++) if (nu != mu) for (int f = 1; f >= 0; f--, k++) {
      long o[3]; plaqlinks(s, mu, nu, f, o); int hasI = -1;
      for (int j = 0; j < 3; j++) if (ltype[o[j]] == 2) hasI = intid[o[j]];
      if (hasI < 0) p0mask[l] |= 1 << k;
      else { int dup = 0; for (int j = 0; j < naffI[l]; j++) if (affI[l*CAP+j] == hasI) dup = 1;
        if (!dup) { if (naffI[l] >= CAP) { fprintf(stderr, "CAP I\n"); exit(1); } affI[l*CAP + naffI[l]++] = hasI; } }
    }
    /* affected paths: own path if endpoint, and paths of affected interiors */
    long cand[CAP+1]; int nc = 0;
    if (linkpath[l] >= 0) cand[nc++] = linkpath[l];
    for (int j = 0; j < naffI[l]; j++) cand[nc++] = linkpath[intlink[affI[l*CAP+j]]];
    for (int j = 0; j < nc; j++) { int dup = 0; for (int k2 = 0; k2 < naffP[l]; k2++) if (affP[l*CAP+k2] == cand[j]) dup = 1;
      if (!dup) { if (naffP[l] >= CAP) { fprintf(stderr, "CAP P\n"); exit(1); } affP[l*CAP + naffP[l]++] = cand[j]; } }
    ltype[l] = (naffI[l] == 0 && linkpath[l] < 0) ? 0 : 1;
  }
  free(linkpath);
  Sint = malloc(sizeof(q_t)*nint); Bint = malloc(sizeof(double)*nint); logZint = malloc(sizeof(double)*nint);
  aint = malloc(sizeof(double)*(long)nint*(NMAX+1)); Pbuf = malloc(sizeof(double)*(NMAX+1));
}

/* ---------- updates ---------- */
static long nacc = 0, ntry = 0;
static double eps = 0.3;
static void heatbath_free(long l) {
  q_t S = fullstaple(l); double n = qn(S);
  if (n < 1e-14) { U[l] = rhaar(); return; }
  q_t Sh = qs(S, 1.0/n);
  U[l] = qm(rvmf(beta*n), qc(Sh));
}
static void overrelax_free(long l) {
  q_t S = fullstaple(l); double n = qn(S); if (n < 1e-14) return;
  q_t Sb = qs(qc(S), 1.0/n);
  U[l] = qm(qm(Sb, qc(U[l])), Sb);
}
static q_t bakS[CAP]; static double bakB[CAP], bakZ[CAP], bakF[CAP], *baka;
static void metro_nonfree(long l) {
  q_t S0 = p0staple(l);
  q_t h = { 1.0, eps*rn(), eps*rn(), eps*rn() }; h = qs(h, 1.0/qn(h));
  if (ru() < 0.5) h = qc(h);
  q_t Uold = U[l], Unew = qm(h, Uold);
  double dlog = beta * (qm(Unew, S0).w - qm(Uold, S0).w);
  int nI = naffI[l], nP = naffP[l];
  for (int j = 0; j < nI; j++) { int I = affI[l*CAP+j]; bakS[j] = Sint[I]; bakB[j] = Bint[I]; bakZ[j] = logZint[I];
    memcpy(baka + (long)j*(NMAX+1), aint + (long)I*(NMAX+1), sizeof(double)*(NMAX+1)); }
  U[l] = Unew;
  for (int j = 0; j < nI; j++) { int I = affI[l*CAP+j]; update_interior(I); dlog += logZint[I] - bakZ[j]; }
  for (int j = 0; j < nP; j++) { long e = affP[l*CAP+j]; bakF[j] = logF[e]; logF[e] = pathlogF(e); dlog += logF[e] - bakF[j]; }
  ntry++;
  if (dlog >= 0 || ru() < exp(dlog)) { nacc++; return; }
  U[l] = Uold;
  for (int j = 0; j < nI; j++) { int I = affI[l*CAP+j]; Sint[I] = bakS[j]; Bint[I] = bakB[j]; logZint[I] = bakZ[j];
    memcpy(aint + (long)I*(NMAX+1), baka + (long)j*(NMAX+1), sizeof(double)*(NMAX+1)); }
  for (int j = 0; j < nP; j++) logF[affP[l*CAP+j]] = bakF[j];
}
static void sweep(void) {
  for (long l = 0; l < NL; l++) { if (ltype[l] == 0) heatbath_free(l); else if (ltype[l] == 1) for (int h = 0; h < 3; h++) metro_nonfree(l); }
  for (int r = 0; r < 2; r++) for (long l = 0; l < NL; l++) if (ltype[l] == 0) overrelax_free(l);
}
static double plaqP0(void) { /* mean Re U_p over P0 plaquettes, via forward plaquettes */
  double s = 0; long c = 0;
  for (long l = 0; l < NL; l++) { if (ltype[l] == 2) continue; long x = l/4; int mu = l%4, k = 0;
    for (int nu = 0; nu < 4; nu++) if (nu != mu) for (int f = 1; f >= 0; f--, k++)
      if (f == 1 && nu > mu && (p0mask[l] >> k & 1)) { s += qm(U[l], staple1(x, mu, nu, 1)).w; c++; } }
  return s / c;
}
static void refresh_all(void) {
  for (int I = 0; I < nint; I++) update_interior(I);
  for (long e = 0; e < npath; e++) logF[e] = pathlogF(e);
}

/* ---------- self test ---------- */
static int selftest(void) {
  rseed(12345); NMAX = 120; Pbuf = malloc(sizeof(double)*(NMAX+1));
  int ok = 1;
  /* logZ vs quadrature */
  double Bs[] = {0.3, 2.0, 7.5, 14.0};
  for (int k = 0; k < 4; k++) { double B = Bs[k], s = 0; int n = 200000;
    for (int i = 0; i < n; i++) { double th = M_PI*(i+0.5)/n; s += exp(B*cos(th))*sin(th)*sin(th); }
    s *= 2.0/M_PI * M_PI/n; double d = fabs(log(s) - logZfun(B));
    printf("logZ B=%.2f quad=%.10f series=%.10f diff=%.2e\n", B, log(s), logZfun(B), d); if (d > 1e-8) ok = 0; }
  /* product density: m factors, compare E[Re(V Q)] direct vs Haar-weighted, and score vs finite difference */
  for (int mm = 2; mm <= 6; mm += 4) for (int k = 0; k < 2; k++) {
    double Bmax = k == 0 ? 1.5 : 9.0; q_t R = rhaar(), Rp = rhaar(), W[8]; double B[8], *a = malloc(sizeof(double)*(NMAX+1)*8);
    for (int i = 0; i < mm; i++) { W[i] = rhaar(); B[i] = Bmax*ru(); acoef(B[i], a + i*(NMAX+1)); }
    for (int n = 0; n <= NMAX; n++) { double p = 1; for (int i = 0; i < mm; i++) p *= a[i*(NMAX+1)+n]; Pbuf[n] = p; }
    q_t T = qc(Rp); for (int i = mm-1; i >= 0; i--) T = qm(T, W[i]);
    q_t Q = rhaar(); long ns = 4000000; double e1 = 0, e2 = 0, e1s = 0, e2s = 0;
    for (long t = 0; t < ns; t++) {
      q_t P = R; for (int i = 0; i < mm; i++) P = qm(P, qm(rvmf(B[i]), qc(W[i]))); P = qm(P, Rp);
      double v = qm(P, Q).w; e1 += v; e1s += v*v;
      q_t Vh = rhaar(); q_t G = qm(qm(qc(R), Vh), T); double F, dF; Fval(Pbuf, G.w, &F, &dF);
      double v2 = qm(Vh, Q).w * F; e2 += v2; e2s += v2*v2; }
    e1 /= ns; e2 /= ns; double s1 = sqrt((e1s/ns - e1*e1)/ns), s2 = sqrt((e2s/ns - e2*e2)/ns);
    double z = fabs(e1-e2)/sqrt(s1*s1+s2*s2);
    printf("m=%d Bmax=%.1f E_direct=%.5f+-%.5f E_char=%.5f+-%.5f z=%.2f\n", mm, Bmax, e1, s1, e2, s2, z); if (z > 4) ok = 0;
    /* score finite difference */
    q_t Vt = rhaar(); q_t G0 = qm(qm(qc(R), Vt), T); double F, dF; Fval(Pbuf, G0.w, &F, &dF);
    q_t H = qm(qm(T, G0), qc(T)); double f = dF/F/sqrt(2.0); double an[3] = {-f*H.x, -f*H.y, -f*H.z};
    for (int a3 = 0; a3 < 3; a3++) { double t = 1e-6; q_t X = {cos(t/sqrt(2.0)),0,0,0}; double sn = sin(t/sqrt(2.0));
      if (a3==0) X.x = sn; else if (a3==1) X.y = sn; else X.z = sn;
      q_t Xm = qc(X); double Fp, Fm, d;
      Fval(Pbuf, qm(qm(qm(qc(R), Vt), X), T).w, &Fp, &d); Fval(Pbuf, qm(qm(qm(qc(R), Vt), Xm), T).w, &Fm, &d);
      double fd = (log(Fp)-log(Fm))/(2*t); double err = fabs(fd-an[a3]);
      printf("  score dir %d analytic=%.8f fd=%.8f err=%.1e\n", a3, an[a3], fd, err); if (err > 1e-5*(1+fabs(fd))) ok = 0; }
    free(a);
  }
  /* small lattice smoke run: detailed-balance sanity -- run, refresh caches, compare cached logF */
  b = 4; L = 2; beta = 1.0; NMAX = 60; build(); baka = malloc(sizeof(double)*CAP*(NMAX+1));
  for (long e = 0; e < npath; e++) Vc[e] = Q1;
  refresh_all(); for (int i = 0; i < 20; i++) sweep();
  double maxd = 0; for (long e = 0; e < npath; e++) { double c0 = logF[e], c1 = pathlogF(e); maxd = fmax(maxd, fabs(c0-c1)); }
  /* also interiors cache */
  for (int I = 0; I < nint; I++) { q_t S = fullstaple(intlink[I]); maxd = fmax(maxd, qn(qadd(S, qs(Sint[I], -1)))); }
  printf("smoke: plaqP0=%.4f acc=%.3f cache maxdiff=%.2e badF=%d\n", plaqP0(), (double)nacc/ntry, maxd, badF);
  if (maxd > 1e-9 || badF) ok = 0;
  printf(ok ? "SELFTEST PASS\n" : "SELFTEST FAIL\n");
  return ok ? 0 : 1;
}

int main(int argc, char **argv) {
  if (argc >= 2 && !strcmp(argv[1], "selftest")) return selftest();
  if (argc < 10) { fprintf(stderr, "usage: pathcov b L beta Vmode ntherm nmeas every seed outprefix\n"); return 2; }
  b = atoi(argv[1]); L = atoi(argv[2]); beta = atof(argv[3]); const char *vm = argv[4];
  int ntherm = atoi(argv[5]), nmeas = atoi(argv[6]), every = atoi(argv[7]); rseed(strtoull(argv[8], 0, 10));
  const char *pre = argv[9];
  NMAX = 40 + (int)(8*beta*6); if (NMAX > 400) NMAX = 400;
  build(); baka = malloc(sizeof(double)*CAP*(NMAX+1));
  for (long e = 0; e < npath; e++) Vc[e] = !strcmp(vm, "haar") ? rhaar() : Q1;
  refresh_all();
  char fn[1024]; snprintf(fn, sizeof fn, "%s.bin", pre); FILE *fo = fopen(fn, "wb");
  snprintf(fn, sizeof fn, "%s.log", pre); FILE *fl = fopen(fn, "w");
  fprintf(fl, "{\"b\":%d,\"L\":%d,\"N\":%d,\"beta\":%.6f,\"V\":\"%s\",\"npath\":%ld,\"nint\":%d,\"NMAX\":%d,\"ntherm\":%d,\"nmeas\":%d,\"every\":%d}\n",
          b, L, N, beta, vm, npath, nint, NMAX, ntherm, nmeas, every); fflush(fl);
  float *buf = malloc(sizeof(float)*npath*3);
  for (int t = 0; t < ntherm + nmeas*every; t++) {
    sweep();
    if (t < ntherm && t % 20 == 19) { double r = (double)nacc/(ntry ? ntry : 1); eps *= r > 0.6 ? 1.15 : (r < 0.4 ? 0.87 : 1.0); nacc = ntry = 0; }
    if (t >= ntherm && (t - ntherm) % every == every - 1) {
      for (long e = 0; e < npath; e++) pathscore(e, buf + e*3);
      fwrite(buf, sizeof(float), npath*3, fo);
    }
    if (t % 50 == 0) { fprintf(fl, "sweep %d plaqP0 %.6f acc %.3f eps %.3f badF %d\n", t, plaqP0(), (double)nacc/(ntry?ntry:1), eps, badF); fflush(fl); }
  }
  fclose(fo); fprintf(fl, "done badF %d\n", badF); fclose(fl);
  return 0;
}
