"""Interval re-check (mpmath iv) of the SU(3) numbers of Theorem groups13 as derived in the text.
The text does not depend on this script."""
from mpmath import iv, mpf
iv.dps = 50
def I(x): return iv.mpf(x)
ln = iv.log
ok = True
def check(cond, msg):
    global ok
    print(('OK   ' if cond else 'FAIL ')+msg); ok &= bool(cond)
# C_0(SU(3)), via the identity in the text
C0a = I(1)/2+I(2)/3*(ln(I(256)/35)+4*ln(I(224)/9)-7*ln(I(95)/96))
C0b = I(1)/2+I(2)/3*(I('29.5')*ln(2)-8*ln(3)+2*ln(5)-I('1.5')*ln(I(50)/49)+7*ln(I(96)/95))
check(abs(C0a.a-C0b.a) < mpf('1e-40'), 'identity for C_0: 63 log2 - log3 - 8 log5 + 3 log7 - 7 log19 = 29.5 log2 - 8 log3 + 2 log5 - 1.5 log(50/49) + 7 log(96/95)')
check(C0a.a >= mpf('10.44720') and C0a.b <= mpf('10.44721'), 'C_0 in [10.44720,10.44721]: %s' % C0a)
# I_8(pi) = int_0^pi sin^7 = 32/35
check(I(2)*6*4*2/(7*5*3) == I(32)/35, 'I_8(pi)=32/35')
C0 = iv.mpf(['10.44720','10.44721'])
H=6**4; M=11**4
def Lq(D): return H*(ln(2)+2*ln(I(D)))+1+4*ln(11)+4*ln(6)
L1=Lq(2400); L5=Lq(23**4-1)
check(abs(L1-I('21090.194')).b<=mpf('1e-3') and abs(L5-I('33424.872')).b<=mpf('1e-3'), '|log q*| within 1e-3 of 21090.194, 33424.872')
def g(beta,rho,L): return beta*rho**2/24-16*ln(beta)-6*C0-L
def bracket(lo,hi,rho,L,name,margin):
    a=g(I(lo),I(rho),L); b=g(I(hi),I(rho),L)
    check(a.b < -margin and b.a > margin, '%s: g(lo) <= %s, g(hi) >= %s' % (name, a.b, b.a))
    check(I(lo) > 384/I(rho)**2, '%s: interval right of the maximum 384/rho^2' % name)
Lw1 = iv.mpf(['21090.193','21090.195']); Lw5 = iv.mpf(['33424.871','33424.873'])
bracket('5.268e21','5.270e21','1e-8',Lw1,'beta*(1e-8,1)',1)
bracket('8.230e21','8.232e21','1e-8',Lw5,'beta*(1e-8,5)',1)
bracket('5.285e23','5.287e23','1e-9',Lw1,'beta*(1e-9,1)',1)
bracket('8.247e23','8.249e23','1e-9',Lw5,'beta*(1e-9,5)',1)
bracket('2.971e22','2.973e22','1e-8',Lw5+H*30*ln(10),'beta*(1e-8,5,1e30)',1)
def eps(beta,L): return (16*ln(I(beta))+6*C0)/L
e_lo=min(eps('5.268e21',Lw1).a,eps('8.230e21',Lw5).a,eps('5.285e23',Lw1).a,eps('8.247e23',Lw5).a)
e_hi=max(eps('5.270e21',Lw1).b,eps('8.232e21',Lw5).b,eps('5.287e23',Lw1).b,eps('8.249e23',Lw5).b)
check(e_lo>mpf('0.026') and e_hi<mpf('0.045'), 'eps_G in (%s, %s)' % (e_lo,e_hi))
# scales j (rho=1e-8, s=5)
c=I('1.01')*I('1.023')
check(16*c**2 <= I('17.0810278') and 24*c**2 <= I('25.6215417'), '16c^2<=17.0810278, 24c^2<=25.6215417')
def Lj(j):
    Hj=8**4 if j<=2 else 9**4
    return Hj*(ln(2)+2*ln(I(23**4-1)))+1+4*ln(I(17*2**j-4))+ln(I(Hj))
check(abs(Lj(1)-I('105605.899')).b<=mpf('1e-3') and abs(Lj(2)-I('105608.929')).b<=mpf('1e-3'), '|log q*_1|, |log q*_2|')
A=6561*(ln(2)+2*ln(I(23**4-1)))+1+8*ln(3)
check(abs(A-I('169133.301')).b<=mpf('1e-3') and (A+4*ln(17)).b<=mpf('169144.64'), 'A and A+4log17<=169144.64')
def gj(beta,j):
    rho0=I('1e-8')/(c*4**j)
    return beta*rho0**2/24-16*ln(beta)-6*C0-Lj(j)
for j,lo,hi in [(1,'4.367e23','4.368e23'),(2,'6.990e24','6.992e24'),(3,'1.785e26','1.787e26')]:
    a=gj(I(lo),j); b=gj(I(hi),j)
    check(a.b< -1 and b.a> 1, 'beta*_%d in [%s,%s]: margins %s %s' % (j,lo,hi,a.b,b.a))
b=I('4.38e22')
t1=b*I('1e-16')/I('25.6215417'); t2=16*ln(b); t3=6*C0; t4=I('169144.64')+68*12*ln(2)
check(t1.a>=mpf('170949.9') and t2.b<=mpf('834.15') and t3.b<=mpf('62.684') and t4.b<=mpf('169710.3'), 'terms at 4.38e22: %s %s %s %s' % (t1.a,t2.b,t3.b,t4.b))
for j in range(1,13):
    v=gj(16**j*b,j); check(v.a>0, '16^-j beta*_j <= 4.38e22 at j=%d (g=%s)' % (j, v.a))
check(abs(C0a-I('10.4472096')).b<=mpf('1e-7'), 'C_0 within 1e-7 of 10.4472096')
check((Lj(1)+64*ln(2)).b<=mpf('105697.7') and (Lj(2)+128*ln(2)).b<=mpf('105697.7'), 'bracket for j<=2 at most 105697.7')
check(all((Lj(j)).b<=(I('169144.64')+4*j*ln(2)).a for j in range(3,13)), '|log q*_j| <= 169144.64+4j log2 for 3<=j<=12')
from fractions import Fraction as Fr
check(Fr(2)**63/3/Fr(5)**8*Fr(7)**3/Fr(19)**7 == Fr(256,35)*Fr(224,9)**4*Fr(96,95)**7, 'factorisation 2^63 3^-1 5^-8 7^3 19^-7')
# (iii), (iv)
check(I('25980.8')**2 >= I('1.5e4')**2*3, 'c_3 <= 25980.8')
c6=I('6.81e5')*8*I('6.75e8'); check(c6.b<=mpf('3.678e15') and c6.a==mpf('3.6774e15'), 'c_6 = 3.6774e15')
check((1440*c6).b<=mpf('5.30e18'), '(b): 480*3*c_6 = %s' % (1440*c6).b)
check((1/(2758*I('7.8e7')*I('25980.8'))).a>=mpf('1.789e-16'), '(a): rho <= 1.789e-16')
check(I('198.0002')/3<=I('66.0001') and I('4.65e6')*8/3==I('1.24e7'), 'Lambda_2 constants')
# joining distance 2s+8
LN1=iv.mpf(['32481.66','32481.68']); LN5=iv.mpf(['38354.06','38354.08'])
check(abs(Lq((4*1+17)**4-1)-I('32481.67')).b<=mpf('1e-2') and abs(Lq((4*5+17)**4-1)-I('38354.07')).b<=mpf('1e-2'), 'joining |log q*|')
bracket('2.569e37','2.571e37','1.78e-16',LN1,'join s=1',1)
bracket('3.014e37','3.016e37','1.78e-16',LN5,'join s=5',1)
bs=I('2.571e37'); r=43*I('1.78e-16')
c2=2880*bs*r**2*(I('66.0001')+I('1.24e7')/bs)**2*45*14**2*9**2*iv.exp(I('0.008'))
check(iv.log10(c2).b<=mpf('22.14'), 'c_2 <= 10^%s' % iv.log10(c2).b)
check((5.30e18 <= 2.569e37), '(b) holds at the joining thresholds')
print('ALL OK' if ok else 'SOME FAILED')
