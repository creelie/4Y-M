"""Interval re-check (mpmath.iv) of the numbers in sfexpand.tex.
The text derives each bound by hand; this script only re-checks them."""
from mpmath import iv, mpf
iv.dps = 40
e = iv.e
def chk(name, ok):
    print(('OK   ' if ok else 'FAIL ') + name)
    assert ok
# Theorem sfglstep(v): log eps(16) < -9.2, eps(t)=96 sqrt(2 e t) e^{-t} (1-e^{-t})^{-96}
t = iv.mpf(16)
loge1 = iv.log(96) + iv.log(2*e*t)/2 - t
print('log eps_1(16) in', loge1)
chk('log eps_1(16) < -9.20278', loge1.b < mpf('-9.20278'))
corr = -96*iv.log(1-iv.exp(-t))
print('-96 log(1-e^-16) in', corr)
chk('-96 log(1-e^-16) < 1e-4', corr.b < mpf('1e-4'))
chk('log eps(16) < -9.2', (loge1+corr).b < mpf('-9.2'))
# hand bounds used in the text
L2 = iv.log(2); L3 = iv.log(3)
chk('7.5 L2 + L3 - 15.5 < -9.20278', (7.5*L2+L3-15.5).b < mpf('-9.20278'))
chk('e^16 > 10^6', iv.exp(16).a > 10**6)
chk('2.7^16 > 7.9e6', (iv.mpf('2.7')**16).a > mpf('7.9e6'))
chk('96 x/(1-x) < 1e-4 at x = 1/7.9e6', (96*(1/iv.mpf('7.9e6'))/(1-1/iv.mpf('7.9e6'))).b < mpf('1e-4'))
chk('7.5*0.69314719 + 1.09861228... upper bounds of thm:lfkp give < -9.20278', (7.5*iv.mpf('0.69314719')+iv.mpf('1.09861229')-15.5).b < mpf('-9.20278'))
# threshold where eps(t)=e^-9 lies below 16
t2 = iv.mpf('15.9')
chk('eps(15.9) < e^-9 (so 16 is not sharp)', (iv.log(96)+iv.log(2*e*t2)/2-t2-96*iv.log(1-iv.exp(-t2))).b < -9)
# Proposition sfglbackground: K_* = 1 + Sbar * S_H
a = iv.mpf(1)/500
x = iv.exp(-a)
c_alpha = iv.mpf(1)/52 - ((iv.exp(3*a)+iv.exp(-3*a))/2-1)*(168+42*iv.sqrt(2))
print('c_alpha in', c_alpha)
chk('c_alpha > 0.0151', c_alpha.a > mpf('0.0151'))
chk('1-x >= 0.999/500', (1-x).a >= mpf('0.999')/500)
inv = iv.mpf(500)/mpf('0.999')
chk('500/0.999 < 500.51', inv.b < mpf('500.51'))
r = iv.mpf('500.51')
chk('r^2 < 250511', (r**2).b < 250511)
chk('250511^2 < 6.2756e10 (so r^4 < 6.2756e10)', (iv.mpf(250511)**2).b < mpf('6.2756e10'))
br = 1 + 384*iv.mpf('6.2756e10') + 16*iv.mpf(250511)
print('1 + 384*6.2756e10 + 16*250511 in', br)
chk('1 + 384 r^4 + 16 r^2 < 2.4100e13', br.b < mpf('2.4100e13'))
SH = iv.mpf('66.3')*45*iv.mpf('2.4100e13')
print('66.3*45*2.4100e13 in', SH)
chk('S_H < 7.20e16', SH.b < mpf('7.20e16'))
Sbar = 168+42*iv.sqrt(2)
chk('Sbar < 227.40', Sbar.b < mpf('227.40'))
chk('Sbar = 48 sqrt2 p_1', abs(Sbar - 48*iv.sqrt(2)*(7/(2*iv.sqrt(2))+iv.mpf(7)/8)).b < mpf('1e-30'))
chk('66.3 > 1/0.0151', (1/iv.mpf('0.0151')).b < mpf('66.3'))
K = 1 + iv.mpf('227.40')*iv.mpf('7.20e16')
print('K_* bound', K)
chk('K_* < 1.64e19', K.b < mpf('1.64e19'))
# blocks at l-infinity distance d: (2d+1)^4-(2d-1)^4 = 64 d^3 + 16 d
for d in range(1, 20):
    assert (2*d+1)**4-(2*d-1)**4 == 64*d**3+16*d
chk('(2d+1)^4-(2d-1)^4 = 64d^3+16d for d<20', True)
print('all checks passed')
