"""Figure of the plaquette blocks: the block plaquettes in the planes {0,1}
and {2,3}, with the links of the line field. Run from this directory:
python3 make_plaqblocks.py. Writes fig_plaqblocks.png (flattened, 256 colours)."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from make_figures import save

plt.rcParams.update({'font.family': 'serif', 'font.size': 10, 'mathtext.fontset': 'cm'})
SHADE, GREY, RED, INK = '#d6e4f5', '#9a9a9a', '#c0392b', '#222222'
ROWS = {'0000': 0, '1000': 0, '0100': 2, '1100': 3, '0010': 3, '1010': 1, '0110': 2, '1110': 1,
        '0001': 1, '1001': 2, '0101': 1, '1101': 3, '0011': 3, '1011': 2, '0111': 0, '1111': 0}


def sigma(x):
    return ROWS[''.join(str(c % 2) for c in x)]


def panel(ax, mu, nu, fixed):
    n = 4
    for a in range(n):
        for b in range(n):
            if a % 2 == b % 2:
                ax.add_patch(Rectangle((a, b), 1, 1, color=SHADE, lw=0))
    for a in range(n + 1):
        for b in range(n + 1):
            x = [0, 0, 0, 0]; x[mu], x[nu] = a, b
            for d, (da, db) in ((mu, (1, 0)), (nu, (0, 1))):
                if a + da > n or b + db > n:
                    continue
                y = list(x)
                for k, v in fixed.items():
                    y[k] = v
                red = sigma(y) == d
                ax.plot([a, a + da], [b, b + db], color=RED if red else GREY,
                        lw=2.6 if red else 0.9, solid_capstyle='round', zorder=3 if red else 2)
            y = list(x)
            for k, v in fixed.items():
                y[k] = v
            out = sigma(y) not in (mu, nu)
            ax.scatter(a, b, s=22, facecolor='white' if out else INK, edgecolor=INK, lw=0.8, zorder=4)
    ax.set_aspect('equal'); ax.set_xlim(-0.4, n + 0.4); ax.set_ylim(-0.6, n + 0.4)
    ax.set_xticks(range(n + 1)); ax.set_yticks(range(n + 1))
    ax.set_xlabel(r'$x_%d$' % mu); ax.set_ylabel(r'$x_%d$' % nu)
    for sp in ax.spines.values():
        sp.set_visible(False)
    ax.tick_params(length=0)


fig, (a1, a2) = plt.subplots(1, 2, figsize=(7.2, 3.4))
panel(a1, 0, 1, {2: 0, 3: 0}); a1.set_title(r'(a) plane $\{0,1\}$, $x_2=x_3=0$', fontsize=10)
panel(a2, 2, 3, {0: 0, 1: 0}); a2.set_title(r'(b) plane $\{2,3\}$, $x_0=x_1=0$', fontsize=10)
fig.subplots_adjust(left=0.06, right=0.99, bottom=0.14, top=0.9, wspace=0.25)
save(fig, 'fig_plaqblocks.png')
