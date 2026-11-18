# Independent BSA-seq analysis used for 4_DG_BULK (high) vs 8_LG_BULK (low),
# MOSCHATA.v1.gatk_gatk.all.vcf.gz. Not part of GoBSAseq.
#
# Input preparation (run before this script):
#
#   V=/media/largeshareolder/DATA/moschata/MERGED_VCFs/v1/gatk_gatk/MOSCHATA.v1.gatk_gatk.all.vcf.gz
#   # 1. subset to the two bulks, per chromosome in parallel; drop sites that are
#   #    hom-ref or missing in both bulks
#   for i in $(seq -w 0 20); do echo Cmo_Chr$i; done | xargs -P 10 -I{} sh -c \
#     "bcftools view -r {} -s 4_DG_BULK,8_LG_BULK --trim-alt-alleles $V 2>/dev/null \
#      | bcftools view -i 'GT[*]!=\"RR\" && GT[*]!=\"mis\"' -Oz -o bulks.{}.vcf.gz"
#   ls bulks.*.vcf.gz | sort > list.txt
#   bcftools concat -f list.txt -Oz -o bulks.vcf.gz && tabix -p vcf bulks.vcf.gz
#   # 2. SNP table
#   bcftools query -i 'TYPE="snp"' -f '%CHROM\t%POS\t%REF\t%ALT\t%QUAL\t%INFO/QD\t%INFO/FS\t%INFO/SOR\t%INFO/MQ\t%INFO/MQRankSum\t%INFO/ReadPosRankSum[\t%GT\t%AD\t%DP]\n' \
#     bulks.vcf.gz > bulks.snps.tsv
#   # 3. pickle it
#   python -c "import pandas as pd; d=pd.read_csv('bulks.snps.tsv',sep='\t',header=None,na_values='.',low_memory=False); \
#     d.columns=['chr','pos','ref','alt','qual','qd','fs','sor','mq','mqrs','rprs','gt1','ad1','dp1','gt2','ad2','dp2']; d.to_pickle('snps.pkl')"
#
# Then:  python indep_bsa.py snps.pkl indep_2Mb 2e6
# Needs: numpy, pandas, scipy.
#
"""Independent BSA-seq analysis (G' after Magwene 2011 + |dSNP| with Takagi CI).
Usage: python indep_bsa.py <snps.pkl> <outdir> [window_bp]
"""
import sys, os
import numpy as np, pandas as pd
from scipy import stats

pkl, out = sys.argv[1], sys.argv[2]
WIN = float(sys.argv[3]) if len(sys.argv) > 3 else 2e6
NH, NL = 20, 16          # plants per bulk (DG = high, LG = low)
MIN_DP, MAX_DP = 20, 200
os.makedirs(out, exist_ok=True)
rng = np.random.default_rng(1)

d = pd.read_pickle(pkl)
n0 = len(d)
log = []
def note(msg):
    print(msg); log.append(msg)

note(f"SNP records in bulk subset: {n0:,}")
d = d[~d.alt.str.contains(',', na=True) & (d.ref.str.len() == 1) & (d.alt.str.len() == 1)]
note(f"biallelic SNVs: {len(d):,}")
q = (d.qual >= 30) & (d.qd.fillna(99) >= 2) & (d.fs.fillna(0) <= 60) & \
    (d.sor.fillna(0) <= 3) & (d.mq.fillna(99) >= 40)
d = d[q]
note(f"after GATK hard filters (QUAL>=30, QD>=2, FS<=60, SOR<=3, MQ>=40): {len(d):,}")

def split_ad(s):
    a = s.str.split(',', expand=True)
    return pd.to_numeric(a[0], errors='coerce'), pd.to_numeric(a[1], errors='coerce')
d['r1'], d['a1'] = split_ad(d.ad1)
d['r2'], d['a2'] = split_ad(d.ad2)
d = d.dropna(subset=['r1', 'a1', 'r2', 'a2'])
for c in ['r1', 'a1', 'r2', 'a2']:
    d[c] = d[c].astype(int)
d['n1'] = d.r1 + d.a1
d['n2'] = d.r2 + d.a2
d = d[(d.n1 >= MIN_DP) & (d.n2 >= MIN_DP) & (d.n1 <= MAX_DP) & (d.n2 <= MAX_DP)]
note(f"bulk AD depth {MIN_DP}-{MAX_DP} in both: {len(d):,}")
# Segregating filter: pooled alt freq across both bulks 0.2-0.8 (an F2 marker
# is ~0.5 in the union of bulks; fixed/near-absent sites are non-informative).
paf = (d.a1 + d.a2) / (d.n1 + d.n2)
d = d[(paf >= 0.2) & (paf <= 0.8)].copy()
note(f"segregating (pooled ALT freq 0.2-0.8): {len(d):,}")

d['si1'] = d.a1 / d.n1
d['si2'] = d.a2 / d.n2
d['dsi'] = d.si1 - d.si2

def gstat(a1, r1, a2, r2):
    obs = np.stack([a1, r1, a2, r2]).astype(float)
    n1, n2 = a1 + r1, a2 + r2
    A, R, N = a1 + a2, r1 + r2, n1 + n2
    exp = np.stack([n1 * A / N, n1 * R / N, n2 * A / N, n2 * R / N])
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.where(obs > 0, obs * np.log(obs / exp), 0.0)
    return 2 * t.sum(0)
d['G'] = gstat(d.a1.values, d.r1.values, d.a2.values, d.r2.values)

def tricube_smooth(pos, vals, win):
    """Nadaraya-Watson tricube smoother, half-width win/2 (Magwene)."""
    h = win / 2
    out = np.empty((len(vals), len(pos)))
    lo = np.searchsorted(pos, pos - h)
    hi = np.searchsorted(pos, pos + h, side='right')
    for i in range(len(pos)):
        x = np.abs(pos[lo[i]:hi[i]] - pos[i]) / h
        w = (1 - x ** 3) ** 3
        sw = w.sum()
        for k, v in enumerate(vals):
            out[k, i] = (w * v[lo[i]:hi[i]]).sum() / sw
    return out

# ---- null simulation for |dSNP| (Takagi-style): F2 pool freq ~ Binom(2n,0.5)/2n,
# reads ~ Binom(depth, p). Smoothing reduces variance by the effective #markers;
# simulate per-marker raw |dsi| null, then CI on window means using the local
# marker count & depths (bootstrap of the window).
def sim_raw_dsi(n1, n2, reps):
    p1 = rng.binomial(2 * NH, 0.5, size=(reps, len(n1))) / (2 * NH)
    p2 = rng.binomial(2 * NL, 0.5, size=(reps, len(n2))) / (2 * NL)
    return rng.binomial(n1, p1) / n1 - rng.binomial(n2, p2) / n2

res = []
autos = [c for c in sorted(d.chr.unique()) if c != 'Cmo_Chr00']
for c in sorted(d.chr.unique()):
    s = d[d.chr == c].sort_values('pos').copy()
    pos = s.pos.values.astype(float)
    # Hampel/outlier-robust G': smooth G; also smooth |dsi|, dsi and depth-wt.
    sm = tricube_smooth(pos, [s.G.values, np.abs(s.dsi.values), s.dsi.values,
                              s.si1.values, s.si2.values], WIN)
    s['Gprime'], s['absdsi_sm'], s['dsi_sm'], s['si1_sm'], s['si2_sm'] = sm
    # markers per window (for CI)
    h = WIN / 2
    s['nSNP'] = np.searchsorted(pos, pos + h, side='right') - np.searchsorted(pos, pos - h)
    res.append(s)
r = pd.concat(res)

# ---- G' null: Magwene/QTLseqr approach. Remove outlier regions with Hampel on
# log G', fit lognormal to the rest (autosomes only), get p, BH-FDR.
A = r[r.chr != 'Cmo_Chr00']
lg = np.log(A.Gprime.values)
med = np.median(lg)
mad = stats.median_abs_deviation(lg, scale='normal')
keep = np.abs(lg - med) <= 5.2 * mad
mu = np.mean(lg[keep]); sd = np.std(lg[keep])
med_fit = np.median(lg[keep])
mode_fit = med_fit  # QTLseqr uses median as location
r['pval'] = stats.lognorm.sf(r.Gprime.values, s=sd, scale=np.exp(med_fit))
def bh(p):
    p = np.asarray(p); n = len(p); o = np.argsort(p)
    q = p[o] * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n); out[o] = np.minimum(q, 1); return out
r['qval'] = np.nan
mask = r.chr != 'Cmo_Chr00'
r.loc[mask, 'qval'] = bh(r.loc[mask, 'pval'].values)
note(f"G' null (lognormal on log G', Hampel-trimmed): median={np.exp(med_fit):.3f}, sdlog={sd:.3f}, kept {keep.mean():.1%}")

# ---- |dSNP| genome-wide 99% threshold from simulation: simulate raw null dsi
# with actual depths, smooth with the same kernel on a few chromosomes'
# real positions; take genome-wide max per replicate -> 95/99% of maxima.
def null_maxima(reps=200):
    mx = np.zeros(reps)
    mx_point = []
    for c in autos:
        s = r[r.chr == c]
        pos = s.pos.values.astype(float)
        # thin to ~ every 20th marker for speed; smoothing variance preserved
        # approximately because weights scale. Use real full set where small.
        step = max(1, len(pos) // 4000)
        idx = np.arange(0, len(pos), step)
        p_, n1, n2 = pos[idx], s.n1.values[idx], s.n2.values[idx]
        # simulate genuine linkage: pool freq is a step function along chrom.
        # simulate F2 haplotypes with Haldane, 1 Morgan per chrom.
        L = max(p_.max(), 1)
        for rep in range(reps):
            def pool(n):
                cnt = np.zeros(len(p_))
                for _ in range(2 * n):
                    k = rng.poisson(1.0)
                    br = np.sort(rng.uniform(0, L, k))
                    start = rng.integers(2)
                    seg = np.searchsorted(br, p_)  # segment index
                    cnt += (start + seg) % 2
                return cnt / (2 * n)
            f1, f2 = pool(NH), pool(NL)
            flip = rng.integers(2, size=len(p_)).astype(bool)  # unknown polarity
            f1 = np.where(flip, 1 - f1, f1); f2 = np.where(flip, 1 - f2, f2)
            dsi = rng.binomial(n1, f1) / n1 - rng.binomial(n2, f2) / n2
            g = gstat(rng.binomial(n1, f1), 0, 0, 0) if False else None
            sm = tricube_smooth(p_, [np.abs(dsi)], WIN)[0]
            mx[rep] = max(mx[rep], sm.max())
    return mx

note("simulating |dSNP| null (200 genomes, linked F2 haplotypes, Haldane 1 M/chrom)...")
mx = null_maxima(200)
thr95, thr99 = np.quantile(mx, [0.95, 0.99])
note(f"|dSNP| genome-wide null maxima: 95%={thr95:.4f}, 99%={thr99:.4f}")
r['absdsi_thr99'] = thr99

r.to_csv(f"{out}/indep_markers.tsv.gz", sep='\t', index=False,
         columns=['chr', 'pos', 'ref', 'alt', 'r1', 'a1', 'r2', 'a2', 'si1', 'si2', 'dsi', 'G',
                  'Gprime', 'pval', 'qval', 'absdsi_sm', 'dsi_sm', 'si1_sm', 'si2_sm', 'nSNP'])

# ---- call regions: G' q<0.01, contiguous; also |dSNP| > 99% thr
def regions(df, flag):
    out = []
    for c, s in df.groupby('chr'):
        s = s.sort_values('pos'); f = flag.loc[s.index].values
        if not f.any(): continue
        edges = np.diff(np.r_[0, f.astype(int), 0])
        st, en = np.where(edges == 1)[0], np.where(edges == -1)[0]
        for a, b in zip(st, en):
            seg = s.iloc[a:b]
            pk = seg.Gprime.idxmax()
            out.append(dict(chr=c, start=int(seg.pos.min()), end=int(seg.pos.max()),
                            nSNP=len(seg), peak_pos=int(s.loc[pk, 'pos']),
                            max_Gprime=round(seg.Gprime.max(), 2), min_q=float(seg.qval.min()),
                            max_absdSNP=round(seg.absdsi_sm.max(), 4),
                            SI_DG_at_peak=round(s.loc[pk, 'si1_sm'], 3),
                            SI_LG_at_peak=round(s.loc[pk, 'si2_sm'], 3)))
    return pd.DataFrame(out)

Ar = r[r.chr != 'Cmo_Chr00']
for name, flag in [('Gprime_q01', Ar.qval < 0.01), ('Gprime_q05', Ar.qval < 0.05),
                   ('absdSNP_99', Ar.absdsi_sm > thr99), ('absdSNP_95', Ar.absdsi_sm > thr95)]:
    reg = regions(Ar, flag)
    reg.to_csv(f"{out}/regions_{name}.tsv", sep='\t', index=False)
    note(f"\n== {name}: {len(reg)} regions")
    if len(reg): note(reg.to_string(index=False))

# per-chrom summary
summ = r.groupby('chr').agg(nSNP=('pos', 'size'), maxGp=('Gprime', 'max'),
                            maxAbsD=('absdsi_sm', 'max'), minQ=('qval', 'min'))
note("\n" + summ.round(4).to_string())
open(f"{out}/log.txt", 'w').write("\n".join(log))
