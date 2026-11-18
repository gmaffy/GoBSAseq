#!/usr/bin/env python3
"""Simulate a two-parent, two-bulk (2p2b) BSA-seq VCF with a known QTL.

Bulk individuals are simulated as F2 plants: each is two gametes from an F1,
with Haldane (Poisson) crossovers along every chromosome. Plants are admitted
to a bulk with a probability that depends on their genotype at the QTL, so the
truth (QTL position and strength) is known exactly. Reads at each marker are
drawn from the pool allele frequency with depth ~ Normal(mean, sd).

This simulator is independent of GoBSAseq's own null-genome code
(stats/nullsim.go) and is meant for testing it, e.g.:

    python3 simulate_bsa.py --strength strong --seed 1 > sim_F2_strongQTL.vcf

Strengths (probability a plant with 0/1/2 low-parent alleles at the QTL
enters the HIGH bulk; the LOW bulk mirrors it):
    strong:   1.0 / 0.0  / 0.0   (complete selection on homozygotes)
    moderate: 1.0 / 0.35 / 0.1
    weak:     1.0 / 0.55 / 0.3
    none:   no QTL (null genome, for calibration checks)
"""
import argparse
import random

STRENGTHS = {
    "strong": [1.0, 0.0, 0.0],
    "moderate": [1.0, 0.35, 0.1],
    "weak": [1.0, 0.55, 0.3],
    "none": None,
}


def gamete(rng, length, morgans):
    """Return (start allele, crossover positions) of a gamete from an F1."""
    xs = []
    t = rng.expovariate(morgans)
    while t < 1:
        xs.append(t * length)
        t += rng.expovariate(morgans)
    return rng.randint(0, 1), xs


def allele_at(g, pos):
    start, xs = g
    return start ^ (sum(1 for x in xs if x < pos) % 2)


def make_bulk(rng, n, chrom, args, high):
    probs = STRENGTHS[args.strength]
    inds = []
    while len(inds) < n:
        g1 = gamete(rng, args.chrom_len, args.morgans)
        g2 = gamete(rng, args.chrom_len, args.morgans)
        if probs is not None and chrom == args.qtl_chrom:
            low_alleles = allele_at(g1, args.qtl_pos) + allele_at(g2, args.qtl_pos)
            p = probs[low_alleles] if high else probs[2 - low_alleles]
            if rng.random() > p:
                continue
        inds.append((g1, g2))
    return inds


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--strength", choices=STRENGTHS, default="strong")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--chroms", type=int, default=5)
    ap.add_argument("--chrom-len", type=int, default=30_000_000)
    ap.add_argument("--morgans", type=float, default=1.0, help="genetic length of each chromosome")
    ap.add_argument("--spacing", type=int, default=10_000, help="mean marker spacing (bp)")
    ap.add_argument("--bulk-size", type=int, default=50)
    ap.add_argument("--depth", type=float, default=60, help="mean bulk depth")
    ap.add_argument("--depth-sd", type=float, default=15)
    ap.add_argument("--qtl-chrom", default="chr2")
    ap.add_argument("--qtl-pos", type=int, default=15_000_000)
    args = ap.parse_args()

    rng = random.Random(args.seed)
    chroms = [f"chr{i + 1}" for i in range(args.chroms)]

    print("##fileformat=VCFv4.2")
    print(f"##simulation=simulate_bsa.py strength={args.strength} seed={args.seed} "
          f"bulk_size={args.bulk_size} qtl={args.qtl_chrom}:{args.qtl_pos} morgans={args.morgans}")
    for c in chroms:
        print(f"##contig=<ID={c},length={args.chrom_len}>")
    print('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    print('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">')
    print('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">')
    print('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">')
    print("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1\tP2\tHB\tLB")

    def sample(ref_n, alt_n):
        gt = "0/0" if alt_n == 0 else ("1/1" if ref_n == 0 else "0/1")
        return f"{gt}:{ref_n},{alt_n}:{ref_n + alt_n}:50"

    for c in chroms:
        hb = make_bulk(rng, args.bulk_size, c, args, high=True)
        lb = make_bulk(rng, args.bulk_size, c, args, high=False)
        pos = 1000
        while pos < args.chrom_len:
            ref, alt = rng.sample("ACGT", 2)
            p1_is_ref = rng.random() < 0.5  # P1 = high parent

            def bulk_ad(inds):
                low = sum(allele_at(g1, pos) + allele_at(g2, pos) for g1, g2 in inds)
                f_low = low / (2 * len(inds))
                f_alt = f_low if p1_is_ref else 1 - f_low
                d = max(5, int(rng.gauss(args.depth, args.depth_sd)))
                a = sum(1 for _ in range(d) if rng.random() < f_alt)
                return d - a, a

            p1 = sample(20, 0) if p1_is_ref else sample(0, 20)
            p2 = sample(0, 20) if p1_is_ref else sample(20, 0)
            hr, ha = bulk_ad(hb)
            lr, la = bulk_ad(lb)
            print(f"{c}\t{pos}\t.\t{ref}\t{alt}\t100\tPASS\t.\tGT:AD:DP:GQ\t{p1}\t{p2}\t"
                  f"{sample(hr, ha)}\t{sample(lr, la)}")
            pos += rng.randint(args.spacing // 2, args.spacing * 3 // 2)


if __name__ == "__main__":
    main()
