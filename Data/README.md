# Simulated test data

Small simulated two-parent, two-bulk (2p2b) F2 datasets with a known truth, for
testing GoBSAseq end to end. They are produced by `simulate_bsa.py`, which is
independent of GoBSAseq's own null-genome simulation (`stats/nullsim.go`).

| File | QTL | Expected result |
|------|-----|-----------------|
| `sim_F2_strongQTL.vcf.gz` | chr2:15,000,000, strong (complete selection) | QTL called on chr2 by CompositeZ, individual stats and BRM; nothing elsewhere |
| `sim_F2_weakQTL.vcf.gz` | chr2:15,000,000, weak | No calls: with 50-plant bulks this QTL is weaker than chance drift at a 1% genome-wide error rate |
| `sim_F2_null.vcf.gz` | none | No calls (checks false positives) |

## Design (all three files)

- 5 chromosomes × 30 Mb, each 1 Morgan (100 cM) long, Haldane crossovers
- One SNP every ~10 kb (≈15,000 markers)
- Samples: `P1` (high parent), `P2` (low parent), `HB` (high bulk), `LB` (low bulk)
- Bulks: 50 F2 plants each; bulk depth ~ Normal(60, 15)
- Parents: homozygous, depth 20
- `FORMAT=GT:AD:DP:GQ`, `QUAL=100`, `FILTER=PASS`, no INFO annotations
  (so the filter is detected as `deepvariant-like`)

Weak-QTL selection: plants with 0/1/2 low-parent alleles at the QTL enter the
high bulk with probability 1.0/0.55/0.3 (the low bulk is the mirror image).

## Run

```bash
gobsaseq -V Data/sim_F2_strongQTL.vcf.gz -P P1,P2 -B HB,LB -S 50,50 -m F2 -o results_strong
```

Results with the default `--seed 1` and `--rep 1000`:

| Dataset | CompositeZ P99 | CompositeZ calls | Individual-stat calls | BRM blocks |
|---------|---------------:|------------------|-----------------------|------------|
| strongQTL | 28.4 | 1 (chr2) | 1 (chr2) | 1 (chr2, peak 15.44 Mb) |
| weakQTL | 29.2 | 0 | 0 | 0 |
| null | 31.2 | 0 | 0 | 0 |

## Regenerate or make new ones

```bash
python3 simulate_bsa.py --strength strong --seed 1 | bgzip -c > sim_F2_strongQTL.vcf.gz
python3 simulate_bsa.py --strength weak   --seed 7 | bgzip -c > sim_F2_weakQTL.vcf.gz
python3 simulate_bsa.py --strength none   --seed 1 | bgzip -c > sim_F2_null.vcf.gz
python3 simulate_bsa.py --help   # bulk size, depth, marker spacing, QTL position, ...
```

Many null datasets with different `--seed` values can be used to check
calibration: with P99 thresholds, about 1% of them should produce any call.
