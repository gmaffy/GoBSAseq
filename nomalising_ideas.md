# Normalising ideas: why CompositeZ is weak and what to do

## Current method

CompositeZ = (|Z_ΔSI| + Z_G′)/√2 (one-bulk: |Z_AFDev| + Z_G1). Each Z is a genome-wide robust Z of the smoothed track, (value − median)/(1.4826·MAD). Its thresholds are genome-wide: quantiles of the maximum CompositeZ over simulated no-QTL genomes.

## Evidence that it is weak

- **In calibration** (plant-level F2 simulation, 20-plant bulks), at the simulated QTL:
  - CompositeZ: 18.8 against a p99 threshold of 18.4 (1.02×);
  - HMM LOD: about 24 against 5 (about 5×);
  - BRM: found the QTL every time, at |z| ≥ about 2.9.
- **On the simulated VCF** the same kind of QTL gave CompositeZ 395 against a threshold of 20. The scale is not stable across datasets.
- **Genome-wide thresholds** are about 15–25 Z units, far from what a "Z" suggests (about 3).

## Why

1. **The scale comes from the data, not a null model.** Robust Z divides by the genome's own median and MAD. Smoothed tracks are dominated by pool-frequency drift that is shared along whole chromosome stretches, so a genome holds only a few independent blocks. The MAD is a noisy, small yardstick for that spread, and null genomes routinely reach large maxima in MAD units. This inflates the genome-wide threshold, so a real QTL must clear a much higher bar.
2. **The scale depends on genome composition.** How much of the genome is null and how noisy it is change the median and MAD, so the same QTL scores very differently in different datasets. CompositeZ values cannot be compared across studies.
3. **The two components are nearly redundant.** G′ is essentially smoothed ΔSI² scaled by depth, and both come from the same reads. Stouffer's method gains power by pooling independent evidence; with two near-copies it gains nothing.
4. **Mismatched shapes.** |Z_ΔSI| is folded (half-normal); Z_G′ is right-skewed and can be negative. Their sum has an awkward null distribution, which the null simulation absorbs only by raising the threshold.

By contrast, BRM standardises the smoothed ΔSI by its model null SD at each marker (pool + smoothed read variance with overdispersion ρ). Its null z has variance about 1 (0.85–0.96), and its genome-wide threshold is about 2.9.

## Suggestions (best first)

### 1. Null-SD standardisation with a correlation-aware Stouffer (recommended)

Keep the goal (one consensus score for ΔSI and G′), but standardise each track by its **null** distribution at each marker instead of the genome's spread:
- **z_ΔSI** = smoothed ΔSI / σ_null (exactly BRM's z: pool variance + `ReadVar0 + ρ·ReadVar1`). One-bulk: z_AFDev the same way.
- **z_G′** = (G′ − E₀[G′]) / SD₀[G′] at the marker: the null mean and SD from the genome-wide null simulation, which already produces G′ at every marker (or an analytic approximation).
- **Combine** with Stouffer's method corrected for their correlation: **ConsensusZ** = (z_ΔSI* + z_G′)/√(2 + 2r), r being their null correlation (estimated from the null genomes). z_ΔSI* is |z_ΔSI|, or the signed value if direction is to be kept.
- **Threshold** from the null simulation, as now (generalized extreme-value fit to genome maxima).

Expected effect: the scale no longer depends on genome composition, values are comparable across datasets, and the threshold should drop to near BRM's (about 3), giving power comparable to BRM's.

### 2. Use BRM's z as the ΔSI–G′ consensus

G′ adds almost nothing beyond ΔSI², so |z_ΔSI| (BRM's z) already is the consensus in practice. Final QTLs would require the HMM LOD + BRM; CompositeZ would be dropped or kept for display.

### 3. Keep CompositeZ but drop it from the final consensus

Report it as now, but final QTLs require the HMM LOD + BRM only. This is the smallest change, but it leaves a weak, dataset-dependent score in the output.

## Naming

- **These are true z-scores:** null-SD standardisation is the textbook (value − null mean)/null SD, closer to it than the current robust Z.
- **Caveat, normality:** z_ΔSI is close to N(0, 1) under the null; z_G′ has the right mean and SD but is right-skewed. Neither should be turned into p-values with normal tables; thresholds are genome-wide maxima from simulation.
- **If robust Z is replaced,** keep the `Z_*` column names and document them as "standardised by the null SD at the marker".
- **If both are kept,** give the new ones a distinct prefix (e.g. `NZ_DeltaSI`, `NZ_Gstat`) to avoid confusion.
- **Rename** CompositeZ to `ConsensusZ`, since its meaning and scale change.

## Checks before adopting

- **Calibration:** null genomes over the ConsensusZ p99 at about 1/100, with and without overdispersion, and in `2b` without `--ref-parent` and one-bulk modes.
- **Power:** many simulated QTLs of varying effect, comparing ConsensusZ, the current CompositeZ, BRM and the HMM LOD; ConsensusZ should beat CompositeZ and roughly match BRM.
- **z_G′ standardisation:** its null mean and SD per marker are stable (enough null genomes), and its correlation r with z_ΔSI is stable across datasets.
- **Final consensus:** re-measure its power once ConsensusZ replaces CompositeZ.
