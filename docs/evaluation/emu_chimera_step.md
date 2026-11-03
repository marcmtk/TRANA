# Should TRANA remove chimeras before emu?

Evaluation on the 2026-07-02 MAB114 run (`trana_full`, 1200–1800 bp, filtlong mean quality ≥94, emu 3.6.2 with the 2026 database, exclusion list on).
Everything below was produced by the scripts in `/home/marc/Projects/2026-07-02/chimera_step_eval/scripts/` (local, not in this repository); tables are the `results_*.tsv` files there. No pipeline output was edited. The recommended check is now in the pipeline as `--check_emu_calls` (off by default; `bin/check_emu_calls.py`, module `EMU_CHECK_CALLS`).

## Recommendation

**Don't add a read-level chimera-removal step before emu.** Instead, add a **post-emu evidence check** that labels a species call *unconfirmed* when none of its reads matches that species' references almost exactly. Default: on, as a label only (nothing is removed), once it has passed at least one more run with known composition; until then, keep it off (see "Still uncertain").

- Rule (prototype `scripts/prototype_near_exact.py`): count the reads emu assigns to the species with P ≥ 0.5. A read is *near-exact* if its best alignment to a reference of that species scores ≥ 99% of the read length (score = matches − mismatches − inserted − deleted bases; about 7 mismatches or 14 indel bases in a 1,450 bp read). A call is confirmed with ≥ 2 near-exact reads.
- Where: after `EMU_ABUNDANCE`, reading emu's own SAM, so emu needs `--keep-files`. It needs no savont, no extra database and no extra container.
- On this run it labels 11 of the 13 unexpected calls: all except AB11 *S. haemolyticus* (real sequence) and AB05 *E. fergusonii* (an *E. coli* tie). It also labels the chimeric FDAARGOS_522 synonym call in AB12. It keeps every expected species, including *S. pneumoniae* in AB12 (15 near-exact reads). For thinned real rare species it detects within 8 points of emu at 10–12 reads and within 2 points from 15 reads (Q3).
- As a pre-emu step, `vsearch --uchime_ref` against the emu database is the best option tested. It works much better on these reads than expected, but it does the job less well than the post-emu rule:
  - it leaves 3 unexpected calls in AB11;
  - it creates a new false *E. fergusonii* call in AB14 by shifting the EM;
  - it raises the number of reads a chimera-prone rare species needs to be called from about 10–12 to about 20.

  If Marc still wants a pre-emu step, this is the one, but default off.

Calls a chimera method can't explain stay: *S. haemolyticus* in AB11 (real sequence, 8 near-exact reads) and *E. fergusonii* in AB05 (an *E. coli* tie, 22 near-exact reads). Those need NTC/blank comparison and tie reporting, not a chimera filter.

## Key facts established

### The reads are much better than assumed
Median read Q in this run is 23–25 per sample, not Q12–20. Only 2–3% of reads are below Q16 (`scripts/quality_dependence.py`). The median residual error of a read against its best reference is 0.75–1.1% in the mixes (`scripts/explore1.py`). That is why reference-based chimera detection works on single reads here.

### Composition (from Marc's message in barbell thread 6b931188, position 345)
AB01 *K. pneumoniae*, AB02 *M. morganii*, AB03 *K. aerogenes*, AB04 *E. cloacae/hormaechei*, AB05 *E. coli*, AB06 *S. epidermidis*, **AB07 *E. faecalis* + *S. gallolyticus* (a two-species mix, not a pure culture)**, AB08 *S. agalactiae*, AB09 *S. aureus*, AB10 *S. pneumoniae*.
AB11: AB01–05. AB12: AB06–10 (including both AB07 species). AB13: *K. pneumoniae*, *K. aerogenes*, *E. faecalis*, *S. gallolyticus*. AB14: *E. coli*, *S. aureus*, *S. epidermidis*. AB15: PBS NTC.

Name mapping used (`scripts/truth.py`):
- *S. pasteurianus* = *S. gallolyticus* subsp. *pasteurianus*.
- *E. hormaechei* and *E. cloacae* complex sp. AR_0002 = the stated *E. cloacae/hormaechei*.
- *Streptococcus* sp. FDAARGOS_520/521/522 = *S. agalactiae*: emu gives 19% of pure AB08's reads to FDAARGOS_522.

What isn't known: the cultures' purity was not checked independently, and there is no species-level truth for contamination such as *S. haemolyticus* in AB11.

### Method used to test filters cheaply
Emu spends almost all of its time in minimap2: AB12's emu log shows 2,296 s of minimap2 in a 39 min task. I ran emu's exact minimap2 command once per sample and stored the alignments (`scripts/align_all.sh`, `sam2table.py`). Emu's scoring and EM are re-implemented in `scripts/emlib.py`. It reproduces the published `rel-abundance` tables for AB05, AB11–14 and the NTC: identical species, counts within 1e-13 (`scripts/validate_em.py`). That makes each filter, spike or thinning a 1–3 s EM rerun.

**Check against real emu on filtered input.** I removed the 2,376 reads uchime_ref flags in AB12 and ran the real emu 3.6.2 container with TRANA's arguments plus `--keep-files` (`scripts/real_emu.sh`). The result matched `emlib`'s prediction exactly: the same 7 species, a largest count difference of 0.0 and 49 EM iterations in both (`scripts/validate_real_filtered.py`). On that run's real SAM, the prototype `prototype_near_exact.py` gave the same near-exact counts per species as the analysis code.

## Q1. Do chimera methods work on these reads?

Test material:
- **Real-read chimeras (`scripts/make_sim.py`):** 6,000 simulated chimeras built from the left part of a real pure-culture read and the right part of another, joined at homologous 16S positions. The noise is real and the truth is known.
- **Non-chimeric reads:** 27,000 pure-culture reads (AB01–06 and AB08–10, 3,000 each). AB07 is a two-species mix and is kept apart.
- **Real mixes:** all reads of AB11–14 and AB05.

| method | chimeras flagged | clean pure-culture reads flagged | notes |
|---|---|---|---|
| `vsearch --uchime_ref`, emu database, verdict Y | 51–83% per mix | 0.00–0.10% per culture (AB07 mix: 8.4%) | needs reads on the database strand; vsearch `--orient` oriented only 54%, so minimap2 was used |
| same, Y or borderline | 83–94% | up to 30% (AB04), **18.5% of pure AB05** | unusable |
| `vsearch --uchime_ref` against the sample's ASVs | 54–85% | not applicable | needs savont, see Q6 |
| read two-parent vs ASVs, ASV-only comparator, gain ≥ 15 | 83% | 0% if the species has an ASV, **6.0% if it has none** | the previous approach |
| same, comparator = best ASV or best database hit, gain ≥ 15 | 76% | 0% either way | gain ≥ 10: 82% / 0% |
| `vsearch --uchime_denovo` (dereplicated reads) | 0 of 1,497 | 143 of 14,881 | reads are 99.3% unique, so abundance skew never applies |
| `vsearch --uchime3_denovo` | 0 | 0 | same reason |
| `vsearch --chimeras_denovo` | 79% | **3.0%** | slow (over 10 min for 16,500 reads) |
| emu `--min-pid` 90–98 | – | – | doesn't remove chimera calls; at 90–98 it adds or reshuffles false calls (*S. lutetiensis*, FDAARGOS) and loses real reads |

Mix AB11 is 5 Enterobacterales, AB12 Staphylococcus/Streptococcus/Enterococcus. The real mixes give a share of reads called chimeric by uchime_ref that fits their make-up: AB11 15.5%, AB12 6.7%, AB14 2.8%, AB13 0.5%, pure AB05 0.03%.

Detection depends on how many sites the parents differ by on the shorter side of the break (two-parent score, published panel):

| parent differences on weaker side | ≤2 | 3–4 | 5–7 | 8–10 | 11–15 | 16–30 | >30 |
|---|---|---|---|---|---|---|---|
| detected at gain ≥ 15 | 0% | 2% | 15% | 75% | 72% | 99% | 99% |

*K. pneumoniae* and *K. aerogenes* ASVs differ at only 26 of about 1,470 sites, and *S. aureus* and *S. epidermidis* at 20. Most chimeras between them are invisible to any method. That largely explains why AB13's controls showed 0/300 before.

### Literature and documentation
No source tests chimera detection on raw Nanopore reads.
- The [VSEARCH manual](https://torognes.github.io/vsearch/) says the de novo commands need abundance annotations, and `--uchime3_denovo` is designed for denoised amplicons. `--uchime_ref` doesn't need abundances, but parents (or close relatives) must be in the database.
- [USEARCH uchime2_ref](https://drive5.com/usearch/manual/cmd_uchime2_ref.html) assumes error-free queries and complete parents. It advises against uchime in OTU/denoising pipelines because of false positives.
- [Hakimzadeh et al. 2025, PeerJ](https://pmc.ncbi.nlm.nih.gov/articles/PMC12684403/) tested PacBio ITS. `chimeras_denovo` and DADA2 flagged about half of empirical reads, most of them likely false positives; `uchime_denovo` was the most conservative.
- [Karst et al. (bioRxiv)](https://www.biorxiv.org/content/10.1101/645903v1.full-text) used ONT UMI consensus and found no detectable chimeras, but didn't test raw reads.
- [Emu (Curry et al. 2022)](https://pmc.ncbi.nlm.nih.gov/articles/PMC9939874/) has no chimera handling. Its built-in floor is "10 reads for anything larger" than 1,000 reads, which is why calls of 10–15 reads appear at all.

The earlier rule's two thresholds aren't independent. Gain = S_two − max(S_single, S_called), and S_two − S_parent equals the support on the other side. So gain ≥ 15 already forces each side to be favoured by ≥ 15, and the "each side ≥ 5" condition never changed a result in the 2,849 reads scored both ways (`scripts/check_prior.py`). The new scorer gives identical scores to the earlier one on all of them.

## Q2. What each option does to emu's calls on the real samples

`scripts/eval_calls.py`, table `results_calls_real.tsv`. Counts are emu's estimated reads.

Unexpected calls, low-level calls and the AB12 FDAARGOS_522 synonym under each option. Every expected call with ≥200 reads is unaffected by every option and is confirmed by the near-exact rule (first line). AB05 has one ASV species group, so the two-parent filters can't flag anything there (n/a).

Expected calls with >=200 reads: 18, confirmed by the near-exact rule: 18

| sample   | species                        | class      |   emu as run |   uchime_ref filter | two-parent, ASV+db, ≥15   | two-parent, ASV only, ≥15   | two-parent, emu refs as parents, ≥10   |   near-exact reads (≤1%) | near-exact rule   |
|:---------|:-------------------------------|:-----------|-------------:|--------------------:|:--------------------------|:----------------------------|:---------------------------------------|-------------------------:|:------------------|
| AB05     | Escherichia fergusonii         | unexpected |           66 |                  71 | n/a                       | n/a                         | n/a                                    |                       22 | confirmed         |
| AB11     | Shigella flexneri              | unexpected |          316 |                 104 | 111                       | 0                           | 38                                     |                        0 | UNCONFIRMED       |
| AB11     | Enterobacter asburiae          | unexpected |           56 |                   0 | 11                        | 0                           | 0                                      |                        1 | UNCONFIRMED       |
| AB11     | Citrobacter freundii           | unexpected |           48 |                  18 | 0                         | 0                           | 0                                      |                        0 | UNCONFIRMED       |
| AB11     | Klebsiella planticola          | unexpected |           38 |                   0 | 0                         | 0                           | 0                                      |                        0 | UNCONFIRMED       |
| AB11     | Pectobacterium parmentieri     | unexpected |           19 |                   0 | 0                         | 0                           | 0                                      |                        0 | UNCONFIRMED       |
| AB11     | Klebsiella quasipneumoniae     | unexpected |           15 |                   0 | 0                         | 0                           | 0                                      |                        1 | UNCONFIRMED       |
| AB11     | Salmonella enterica            | unexpected |           14 |                   0 | 0                         | 0                           | 0                                      |                        0 | UNCONFIRMED       |
| AB11     | Staphylococcus haemolyticus    | unexpected |           14 |                  14 | 14                        | 14                          | 14                                     |                        8 | confirmed         |
| AB12     | Streptococcus pneumoniae       | expected   |           50 |                  29 | 28                        | 27                          | 35                                     |                       15 | confirmed         |
| AB12     | Streptococcus sp. FDAARGOS_522 | expected   |           43 |                 117 | 141                       | 153                         | 149                                    |                        0 | UNCONFIRMED       |
| AB12     | Staphylococcus haemolyticus    | unexpected |           29 |                   0 | 0                         | 0                           | 0                                      |                        0 | UNCONFIRMED       |
| AB12     | Streptococcus equinus          | unexpected |           25 |                   0 | 0                         | 0                           | 0                                      |                        0 | UNCONFIRMED       |
| AB13     | Klebsiella variicola           | unexpected |           16 |                   0 | 0                         | 0                           | 0                                      |                        0 | UNCONFIRMED       |
| AB14     | Shigella flexneri              | unexpected |           30 |                   0 | 0                         | 0                           | 0                                      |                        0 | UNCONFIRMED       |
| AB14     | Escherichia fergusonii         | unexpected |            0 |                  24 | 27                        | 27                          | 26                                     |                          |                   |

Summary of the 13 unexpected calls in emu's output, plus AB14's new *E. fergusonii*:

| option | unexpected calls removed or labelled | expected calls lost | new false calls |
|---|---|---|---|
| near-exact rule after emu (label) | 11 of 13 (all except AB11 *S. haemolyticus*, AB05 *E. fergusonii*) | 0 (also labels the AB12 FDAARGOS_522 name; *S. agalactiae* itself stays confirmed) | 0, because emu's result is unchanged |
| uchime_ref before emu | 9 of 13 (AB11 keeps *S. flexneri* 104, *C. freundii* 18, *S. haemolyticus*) | 0 | AB14 *E. fergusonii* 24 |
| two-parent, ASV + database, ≥15, before emu | 9 of 13 (AB11 keeps *S. flexneri* 111, *E. asburiae* 11, *S. haemolyticus*) | 0 | AB14 *E. fergusonii* 27 |
| two-parent, ASV only, ≥15, before emu | 11 of 13 (AB11 keeps *S. haemolyticus*) | 0 in the real samples, but see spikes | AB14 *E. fergusonii* 27 |

How the filtered EM goes wrong:
- **New calls in AB14:** removing reads shifts emu's EM between tied references. *S. flexneri* disappears from AB14 and *E. fergusonii* appears instead.
- **AB12 synonym:** FDAARGOS_522 grows from 43 to 117–153 reads and *S. agalactiae* falls from 1,060 to 421–525.
- **Reads lost:** every read filter also removes real reads of the main species. uchime_ref removes 15.5% of AB11, 6.7% of AB12, 2.8% of AB14, 0.5% of AB13 and 0.03% of AB05.


## Q3. Sensitivity

### Thinning real rare species (their real chimeras kept)
`scripts/thin.py`, `thin_uchime.py`. Reads emu assigns to a true species are cut to n, 40 replicates, and emu's EM is rerun. The thinned species keeps the chimera share it really had in that PCR, which a spike of pure-culture reads can't reproduce.

Share of 40 replicates in which the species is still called (emu), still called and confirmed by the near-exact rule (≥2 reads), or called after the uchime_ref read filter or the two-parent read filter (gain ≥15, database comparator) is applied before emu:

| sample, species (its chimera share) | reads kept | emu | + near-exact rule | uchime_ref before emu | two-parent filter before emu |
|---|---|---|---|---|---|
| AB12 *Streptococcus pneumoniae* (44%) | 10 | 78% | 75% | 0% | 0% |
| AB12 *Streptococcus pneumoniae* (44%) | 12 | 100% | 95% | 0% | 0% |
| AB12 *Streptococcus pneumoniae* (44%) | 15 | 100% | 100% | 30% | 20% |
| AB12 *Streptococcus pneumoniae* (44%) | 20 | 100% | 100% | 85% | 95% |
| AB12 *Streptococcus pneumoniae* (44%) | 25 | 100% | 100% | 100% | 100% |
| AB12 *Enterococcus faecalis* (46%) | 10 | 2% | 2% | 0% | 0% |
| AB12 *Enterococcus faecalis* (46%) | 12 | 100% | 95% | 0% | 2% |
| AB12 *Enterococcus faecalis* (46%) | 15 | 100% | 100% | 22% | 8% |
| AB12 *Enterococcus faecalis* (46%) | 20 | 100% | 100% | 78% | 72% |
| AB12 *Enterococcus faecalis* (46%) | 25 | 100% | 100% | 98% | 80% |
| AB12 *Streptococcus pasteurianus* (42%) | 10 | 78% | 70% | 0% | 0% |
| AB12 *Streptococcus pasteurianus* (42%) | 12 | 100% | 92% | 5% | 8% |
| AB12 *Streptococcus pasteurianus* (42%) | 15 | 100% | 98% | 52% | 38% |
| AB12 *Streptococcus pasteurianus* (42%) | 20 | 100% | 100% | 92% | 82% |
| AB12 *Streptococcus pasteurianus* (42%) | 25 | 100% | 100% | 100% | 100% |
| AB13 *Enterococcus faecalis* (5%) | 10 | 2% | 2% | 0% | 0% |
| AB13 *Enterococcus faecalis* (5%) | 12 | 100% | 100% | 90% | 90% |
| AB13 *Enterococcus faecalis* (5%) | 15 | 100% | 100% | 100% | 100% |
| AB13 *Enterococcus faecalis* (5%) | 20 | 100% | 100% | 100% | 100% |
| AB13 *Enterococcus faecalis* (5%) | 25 | 100% | 100% | 100% | 100% |
| AB13 *Streptococcus pasteurianus* (6%) | 10 | 90% | 88% | 62% | 55% |
| AB13 *Streptococcus pasteurianus* (6%) | 12 | 100% | 100% | 100% | 100% |
| AB13 *Streptococcus pasteurianus* (6%) | 15 | 100% | 100% | 100% | 100% |
| AB13 *Streptococcus pasteurianus* (6%) | 20 | 100% | 100% | 100% | 100% |
| AB13 *Streptococcus pasteurianus* (6%) | 25 | 100% | 100% | 100% | 100% |
| AB11 *Morganella morganii* (17%) | 10 | 92% | 90% | 18% | 18% |
| AB11 *Morganella morganii* (17%) | 12 | 100% | 98% | 72% | 72% |
| AB11 *Morganella morganii* (17%) | 15 | 100% | 98% | 95% | 98% |
| AB11 *Morganella morganii* (17%) | 20 | 100% | 100% | 100% | 100% |
| AB11 *Morganella morganii* (17%) | 25 | 100% | 100% | 100% | 100% |
| AB11 *Klebsiella pneumoniae* (13%) | 10 | 5% | 5% | 10% | 5% |
| AB11 *Klebsiella pneumoniae* (13%) | 12 | 78% | 78% | 52% | 65% |
| AB11 *Klebsiella pneumoniae* (13%) | 15 | 100% | 100% | 98% | 92% |
| AB11 *Klebsiella pneumoniae* (13%) | 20 | 100% | 100% | 100% | 100% |
| AB11 *Klebsiella pneumoniae* (13%) | 25 | 100% | 100% | 100% | 100% |
| AB14 *Staphylococcus aureus* (10%) | 10 | 80% | 80% | 57% | 28% |
| AB14 *Staphylococcus aureus* (10%) | 12 | 98% | 98% | 98% | 92% |
| AB14 *Staphylococcus aureus* (10%) | 15 | 100% | 100% | 100% | 100% |
| AB14 *Staphylococcus aureus* (10%) | 20 | 100% | 100% | 100% | 100% |
| AB14 *Staphylococcus aureus* (10%) | 25 | 100% | 100% | 100% | 100% |

Emu itself can't call fewer than 10 reads in these samples (its floor is 10/n of mapped reads). Chimera-prone rare species are the AB12 streptococci and *E. faecalis*, where about 45% of reads score as chimeric, mostly with the staphylococci. For them, both read filters move the limit from about 10–12 reads to about 20. Where the chimera share is low (AB13, AB14), the filters cost little. The near-exact rule costs nothing measurable at 12 reads or more.

*S. agalactiae* in AB12 and *E. coli* in AB14 were also thinned (`results_thin.tsv`). Their detection is limited by emu splitting reads between tied references (FDAARGOS_520/522; *Shigella*), not by chimeras, so they aren't in the table.


### Spiking clean reads of a species absent from the mix
`scripts/spike.py`: pure-culture reads of a species absent from the mix, n = 5–100, 15 replicates.

15 replicates per level; recipients AB11–14; 19 donor/recipient pairs (`results_spike.tsv`). Share of replicates where the spiked species is called:

|   n | called   | p12   | f_uchime_db   | f_dbcomp15   | f_asv15   |
|----:|:---------|:------|:--------------|:-------------|:----------|
|   5 | 0%       | 0%    | 0%            | 0%           | 0%        |
|   8 | 6%       | 5%    | 6%            | 6%           | 6%        |
|  10 | 56%      | 56%   | 65%           | 65%          | 62%       |
|  12 | 97%      | 97%   | 97%           | 97%          | 92%       |
|  15 | 99%      | 99%   | 100%          | 100%         | 94%       |
|  20 | 100%     | 100%  | 100%          | 100%         | 95%       |
|  30 | 100%     | 100%  | 100%          | 100%         | 99%       |
|  50 | 100%     | 100%  | 100%          | 100%         | 100%      |
| 100 | 100%     | 100%  | 100%          | 100%         | 100%      |

- uchime_ref removes 0.0–0.1% of spiked clean reads, and the database-comparator two-parent filter 0%. Neither loses a clean rare species.
- The ASV-only two-parent filter removes 4% of spiked reads on average, but **64% of clean *E. coli* reads spiked into AB12**, which has no Enterobacterales ASV. *E. coli* is then not called at 12–15 reads and is called in only 7% of replicates at 20 reads (80% at 30), against 100% from 12 reads up with every other option.
- At 10 reads some filters call the spike more often than emu does. Removing reads lowers emu's 10/n floor; it isn't a gain in sensitivity.
- The near-exact rule follows emu exactly except *E. hormaechei* spiked at 8 reads into AB13 (93% vs 100%).

Per pair: `report_tables.md`, Spikes.


## Q4. Thresholds

`scripts/thresholds.py`, `results_thresholds.tsv`. Share of reads flagged:

Two-parent read scoring, published savont ASVs as parents. "Clean reads" are pure-culture reads, either of a species that has ASVs in the panel or of one that has none (the situation of a rare real species):

| comparator | gain ≥ | simulated chimeras | clean reads, species in panel | clean reads, species not in panel | AB07 (2-species PCR) |
|---|---|---|---|---|---|
| ASVs only | 5 | 93.5% | 1.99% | 29.02% | 43.6% |
| ASVs only | 10 | 88.5% | 0.00% | 11.02% | 8.1% |
| ASVs only | 15 | 82.9% | 0.00% | 6.03% | 5.5% |
| ASVs only | 20 | 74.0% | 0.00% | 3.13% | 4.8% |
| ASVs only | 25 | 69.8% | 0.00% | 0.95% | 4.5% |
| ASVs only | 30 | 65.9% | 0.00% | 0.51% | 4.2% |
| ASVs only | 40 | 61.1% | 0.00% | 0.17% | 3.6% |
| ASVs only | 60 | 52.8% | 0.00% | 0.01% | 2.8% |
| ASVs or database | 5 | 87.8% | 0.12% | 0.00% | 5.2% |
| ASVs or database | 10 | 81.7% | 0.00% | 0.00% | 4.5% |
| ASVs or database | 15 | 75.6% | 0.00% | 0.00% | 4.3% |
| ASVs or database | 20 | 69.2% | 0.00% | 0.00% | 4.2% |
| ASVs or database | 25 | 65.8% | 0.00% | 0.00% | 4.1% |
| ASVs or database | 30 | 63.1% | 0.00% | 0.00% | 3.8% |
| ASVs or database | 40 | 59.5% | 0.00% | 0.00% | 3.4% |
| ASVs or database | 60 | 52.0% | 0.00% | 0.00% | 2.7% |

- With an ASV-only comparator, gain ≥ 15 removes 6% of a rare species' clean reads when that species has no ASV. Only gain ≥ 40–60 makes this negligible, and by then 39–48% of chimeras are missed.
- Adding the database's best single reference as a comparator brings that to 0% at every gain tested (0.12% at 5), so a lower gain can be used. Gain ≥ 10 catches 82% of chimeras.
- So ≥ 15 is defensible only with a database comparator. ≥ 5 per side is redundant, as shown in Q1.
- None of these thresholds fixes what matters at the call level: the read loss for chimera-prone rare species shown in Q3.

Near-exact reads (score ≥ (1 − T) × read length to the best reference), the basis of the post-emu rule:

| T | clean pure-culture reads | simulated chimeras |
|---|---|---|
| 0.25% | 9.6% | 0.3% |
| 0.5% | 31.5% | 1.6% |
| 0.75% | 51.9% | 4.0% |
| 1% | 64.7% | 7.0% |
| 1.5% | 82.9% | 15.8% |
| 2% | 90.1% | 22.9% |
| 3% | 95.7% | 33.1% |

At T = 1%, 65% of clean reads are near-exact against 7% of simulated chimeras. Chimeras that pass have their break near an end or nearly identical parents. A chimera-driven call of 15–60 reads can therefore carry 1–4 near-exact reads by chance. In the 11 unexpected calls the rule labels, the count was 0 for 9 and 1 for 2 (*E. asburiae*, *K. quasipneumoniae*), hence ≥ 2. A 12-read real species with half its reads chimeric expects about 0.5 × 12 × 0.65 ≈ 4.

T = 0.5% also separates the real calls (every chimera-driven call has 0) but loses some sensitivity in thinning (`results_thin.tsv`, column perfect05).


## Q5. Cost

Times are wall clock on this 16-core machine, which other threads were also using, so they are upper bounds. Emu figures come from the pipeline trace (12 threads).

| step | per sample | memory | needs |
|---|---|---|---|
| emu (for scale) | 39–104 min, 12 threads (8–21 CPU-hours); in AB12, 98% of it is minimap2 | 2.9–4.9 GB | – |
| **near-exact rule after emu** | **4 s, 1 thread** (AB12 real SAM) | 24 MB | emu `--keep-files`: a 0.4–0.6 GB SAM per sample, which can be deleted after the check; a small python script (stdlib only, runs in the emu container); a new local module |
| uchime_ref before emu | 7–9 min, 4 threads (about 22–36 CPU-min, 2–7% of emu) + orientation 51 s | 0.7 GB | vsearch container (not in TRANA today); an orientation step (vsearch `--orient` failed on 46% of reads, so minimap2/mappy against a few references); 2 new local modules |
| two-parent read scoring before emu | 1–3 min, 4 threads + savont ASV (13–20 s) + classify | 0.3 GB (savont 1.1 GB) | savont would have to run before emu instead of in parallel; a custom numba script; the database-comparator version also needs emu's alignments, so it can't run before emu without a second database mapping |
| `--chimeras_denovo` | over 10 min for 16,500 reads | – | rejected (3% false flags) |

Complexity:
- **Near-exact rule:** one script reading files emu already writes.
- **uchime_ref:** two new processes and a container. It changes emu's input, so every abundance changes, and the run report must say which reads were removed.
- **Two-parent filter:** ties emu's input to savont's output, so it inherits savont's bugs and missing ASVs (Q6).

## Q6. Dependence on savont

Three parent panels were compared:
- **published:** savont 11a4c43 with the reversed-key chimera bug. I reran savont from `savont-gtdb` and it reproduced the published ASVs exactly.
- **patched:** a snapshot of the uncommitted fix in worktree `chimera-parent-similarity`, built in `savont_patched/`, `src/chimera.rs` sha256 5304226d…. It is work in progress by another thread, not a released fix.
- **emupanel:** no savont at all. For each species emu calls at ≥1%, the 6 database references that take most of that sample's primary alignments.

| | AB11 | AB12 | AB13 | AB14 |
|---|---|---|---|---|
| ASVs, published | 17 | 12 | 14 | 13 |
| ASVs, patched | 18 | 18 | 15 | 15 |

The fix changes little here:
- The patched panel gives the same call-level outcome as the published one in every mix (`results_calls_real.tsv`, `asv2p_*[patched]`).
- It catches 4–6 points fewer chimeras than the published panel at gain ≥ 15 (asv-only: 76.9% vs 82.9%), with the same 0% of in-panel clean reads (`report_tables.md`, Thresholds).
- The savont-free emupanel performs about as well (ASV-or-database comparator, gain ≥ 10: 80% of chimeras, 0% of clean reads flagged). So a read-level two-parent filter doesn't need savont at all.

What savont can't supply is a parent that never became an ASV, which is the normal state of a rare species (`scripts/missing_parent.py`). I removed one species' ASVs from the panel and re-scored the simulated chimeras that have that species as a parent:

| | gain ≥ 10, ASV only | gain ≥ 10, ASV or database | gain ≥ 15, ASV or database |
|---|---|---|---|
| all parents in panel | 89% | 82% | 76% |
| one parent missing | 66% | 46% | 43% |

The worst case is AB14 *E. coli* × staphylococcus chimeras with the *E. coli* ASVs removed: 1–2% detected.

Its counterpart is the 6% of a rare species' clean reads that the ASV-only comparator removes when that species has no ASV (Q4). A savont-dependent pre-emu filter is weakest exactly where the low-level calls are. `--uchime_ref` and the near-exact rule both use the full emu database and don't depend on savont.

## Still uncertain

- **One run.** The ≥ 2 near-exact reads at ≤ 1% threshold was chosen on this run's 11 labelled unexpected calls, and the margin is thin: two of them (*E. asburiae*, *K. quasipneumoniae*) had 1 near-exact read. It needs checking on at least one more run with known composition before any clinical use, ideally a run with lower-quality reads.
- **Read quality.** The rule depends on how many reads are near-exact, and that falls steeply with read quality. Of clean reads, 58–92% are near-exact at Q22 and above, but only 8–12% at Q12–20 (`results_quality_dependence.tsv`). On a markedly worse run (older chemistry or a faster basecalling model), real rare species would often be labelled unconfirmed. A run-level guard is needed and is not built: report the near-exact share of the main species, and don't apply the rule below about 20%.
- **Species far from every reference.** Calls of organisms without a reference within about 1% will be labelled unconfirmed. In the leave-species-out test, *M. morganii* reads with *Morganella* removed from the database went to *Providencia rettgeri* and were all unconfirmed, which is the desired outcome. A real strain with an unusual 16S would be labelled the same way. This run's strains are 0–5 mismatches from a reference, so this wasn't seen here.
- **Close relatives.** The rule can't tell apart species whose 16S is nearly identical. *S. pneumoniae* reads with *S. pneumoniae* removed were confirmed as *S. mitis*. Ties such as *E. coli*/*E. fergusonii* in AB05 pass. That needs tie-aware reporting, a separate problem.
- **Not tested here:**
  - **Low-biomass samples.** Chimera formation is known to depend on template amount and cycle number, and only a 15-read NTC was looked at. In the NTC the rule confirms *C. acnes* (5 near-exact reads) and *S. agalactiae* (4), so it doesn't replace NTC/blank comparison.
  - **Other read lengths and kits.**
  - **Savont's patched chimera filter as it will finally be merged.** I used a work-in-progress snapshot.
- **Undetectable chimeras.** Chimeras with ≤ 4 parent differences on one side are almost never detected by any method (≤ 2% at gain 15). Between the closest pairs here (*K. pneumoniae*/*K. aerogenes* 26 sites, *S. aureus*/*S. epidermidis* 20 sites) such chimeras are common. They can only produce calls inside those species complexes.
- **In the pipeline, but only as a separate table.** `--check_emu_calls` adds the prototype's rule as `bin/check_emu_calls.py`, together with a run-level guard: when under 20% of the main species' reads are near-exact, every call in the sample is `not_assessed`. On this run the main species' share was 0.42–0.66. The 20% cut-off follows from the quality-dependence table above but hasn't been tested on a poor run. The labels are in `results/emu_call_check.tsv`, the classifier comparison (`emu_call_check`) and a card in the master report.
  - Checked with a pipeline run of AB12 and the NTC from raw reads (`pipeline_check/run.sh`). Emu's counts were identical to the full run. Every label matched this evaluation: in AB12, *S. pneumoniae* is confirmed with 15 near-exact reads, and FDAARGOS_522, *S. equinus* and *S. haemolyticus* are unconfirmed. In the NTC, *C. acnes* (5) and *S. agalactiae* (4) are confirmed and *S. hominis* is unconfirmed. The check took 6 s and 28 MB for AB12.
  - The near-exact counts of four large AB12 calls were 1–17 reads higher than in the tables above. The pipeline counts the 35 reads that score exactly 99% of their length, which the evaluation's float comparison (1 − score/L ≤ 0.01) dropped. Read by read, everything else agreed (`scripts/diag_scores.py`). Across AB05 and AB11–14 this changes no label, and no count below 5 (`scripts/boundary_effect.py`).
  - With the check off, emu runs without `--keep-files` as before, and nothing else changes (NTC run, `pipeline_check/off`).
  - Emu's SAM is converted to BAM (pysam, in the emu container) inside the emu task, and the SAM is removed unless `--keep_files` is set. For AB12 the BAM is 70 MB. On a 4,000-read AB12 subset, the SAM was 51.6 MB and the BAM 8.8 MB, and the check gave identical output from the SAM and the BAM (`scripts/sam_vs_bam.sh`, `scripts/bam_equivalence.sh`). A second pipeline run with the BAM gave the same labels as the first, and the check took 7 s and 46 MB for AB12. The BAM stays in the work folder until that is cleaned; Nextflow can't delete a task's output once a later task has read it.
- **The EM copy.** `emlib.py` reproduced emu exactly on the 5 published samples and the NTC. Its agreement on filtered input was checked once against a real emu run (AB12, uchime-filtered; see "Method used to test filters cheaply").

## How to rerun

All paths are under `/home/marc/Projects/2026-07-02/chimera_step_eval/`. The environment is `pixi.toml`/`pixi.lock` (vsearch 2.32.0, minimap2 2.30, numba, mappy). The full run's `results/` and `work/` are only read.

1. `scripts/align_all.sh <samples>`: emu's minimap2 command, about 40–100 min per sample at 12 threads. Then `scripts/queue_after_ab05.sh` for the simulated chimeras and pure-culture subsets (`donors/`, made with `seqkit sample -s 11 -n 3000`).
2. `scripts/validate_em.py`: checks `emlib.py` against the published emu tables.
3. `scripts/make_sim.py 1500`, then `scripts/split_sim.py`: real-read chimeras.
4. `scripts/run_savont.sh` (baseline, patched), then `scripts/make_panels.py`: the parent panels.
5. Scoring and uchime: `scripts/run_score_published.sh`, `run_score_panel.sh`, `score_sets.sh`, `run_uchime.sh`, `denovo_test.sh`.
6. Analyses: `eval_calls.py`, `postrules.py`, `thin.py`, `thin_uchime.py`, `spike.py`, `thresholds.py`, `missing_parent.py`, `quality_dependence.py`, `loso.py`, `pure_subsets_rule.py`, `sam_size.py`.
7. Real emu check: `scripts/real_emu.sh`, then `validate_real_filtered.py` and `prototype_near_exact.py`.
8. `scripts/report_tables.py`: writes `report_tables.md`, from which the tables here are copied.

Deleted after use to save space, and regenerated by the scripts above:
- the oriented reads for uchime;
- savont's temp folders and the patched build directory (the binary is kept in `savont_patched/bin`);
- the leave-species-out database copy;
- the real emu SAM (397 MB) and its filtered input.
