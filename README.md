# Danio AluI satellite atlas

Reproducible sequence, structural and methylation analyses of AT-rich and GC-rich AluI satellite regions. The main reference is GRCz12tu **GCF_049306965.1**. Family-assigned methylation, independently anchored array estimates and regional reference-coordinate tracks use different read sets and weighting schemes.

## Reproduction

Install the Python packages in `requirements.txt`; external tools and input accessions are listed in `metadata/software.json` and `metadata/public_inputs.tsv`. Plotting uses an installed CJK font (Droid Sans Fallback) plus DejaVu Sans.

```sh
python code/figures/render_fig1.py
python code/figures/render_fig2.py
python code/figures/render_fig3.py
python code/figures/render_congeners.py
python data/sequence_sensitivity/code/plot_sequence_sensitivity.py
python reproduce.py verify
python reproduce.py metrics
python code/regional_methylation/reproduce_windows.py
python code/regional_display/reproduce.py --output /tmp/danio-regional-reproduction
python code/genome_methylation/render_regional_500kb.py
python reproduce.py plot aesculapii
python reproduce.py plot kyathit
python code/figures/render_sensitivity.py
```

The regional reproduction command reconstructs CpG-equal estimates and distinct window read counts from deposited valid read–CpG observations, verifies the 10 kb baseline and all compared scales, and writes to a separate directory. Add `--plots` to render the structural panels there. It does not realign reads or relax alignment or modification filters. Raw-input decoding and the targeted competing-alignment example remain in `code/regional_methylation/`.

The full parameter atlases can be redrawn after the regional reconstruction:

```sh
DANIO_COMPARISON_DIR=/tmp/danio-regional-reproduction python code/regional_display/summarize_plot.py
```

The deposited `results/` tables include the window comparison. For reference-CpG census checks, `code/regional_display/reference_coverage.py` accepts `DANIO_REFERENCE_FASTA` pointing to the original GCF_049306965.1 FASTA. The window comparison uses 200 kb flanks (S54, S55 and S8). Figure 4 and S4 use the 500 kb flank data (S65 and S66).

## Regional tracks: 500 kb flanks

`data/regional_500kb/windows.tsv.gz` supplies Fig. 4 and S4, with 500 kb on each side of an array: 20 kb windows / 5 kb steps, at least 10 valid reference CpGs and 3 distinct read IDs per window. At least 5 reads marks a subset of the **same** CpG-equal mean. Missing windows break the line; measured zero remains valid. Isolated windows are retained for each displayed condition. Coverage strips show CpG counts and read counts; spatial flags identify a longest unobserved interval exceeding half the actual window or no observed CpG in its central half. These are display sensitivity choices, not biological standards.

Regional tracks describe methylation within each array and its adjacent sequence. Long missing intervals limit the continuity of these profiles. In the chromosome 21 AT-rich array, local decreases were observed in some individuals and differed in position.

The 10 kb baseline remains in `data/regional_methylation/windows.tsv`, corresponding to S31/S51. Its zero-placeholder plotting field is never treated as measured zero. S53 distinguishes model-name-limited observations from the total valid observations of that baseline window. Family and independently anchored estimates (Fig. 3) retain their original data.

The 500 kb table contains 139158 array, sample and window records. Overlapping windows share CpGs and reads. S65 reports window values and S66 reports display counts. S57 and S62 use the predefined 200 kb backgrounds. The structural matrices in `data/regional_500kb/structure/500kb/` use the same TU reference and ModDotPlot parameters.

## Sequence and comparison records

| Location | Contents |
|---|---|
| `data/five_part/` | Discovery, annotated intervals, array organization and related-species comparisons |
| `data/five_part/16_selected_danio_comparison/consensuses.fa` | Local consensus sequences; IDs correspond to supplementary table S9 |
| `data/cypriniformes_blast/` | Queries, taxonomy and the four Cypriniformes search controls |
| `data/sequence_sensitivity/` | Fixed input coordinates, eligible-copy pools, normalization rules and sensitivity results |
| `data/manuscript_evidence/` | Alternation permutations, copy-filter summaries, congener common-region comparisons and methylation denominators/threshold summaries |
| `data/template_controls/` | Four seed configurations, rerun coverage, block coordinates and subtype comparisons |

In `data/template_controls/`, `coverage_summary.tsv` and `subtype_rotation_summary.json` summarize the rerun controls. Assembly/core directories retain `family_bases.tsv`, `receipt.json`, and the `two_verified.json`, `three_verified.json`, `reverse_verified.json`, and `rot92_verified.json` parameter/validation records. `native_output_index.tsv` maps complete native decomposition tables to repository paths; `subtype_blocks.tsv.gz` retains block-level subtype results. These reruns test boundaries after decomposition, whereas sequence-normalization sensitivity keeps the original decomposition coordinates fixed.

The saved commands use paths relative to each assembly/core directory. `working.fa` is the reference tile defined in `receipt.json`. The four template FASTA files are in `data/template_controls/`.

Source fields prefixed with `archive:` identify original evidence records. The tables and scripts listed above provide the deposited data.

The source scripts in `code/analysis_source/` accept `DANIO_ASSEMBLY_RUN_ROOT`, `DANIO_ANNOTATION_ROOT` and `DANIO_REPEAT_EVIDENCE_ROOT` for external analysis inputs. Their default directories are `inputs/assembly_runs/`, `inputs/annotation/` and `inputs/repeat_evidence/`. To repeat the chromosome 11 record audit using the deposited files, run `python data/manuscript_evidence/chr11_and_units/audit_existing_records.py --output /tmp/danio-chr11-audit`.

## Files and statistical units

| Location | Contents |
|---|---|
| `data/array_study/` | Array membership, structure matrices and sequence metrics |
| `data/regional_display/per_run/` | Unique valid read–reference-CpG observations, original model metadata and reconstruction receipts |
| `data/regional_display/results/` | All five window schemes, direct CpG-equal backgrounds, continuity and local exploratory comparisons |
| `data/regional_display/regional_array_coverage_20kb.tsv` | Array/sample summaries; CpGs, observations and read IDs deduplicated within each plot range |
| `data/regional_display/model_metadata_limited_windows_10kb.tsv` | S53 model-name-limited subset and total observations |
| `data/regional_display/gap_*.tsv` | Target and adjacent-control diagnostic counts, inventories and competition summaries |
| `code/regional_display/` | Aggregation and independent parameter-comparison algorithms |
| `code/figures/` | Equal-scale triangular structure plots and regional tracks |
| `metadata/` | Versions, inputs, filtering conditions and data dictionaries |

Coordinates ending in `start0/end0` are 0-based half-open; the publication uses 1-based inclusive coordinates. `mean` is the mean of valid reference-CpG fractions, while `weighted` pools read–CpG observations. Window read counts concern the whole window, not equal depth at each CpG. Array-internal/background means are calculated directly from CpGs; `displayed_window_equal_mean_pct` instead averages displayed window estimates. Overlapping windows share observations.

A read–CpG observation is not a read, and a read is not a biological replicate. Returned high-scoring alignments provide a restricted coordinate-consistency test, not exhaustive uniqueness. The long-gap diagnostic describes retained alignments; it does not establish physical sequencing coverage. Alternative positions not checked against reference CG are called candidate reference coordinates. Regional variation motivates CENP-A/CENH3 localization rather than defining a functional centromere by curve shape alone.

Structural panels use u=(s+t)/2 and v=|s−t|/2 with equal physical coordinate scales (full triangle base:height 2:1). Each array retains its own full range; congeners have background and detail views. Matrices and array coordinates are preserved.

The summary plotters read frozen tables; they do not rerun genome annotation or methylation analysis. `render_congeners.py` redraws S1–S3 using the assembly chromosome inventory in `data/figure_source/congener_chromosomes.tsv` (the same GCA_903798145.1 and GCA_903798195.1 sources as the original maps). Its size rank uses supported satellite bases; the single 424-bp AT-rich record is plotted at 0.424 kb and integer rank 1. Figure S7 coverage uses all eligible copies divided by the satellite bases in each region; distance estimates use the separately sampled copies.

## Genome-wide TL forebrain methylation

The six TL adult forebrain datasets ERR16782314–ERR16782319 are projected onto the full 26-sequence TU GCF_049306965.1 reference. Nuclear chromosome statistics cover 25 sequences; mitochondrial statistics remain separate. Chromosome-wide projection is distinct from Fig. 3A family-assigned read observations and from independently anchored array estimates. It uses reference-CG agreement, the existing MAPQ, modification-state and returned-candidate coordinate filters. Returned candidates are not an exhaustive enumeration of potential matches.

`data/genome_methylation/` contains individual summaries, 200 kb non-overlapping windows, denominators and display masks. `parameters.json` records the reference, samples, software versions and filters. `column_dictionary.json` links the source columns to S59–S64. Figure 5 shows chromosomes 3 and 21. S9 covers all 25 nuclear chromosomes.

```sh
python code/genome_methylation/verify.py --output /tmp/danio-genome-check
DANIO_FIGURE_OUTPUT=/tmp/danio-genome-check/figures python code/genome_methylation/render_main.py
DANIO_ATLAS_OUTPUT=/tmp/danio-genome-check/atlases python code/genome_methylation/render_atlas.py 200
```

Black lines average valid reference-CpG fractions equally. Orange dashed lines pool valid read–CpG observations; C, 5mC and 5hmC form the denominator. A window is displayed with at least 10 effective reference CpGs and three distinct reads. Isolated valid windows retain dots/crosses, while missing windows break the line. Red ticks mark a longest unobserved interval exceeding half the actual window or no valid CpG in the central half. Coverage strips show the fraction of reference CpGs with a valid observation. The 75%-opacity background joins short same-family interruptions only for drawing; original fragment coordinates and statistical masks remain unchanged.

The 200 kb tables provide CpG and observation summaries, distinct-read bounds, and exact read counts for the 95 windows that required them to determine three-read support. The commands above reproduce figures and check the deposited summaries. Source read accessions and projection settings are listed in `data/genome_methylation/parameters.json`.

In the callable set, GC-rich satellites have lower methylation than AT-rich satellites and sequence outside the two-family catalogue across all six individuals. This supports a hypothesis of differentiated chromatin states and a possible association of GC-rich regions with active centromeric chromatin. It does not establish functional inactivity of AT-rich arrays. TU reference structure and TL forebrain methylation come from different genetic backgrounds; protein occupancy and cytological localization in matched material are needed to assign centromere function.
