# Danio AluI satellite atlas

Reproducible sequence, structural and methylation analyses of AT-rich and GC-rich AluI satellite regions. The main reference is GRCz12tu **GCF_049306965.1**. Family-assigned methylation, independently anchored array estimates and regional reference-coordinate tracks use different read sets and weighting schemes.

## Reproduction

Install the Python packages in `requirements.txt`; external tools and input accessions are listed in `metadata/software.json` and `metadata/public_inputs.tsv`. Plotting uses an installed CJK font (Droid Sans Fallback) plus DejaVu Sans.

```sh
python reproduce.py verify
python reproduce.py metrics
python code/regional_methylation/reproduce_windows.py
python code/regional_display/reproduce.py --output /tmp/danio-regional-reproduction
python reproduce.py plot tu
python reproduce.py plot aesculapii
python reproduce.py plot kyathit
python code/figures/render_sensitivity.py
```

The regional reproduction command reconstructs CpG-equal estimates and distinct window read counts from deposited valid read–CpG observations, verifies the 10 kb baseline and all compared scales, and writes to a separate directory. Add `--plots` to render the structural panels there. It does not realign reads or relax alignment or modification filters. Raw-input decoding and the targeted competing-alignment example remain in `code/regional_methylation/`.

The full parameter atlases can be redrawn after the regional reconstruction:

```sh
DANIO_COMPARISON_DIR=/tmp/danio-regional-reproduction python code/regional_display/summarize_plot.py
```

The deposited `results/` tables include the complete comparison. For reference-CpG census checks, `code/regional_display/reference_coverage.py` accepts `DANIO_REFERENCE_FASTA` pointing to the original GCF_049306965.1 FASTA. The formal figure script uses the deposited chosen-window table; it does not alter the input observations.

The mapping option `--secondary=yes` begins with two ASCII hyphens.

## Regional tracks

`data/regional_display/regional_windows_20kb_step5kb.tsv.gz` supplies Fig. 4 and S4: 20 kb windows / 5 kb steps, at least 10 valid reference CpGs and 3 distinct read IDs per window. At least 5 reads marks a subset of the **same** CpG-equal mean. Missing windows break the line; measured zero remains valid. Isolated windows are retained for each displayed condition. Coverage strips show CpG counts and read counts; spatial flags identify a longest unobserved interval exceeding half the actual window or no observed CpG in its central half. These are display sensitivity choices, not biological standards.

Regional tracks describe methylation within each array and its adjacent sequence. Long missing intervals limit the continuity of these profiles. In the chromosome 21 AT-rich array, local decreases were observed in some individuals and differed in position.

The 10 kb baseline remains in `data/regional_methylation/windows.tsv`, corresponding to S31/S51. Its zero-placeholder plotting field is never treated as measured zero. S53 distinguishes model-name-limited observations from the total valid observations of that baseline window. Family and independently anchored estimates (Fig. 3) retain their original data.

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
