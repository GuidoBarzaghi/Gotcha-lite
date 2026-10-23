# Gotcha-lite

A fork of [**landau-lab/Gotcha**](https://github.com/landau-lab/Gotcha) — the R package for
[GoT-ChA](https://www.biorxiv.org/content/10.1101/2022.05.11.491515v1) (Genotyping of Targeted loci with
single-cell Chromatin Accessibility), by Franco Izzo and Robert Myers.

This fork branches from the upstream `development` branch at commit
[`445546e`](https://github.com/landau-lab/Gotcha/commit/445546eee600030f550f8167ffd68646420286ba) and changes
three files. Everything else — the method, the remaining functions, the pipeline vignette — is upstream's;
see the [original README](https://github.com/landau-lab/Gotcha#readme) and
[`GoTChA_pipeline.Rmd`](GoTChA_pipeline.Rmd) for the package overview and usage.

```r
devtools::install_github("GuidoBarzaghi/Gotcha-lite")
```

The changes target the two steps that dominate wall time and memory on large GoT-ChA libraries processed on a
Slurm cluster: fastq filtering and whitelist barcode matching.

## Changes relative to `landau-lab/Gotcha@development`

### `R/FastqFiltering.R` — re-runs no longer append, and memory no longer scales with library size

- **Append-on-rerun fixed.** Upstream created `Filtered/<chunk>/` with a plain `mkdir` and wrote with
  `writeFastq(..., mode = "a")`. Re-running the function over a directory that already held filtered output
  appended a second copy of every read to the existing files, silently inflating read counts downstream. The
  chunk directories are now cleared (`unlink(..., recursive = TRUE)`) and recreated with `mkdir -p` before
  writing, and the fastqs are written with `mode = "w"`.
- **Streaming instead of load-everything.** Upstream read every unfiltered fastq into a list (`orig.files`),
  then built a second list of filtered copies, holding both in memory at once. Each fastq is now read,
  filtered, written and released inside a single `mclapply` pass; only a two-column per-file read count
  (`before.filtering`, `after.filtering`) is accumulated.
- **Metrics computed per file.** The retained-read percentage is derived from those per-file counts, rather
  than from `length(orig.files[[1]])`, which assumed every chunk held the same number of reads.

### `R/MatchBarcodes.R` — new function: segment-indexed Hamming matching against the whitelist

Matching each cell barcode by computing `stringdist()` against the entire whitelist is O(reads x whitelist)
and is the bottleneck of `MutationCalling` on a full library. `MatchBarcodes()` replaces it with a
pigeonhole-principle prefilter: a barcode of length L is split into `k = max_dist + 1` consecutive segments,
the whitelist is indexed by the exact sequence of each segment, and the Hamming distance is computed only
against whitelist entries that share at least one exact segment with the query. Any barcode within
`max_dist` mismatches must share one, so the prefilter is exact, not heuristic.

- `max_dist` defaults to `1`, matching the cellranger-atac convention.
- A match is returned only when it is unique at the smallest observed distance; ties and misses return
  `NA_character_`.
- Returns the input data frame with a `WhiteListMatch` column added. Requires columns `CB`, `R1`, `R2`.
- `ncores > 1` splits the barcodes across `mclapply` workers and restitches them in order; `ncores = 1`
  takes a serial path that reports progress every 10^6 barcodes.
- The function is exported (the package uses `exportPattern`) but has no `.Rd` page yet.

### `R/MutationCalling.R` — barcode matching delegated to `MatchBarcodes()`

The inline per-barcode matching block is replaced by a single `MatchBarcodes(max_dist = 1, ncores = ncores)`
call. Two behavioural consequences worth knowing:

- **Mismatch tolerance is 1, not 2.** Upstream's nested conditionals accepted a unique match at Hamming
  distance 2 when none existed at 0 or 1.
- **Unmatched reads are `NA`, not sentinel strings.** Upstream wrote `"No Match"` / `"Too many matches"` into
  `WhiteListMatch` and dropped those rows before genotyping. Here they carry `NA` and are *not* dropped, so
  the per-cell count table returned by `MutationCalling()` contains one aggregate row with
  `WhiteListMatch = NA` holding the unassigned reads. Filter it before merging into scATAC-seq metadata.
  The reported barcode-matching percentage uses `sum(!is.na(WhiteListMatch))`.

## Unchanged

No change to the genotype-calling logic, noise correction (`AddGenotypingArchr`, `FilterGenotyping`,
`GotchaLabeling`), batch submission (`BatchMutationCalling`, `MergeMutationCalling`), or downstream analysis
(`DiffLMM`, `DiffPeaks`). Dependencies in `DESCRIPTION` are upstream's.

## Citation and license

Cite the original method and package, not this fork:
[GoT-ChA preprint](https://www.biorxiv.org/content/10.1101/2022.05.11.491515v1) and
[landau-lab/Gotcha](https://github.com/landau-lab/Gotcha). License is upstream's CC BY-NC 4.0
(see [LICENSE.txt](LICENSE.txt)).
