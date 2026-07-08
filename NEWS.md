# Version 1.1.0 (2026-07-08):

## New features
* `read_metal()` now accepts a z-score column via `zscore_col`, computing a
  numerically stable two-sided -log10(p) that never underflows to `Inf` for
  extremely significant values (#9).
* `read_metal()` now accepts a pre-computed `-log10(p)` column via `logp_col`,
  and p-values that underflow to `0` are floored (with a warning) instead of
  producing `Inf` (#22).
* `locuscompare()` gains an `ld` argument for supplying your own LD table
  (`SNP_A` / `SNP_B` / `R2`), bypassing the reference database -- useful for
  custom or multi-ancestry reference panels (#12, #20).

## Bug fixes
* The lead SNP is now reliably drawn as the purple diamond; a row-ordering bug
  in `assign_color()` could color the wrong SNP (#30).
* Fixed the broken `get_lead_snp()` / `get_position()` examples and the
  malformed `DESCRIPTION` `Authors@R` field so the package passes `R CMD check`
  (#13).

## Infrastructure
* The LD and position database has moved to a new server; reinstall to pick up
  the new connection settings (#35).

# Version 1.0.0 (2019-05-02):
Much faster LD calculation (querying a MySQL database).

# Version 0.1.0:
Calculating LD on the fly using PLINK and 1000 Genomes files.
