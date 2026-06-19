# Create `liu_2025_draft_release`: a fresh paper-release repo with no history

## Context

The current repo `~/github/Code_Liu_2025` has accumulated ~22 GB of LFS objects across 1,182 LFS-tracked files spanning 955 commits — mostly regenerable HTML reports, Claude Code session traces, and superseded analyses. Pushes now fail with a GitHub LFS budget error.

Rather than rewrite history with `git filter-repo`, the chosen approach is to **start a new repo with no history** for the paper release, containing only the files the public needs. This avoids GitHub's orphaned-LFS-object retention problem entirely.

The resulting repo will be small enough that **no LFS is needed** — plain git.

## Destination

`~/github/liu_2025_draft_release/` — new directory, will be initialized as a fresh git repo and pushed to `github.com/steverozen/liu_2025_draft_release` (new empty repo to be created on GitHub).

## What gets copied

The working tree of `~/github/Code_Liu_2025` (current commit on `starting.to.clean.up`), filtered through `rsync` with the exclusions listed below.

## Exclusion list (rsync `--exclude` patterns)

### Already-decided drops (user-confirmed)
- `.git/` — fresh history
- `.claude-trace/` — Claude session logs (~160 MB)
- `.claude/` — local Claude Code settings (not for release)
- `.pixi/` — Python virtual environment, contains ~3 GB GRCh37 references
- `code_for_internal_exploration/msi_study/clip_study/` — ~5 GB of regenerable HTML reports
- `Manuscript_data/Mo_CAP9_analysis_remove_for_release/` — name literally says "remove for release" (~362 MB)
- `old_code/` — superseded code (4.7 MB)
- `vignette/figure/separate_pages/` — regenerable plots (154 MB)
- `vignette/figure/parallel_plots/` — regenerable plots (227 MB)
- `vignette/figure/cache/` — render cache
- `vignette/repertoire_of_indel_mutational_signatures.html` — 109 MB, hosted on Zenodo
- `vignette/repertoire_of_indel_mutational_signatures.pdf` — 82 MB, hosted on Zenodo
- `prev_Manuscript_data/` — stale superseded data (87 MB)
- `tmp/` — transient (86 MB)
- `ID15_ID16/presumed-not-used/` — name says it (LFS PDFs)
- `Assignment_Check/` — sanity-check scratchwork (13 MB)
- `plot_output/old_analysis/` — superseded analysis
- `test_data/` — test catalogs (88 KB but per user)
- `test_output/` — test run outputs (14 MB)

### Standard hygiene drops
- `some_sup_tables/rosetta_stone_full.csv` (2.7 GB — already gitignored)
- `some_sup_tables/rosetta_stone_full_cap9.csv` (1.9 GB — already gitignored)
- `.Rproj.user/`, `.Rhistory`, `.RData`, `.Ruserdata`, `Rplots.pdf` (already gitignored)
- `.DS_Store`, `__pycache__/`, `*.pyc`

### Kept (explicitly per user)
- `code_for_internal_exploration/` — except `msi_study/clip_study/`
- `sigFstudy/`
- `plot_output/` — except `old_analysis/`
- `ID15_ID16/` — except `presumed-not-used/`
- `Manuscript_data/` — except `Mo_CAP9_analysis_remove_for_release/`
- All other code/data directories (`code/`, `script/`, `signature_comparisons/`, `vignette/` minus the items above, `deep_analysis/`, `Degasperi_2022/`, `Code_for_extraction_assignment/`, `some_figures/`, `output/`, etc.)
- `CLAUDE.md` (useful for collaborators)

## Steps

1. **Read-only audit before copy** — print final inclusion list size and a count of files >50 MB after exclusions, to confirm no single file exceeds the 100 MB GitHub limit. Critical files to verify pass this check:
   - `vignette/figure/standalone_data/` (144 KB — fine)
   - Largest in `Manuscript_data/` after exclusions: `finalized_cap9/` (12 MB)
   - Largest expected: nothing >50 MB

2. **Create destination + rsync**
   ```bash
   mkdir -p ~/github/liu_2025_draft_release
   rsync -av --exclude-from=/tmp/release_excludes.txt \
     ~/github/Code_Liu_2025/ ~/github/liu_2025_draft_release/
   ```
   where `/tmp/release_excludes.txt` contains the patterns above.

3. **Write `.gitignore` in new repo** — carry over the originals plus new entries:
   ```
   .Rproj.user
   .Rhistory
   .RData
   .Ruserdata
   Rplots.pdf
   .claude-trace/
   .claude/
   .pixi/
   tmp/
   vignette/figure/cache/
   vignette/figure/parallel_plots/
   vignette/figure/separate_pages/
   vignette/repertoire_of_indel_mutational_signatures.html
   vignette/repertoire_of_indel_mutational_signatures.pdf
   some_sup_tables/rosetta_stone_full.csv
   some_sup_tables/rosetta_stone_full_cap9.csv
   .DS_Store
   __pycache__/
   *.pyc
   ```

4. **`git init` + initial commit**
   ```bash
   cd ~/github/liu_2025_draft_release
   git init -b main
   git add .
   git commit -m "Initial commit: paper-release snapshot of Code_Liu_2025"
   ```

5. **Verify** — before pushing, run:
   ```bash
   du -sh ~/github/liu_2025_draft_release       # expect <200 MB
   find ~/github/liu_2025_draft_release -type f -size +50M -not -path '*/.git/*'   # expect empty
   git -C ~/github/liu_2025_draft_release log --oneline       # one commit
   git -C ~/github/liu_2025_draft_release ls-files | wc -l    # sanity
   ```

6. **Create GitHub repo + push** — interactive step the user runs:
   ```bash
   gh repo create steverozen/liu_2025_draft_release --private --source=. --remote=origin --push
   ```
   (or via web UI + `git remote add origin … && git push -u origin main`)

## What this plan does NOT touch

- The existing `~/github/Code_Liu_2025` repo and its `github.com/liumoLM/Code_Liu_2025` remote are left untouched. The LFS budget error on that repo remains a separate concern — to be resolved by the owner of `liumoLM` purchasing an LFS data pack, OR by eventually deleting that repo entirely once the new release repo is established.

## Verification (end-to-end)

After step 5, the user should:
- `xdg-open ~/github/liu_2025_draft_release` and browse the tree to confirm no surprises.
- Spot-check that `Manuscript_data/` contains everything needed (per CLAUDE.md: `Liu_et_al_*` signature/spectra TSVs, `sample_info.tsv`, `finalized_cap9/`, `COSMIC_v3.5_ID_GRCh37_signatures.tsv`, `Koh_signatures.tsv`, etc.).
- Optionally try rendering one vignette (`vignette/vignette.qmd` or similar) from inside the new repo to confirm relative paths still resolve.
- Confirm `git status` is clean and `git log` shows exactly one commit.
