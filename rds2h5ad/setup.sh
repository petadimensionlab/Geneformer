#!/usr/bin/env bash
# ============================================================================
# ADPD Geneformer pipeline — environment setup
#
# Builds everything needed to run:
#   01_inspect.R  02_export.R  03_map_and_build_h5ad.py
#
# Reference: analysis/mouse2human.md
#   - R >= 4.4 (R 4.6 on this machine, via Homebrew) + Seurat 5, managed by a renv project
#   - Python 3.12 venv (anndata / numpy / pandas / scipy — steps 1-3 are CPU-only;
#     torch/geneformer for steps 4-6 are NOT installed here)
#   - user-level install, no root, idempotent (safe to re-run)
#
# Note on R: the md asks for R 4.4; per project owner R 4.6 is acceptable.
# Posit PPM serves source packages on macOS, so the Seurat stack is compiled
# locally (~30 min first run; needs Apple clang + gfortran).
#
# Usage:
#   bash setup.sh [TOOLS_DIR]             # default: this directory (rds2h5ad)
#
# Afterwards, run the pipeline per tissue (see run.sh):
#   ./run.sh PD_LN                        # 01 -> 02 -> 03 on rds_data/PD_LN
# ============================================================================
set -euo pipefail

log()  { printf '\n\033[1;34m[setup]\033[0m %s\n' "$*"; }
die()  { printf '\n\033[1;31m[setup][ERROR]\033[0m %s\n' "$*" >&2; exit 1; }
start() { printf '[setup] %s ... ' "$*"; }
ok()    { printf 'done (%ss)\n' "$(( $(date +%s) - _t0 ))"; }
t0()    { _t0=$(date +%s); }

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"           # .../Geneformer/rds2h5ad
REPO_ROOT="$(cd "$PIPE_DIR/.." && pwd)"                            # .../Geneformer
ADPD_TOOLS="${1:-$PIPE_DIR}"
PPM_REPO="https://packagemanager.posit.co/cran/latest"             # source packages on macOS (no mac binaries)
R_PKGS=(Seurat SeuratObject Matrix jsonlite)

mkdir -p "$ADPD_TOOLS"
log "ADPD_TOOLS = $ADPD_TOOLS"

# ---------------------------------------------------------------------------
# 1. R base (>= 4.4 from PATH; Apple clang + gfortran required for source builds)
# ---------------------------------------------------------------------------
R_BIN_CAND="$(command -v Rscript || true)"
[ -n "$R_BIN_CAND" ] || die "Rscript not found on PATH (need R >= 4.4, e.g. 'brew install r')"
R_BIN="$(dirname "$R_BIN_CAND")"
R_VER="$("$R_BIN/Rscript" --version 2>&1 | sed -n 's/.*version \([0-9][0-9.]*\).*/\1/p' | head -1)"
log "using R $R_VER at $R_BIN"
awk -v v="$R_VER" 'BEGIN { split(v, a, "."); exit !(a[1] > 4 || (a[1] == 4 && a[2] >= 4)) }' \
  || die "R >= 4.4 required, found $R_VER"
command -v clang >/dev/null 2>&1  || die "clang not found (xcode-select --install)"
command -v gfortran >/dev/null 2>&1 || die "gfortran not found (brew install gcc) — needed by deldir/dotCall64/irlba"

# Homebrew's R bottle bakes the gcc runtime paths of its build time into
# Makeconf (FLIBS/LIBS). When the installed gcc moves on (R built against
# gcc-15, gcc since upgraded to 16), the versioned target lib dir vanishes
# and every C/C++ package fails at link with "library 'emutls_w' not found".
# Repair: symlink the expected version dir to the one actually present.
R_ETC="$(Rscript --no-echo -e 'cat(R.home("etc"))' 2>/dev/null | tail -1 | tr -d "'")"
GCC_VERDIR="$(grep -m1 -oE '(/opt[^ ]*gcc/lib/gcc/current/gcc/aarch64-apple-darwin[0-9.]+/[0-9]+)' "$R_ETC/Makeconf" 2>/dev/null || true)"
if [ -n "$GCC_VERDIR" ] && [ ! -e "$GCC_VERDIR" ]; then
  GCC_PARENT="$(dirname "$GCC_VERDIR")"
  GCC_ACTUAL="$(ls "$GCC_PARENT" 2>/dev/null | grep -E '^[0-9]+$' | sort -n | tail -1)"
  [ -n "$GCC_ACTUAL" ] || die "gcc runtime dir $GCC_VERDIR missing and no versioned dir found in $GCC_PARENT"
  ln -s "$GCC_ACTUAL" "$GCC_VERDIR"
  log "repaired gcc runtime lib path: $GCC_VERDIR -> $GCC_ACTUAL (Homebrew R/gcc version drift)"
fi

# ---------------------------------------------------------------------------
# 2. renv project in $ADPD_TOOLS (Seurat 5 + friends; source from Posit PPM)
# ---------------------------------------------------------------------------
t0; start "installing renv + Seurat stack into $ADPD_TOOLS (source packages from Posit PPM, ~30 min)"
( cd "$ADPD_TOOLS" && "$R_BIN/Rscript" -e '
local({
  ppm  <- "https://packagemanager.posit.co/cran/latest"
  pkgs <- c("Seurat", "SeuratObject", "Matrix", "jsonlite")
  if (!requireNamespace("renv", quietly = TRUE))
    install.packages("renv", repos = ppm, quiet = TRUE)
  if (!file.exists(file.path(getwd(), "renv", "library")))
    renv::init(repos = ppm)
  renv::activate()
  options(repos = c(CRAN = ppm))
  lib  <- renv::paths$lib()
  have <- tryCatch(rownames(installed.packages(lib.loc = lib)),
                   error = function(e) character())
  need <- setdiff(pkgs, have)
  if (length(need)) {
    cat("installing:", paste(need, collapse = ", "), "\n")
    renv::install(need, prompt = FALSE)
  } else {
    cat("all packages already present in renv library\n")
  }
  # hard verification: load everything in a clean library context
  for (p in pkgs) suppressPackageStartupMessages(library(p, character.only = TRUE))
  try(renv::snapshot(prompt = FALSE), silent = TRUE)
  cat("VERIFY-OK R", R.version.string, "Seurat", as.character(packageVersion("Seurat")),
      "SeuratObject", as.character(packageVersion("SeuratObject")), "\n")
})
')
ok
log "renv project ready: $ADPD_TOOLS/renv"

# renv::init writes a .Rprofile that sources "renv/activate.R" by RELATIVE
# path, which only works when R starts inside the project dir. Overwrite it
# with a cwd-independent version: resolve activate.R via RENV_PROJECT (set by
# env.sh), so any R/Rscript invocation in a sourced shell activates the project.
cat > "$ADPD_TOOLS/.Rprofile" <<'EOF'
local({
  prj <- Sys.getenv("RENV_PROJECT", unset = "")
  act <- if (nzchar(prj)) file.path(prj, "renv", "activate.R")
         else file.path(getwd(), "renv", "activate.R")
  if (file.exists(act)) source(act)
})
EOF
log "wrote $ADPD_TOOLS/.Rprofile (cwd-independent renv activation)"

# ---------------------------------------------------------------------------
# 3. Python 3.12 venv (anndata / numpy / pandas / scipy)
# ---------------------------------------------------------------------------
PY_VENV="$ADPD_TOOLS/pyenv"
if [ -x "$PY_VENV/bin/python" ] && \
   "$PY_VENV/bin/python" -c 'import anndata, numpy, pandas, scipy' >/dev/null 2>&1; then
  log "python venv already complete — skipping"
else
  t0; start "creating python venv at $PY_VENV"
  if command -v uv >/dev/null 2>&1; then
    [ -x "$PY_VENV/bin/python" ] || uv venv --python 3.12 "$PY_VENV" >/dev/null
    uv pip install --python "$PY_VENV/bin/python" anndata numpy pandas scipy >/dev/null
  else
    PY312="$(command -v python3.12 || command -v python3.11 || command -v python3.10 || true)"
    [ -n "$PY312" ] || die "uv not found and no python3.10-3.12 on PATH"
    [ -x "$PY_VENV/bin/python" ] || "$PY312" -m venv "$PY_VENV"
    "$PY_VENV/bin/python" -m pip install --quiet anndata numpy pandas scipy
  fi
  ok
fi
"$PY_VENV/bin/python" -c 'import anndata, numpy, pandas, scipy; \
print("VERIFY-OK", __import__("sys").version.split()[0], "| anndata", anndata.__version__, \
"| numpy", numpy.__version__, "| pandas", pandas.__version__, "| scipy", scipy.__version__)'

# ---------------------------------------------------------------------------
# 4. ref/ symlink — 03_map_and_build_h5ad.py looks for
#    Path(__file__).parent.parent / "ref"  i.e. <repo root>/ref,
#    while the MGI table lives in analysis/ref/
# ---------------------------------------------------------------------------
REF_LINK="$REPO_ROOT/ref"
if [ -d "$REPO_ROOT/analysis/ref" ] && [ ! -e "$REF_LINK" ]; then
  ln -s "$REPO_ROOT/analysis/ref" "$REF_LINK"
  log "symlinked $REF_LINK -> analysis/ref (so 03 finds the MGI homology table)"
elif [ -e "$REF_LINK" ]; then
  log "ref location already present at $REF_LINK — leaving as is"
fi

# ---------------------------------------------------------------------------
# 5. Locate a .rds (rds_data/<tissue>/ preferred; input/ as fallback)
# ---------------------------------------------------------------------------
RDS="$(find "$PIPE_DIR/rds_data" -mindepth 2 -maxdepth 3 -name '*.rds' 2>/dev/null | head -1 || true)"
[ -n "$RDS" ] || RDS="$(find "$REPO_ROOT/input" -maxdepth 3 -name '*.rds' 2>/dev/null | head -1 || true)"
if [ -n "$RDS" ]; then
  log "rds: $RDS ($(du -h "$RDS" | cut -f1))"
else
  log "WARNING: no .rds found in rds_data/ or input/ — set RDS manually"
fi

# ---------------------------------------------------------------------------
# 6. Write env.sh (source it in each session)
# ---------------------------------------------------------------------------
GENEFORMER_DIR_CAND="$( [ -d "$REPO_ROOT/geneformer_hf/geneformer" ] && echo "$REPO_ROOT/geneformer_hf" || true )"
[ -n "${GENEFORMER_DIR_CAND:-}" ] || log "WARNING: geneformer_hf checkout not found — set GENEFORMER_DIR manually"

cat > "$ADPD_TOOLS/env.sh" <<EOF
# Generated by analysis/setup.sh on $(date +%Y-%m-%dT%H:%M:%S%z)
export ADPD_TOOLS="$ADPD_TOOLS"
export R_BIN="$R_BIN"
export RENV_PROJECT="\$ADPD_TOOLS"
export R_PROFILE_USER="\$ADPD_TOOLS/.Rprofile"
export RSCRIPT="$R_BIN/Rscript"
export PYTHON="$PY_VENV/bin/python"
export GENEFORMER_DIR="${GENEFORMER_DIR_CAND:-}"
export GENEFORMER_MODEL="Geneformer-V2-104M"
export RDS="${RDS:-}"
export ADPD_ROOT="${ADPD_ROOT:-$PIPE_DIR/ws}"
export ADPD_PREFIX="\${ADPD_PREFIX:-PD_LN}"
EOF
log "wrote $ADPD_TOOLS/env.sh"

# ---------------------------------------------------------------------------
# 7. End-to-end verification of the RSCRIPT entry point
# ---------------------------------------------------------------------------
t0; start "verifying RSCRIPT (renv + Seurat load)"
# Rscript -e applies an extra level of string unescaping, so run the check
# from a temp file (file mode parses escapes standardly).
VERIFY_R="$(mktemp "${TMPDIR:-/tmp}/adpd_verify_XXXXXX.R")"
cat > "$VERIFY_R" <<'EOF'
mm <- unlist(strsplit(sub(".*R version ([0-9.]+).*", "\\1", R.version.string), "\\."))
v <- as.numeric(paste(mm[1:2], collapse = "."))
stopifnot(v >= 4.4)
stopifnot(any(grepl("renv", .libPaths())))
suppressPackageStartupMessages({ library(Seurat); library(SeuratObject); library(Matrix); library(jsonlite) })
cat("RSCRIPT-OK", R.version.string, "| Seurat", as.character(packageVersion("Seurat")), "\n")
EOF
( cd /tmp && \
  ADPD_TOOLS="$ADPD_TOOLS" \
  RENV_PROJECT="$ADPD_TOOLS" \
  R_PROFILE_USER="$ADPD_TOOLS/.Rprofile" \
  "$R_BIN/Rscript" "$VERIFY_R" )
rm -f "$VERIFY_R"
ok

# ---------------------------------------------------------------------------
# 8. Manifest + summary
# ---------------------------------------------------------------------------
cat > "$ADPD_TOOLS/manifest.json" <<EOF
{
  "created": "$(date -u +%Y-%m-%dT%H:%M:%SZ)",
  "adpd_tools": "$ADPD_TOOLS",
  "r": {"binary": "$R_BIN", "version": "$R_VER"},
  "renv": {"project": "$ADPD_TOOLS", "lock": "$ADPD_TOOLS/renv.lock", "repos": "$PPM_REPO"},
  "python": {"venv": "$PY_VENV", "version": "$("$PY_VENV/bin/python" --version 2>&1 | sed 's/Python //')"},
  "geneformer_dir": "${GENEFORMER_DIR_CAND:-}",
  "rds": "${RDS:-}",
  "mgi_sha256": "$( [ -f "$REPO_ROOT/analysis/ref/mgi_homology.rpt.sha256" ] && cat "$REPO_ROOT/analysis/ref/mgi_homology.rpt.sha256" || echo none )"
}
EOF

cat <<EOF

============================================================================
 Environment ready.

 Next session:
   source "$ADPD_TOOLS/env.sh"

  Run the pipeline (from $PIPE_DIR):
    ./run.sh PD_LN                 # 01 -> 02 -> 03 on rds_data/PD_LN
    ./run.sh <tissue> <sample_col> <ct_col> [max_per_group] [assay]
    # outputs land in \$ADPD_ROOT/<tissue>/{logs,export,h5ad}

 Notes:
   - Steps 1-3 are CPU-only; steps 4-6 (torch/GPU, fine-tuning) are NOT set up here.
   - The repo .venv (python 3.12 + geneformer + torch) from previous runs is untouched.
   - RDS was auto-detected; override with: export RDS=/path/to/file.rds
============================================================================
EOF
