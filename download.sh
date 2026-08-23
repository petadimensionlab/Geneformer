#!/usr/bin/env bash
#
# download.sh — Geneformer V2-104M に必要なファイルをダウンロードする。
#
# 取得内容:
#   - コード + V2-104M モデルディレクトリ (git-lfs または直接 HTTPS)
#     - Geneformer-V2-104M/model.safetensors  (417 MB 重み)
#     - Geneformer-V2-104M/config.json
#     - Geneformer-V2-104M/generation_config.json
#     - Geneformer-V2-104M/training_args.bin
#   - パッケージ内のトークン/メディアン辞書 (tokenize / inference に必須)
#     - geneformer/token_dictionary_gc104M.pkl
#     - geneformer/gene_median_dictionary_gc104M.pkl
#     - geneformer/ensembl_mapping_dict_gc104M.pkl
#     - geneformer/gene_name_id_dict_gc104M.pkl
#
# 使い方:
#   ./download.sh                    # geneformer_hf/ に展開 (git-lfs 優先)
#   ./download.sh --model V2-104M    # 対象モデルを指定
#   ./download.sh --force-https      # git-lfs を使わず直接 HTTPS で取得
#   ./download.sh --all              # 全モデル (V1-10M, V2-104M, V2-104M_CLcancer, V2-316M) を取得
#
# 環境変数:
#   GENEFORMER_DIR   出力先ディレクトリ (既定: <script_dir>/geneformer_hf)
#   GF_REPO_URL      Hugging Face リポジトリ (既定: ctheodoris/Geneformer)
#   HF_MIRROR        Hugging Face ミラー (例: hf-mirror.com)。未指定なら huggingface.co
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENEFORMER_DIR="${GENEFORMER_DIR:-$SCRIPT_DIR/geneformer_hf}"
GF_REPO_URL="${GF_REPO_URL:-https://huggingface.co/ctheodoris/Geneformer}"
HF_MIRROR="${HF_MIRROR:-huggingface.co}"

MODEL="${MODEL:-V2-104M}"
DOWNLOAD_ALL=false
FORCE_HTTPS=false

usage() {
  sed -n '2,28p' "$0"
  exit 0
}

# ---- 引数パース ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --model)  MODEL="$2"; shift 2 ;;
    --all)    DOWNLOAD_ALL=true; shift ;;
    --force-https) FORCE_HTTPS=true; shift ;;
    -h|--help) usage ;;
    *) echo "不明な引数: $1"; usage ;;
  esac
done

# モデルディレクトリ名の解決
resolve_model_dir() {
  case "$1" in
    V1|V1-10M|Geneformer-V1-10M)            echo "Geneformer-V1-10M" ;;
    V2|V2-104M|Geneformer-V2-104M)          echo "Geneformer-V2-104M" ;;
    V2-CLcancer|Geneformer-V2-104M_CLcancer) echo "Geneformer-V2-104M_CLcancer" ;;
    V2-316M|Geneformer-V2-316M)             echo "Geneformer-V2-316M" ;;
    *) echo "不明なモデル: $1" >&2; exit 1 ;;
  esac
}

MODEL_DIR="$(resolve_model_dir "$MODEL")"

# 必要なツールの確認
check_prereq() {
  local tool="$1"
  if ! command -v "$tool" >/dev/null 2>&1; then
    echo "エラー: '$tool' が見つかりません。インストールしてください。" >&2
    return 1
  fi
}

echo "===== Geneformer ダウンロード ====="
echo "出力先       : $GENEFORMER_DIR"
echo "対象モデル   : $MODEL_DIR"
echo "リポジトリ   : $GF_REPO_URL"
echo "ダウンロード : $([ "$FORCE_HTTPS" = true ] && echo '直接 HTTPS' || echo 'git-lfs (フォールバック HTTPS)')"
echo ""

mkdir -p "$GENEFORMER_DIR"

# ---------------------------------------------------------------- git-lfs 経由
clone_and_lfs_pull() {
  check_prereq git || return 1
  check_prereq git-lfs || { echo "  git-lfs なし → HTTPS フォールバック" >&2; return 1; }

  echo "[1/2] リポジトリをクローン (git-lfs skip) ..."
  if [ ! -d "$GENEFORMER_DIR/.git" ]; then
    GIT_LFS_SKIP_SMUDGE=1 git clone "$GF_REPO_URL" "$GENEFORMER_DIR"
  else
    echo "  既存のクローンを使用: $GENEFORMER_DIR"
    git -C "$GENEFORMER_DIR" fetch --all 2>/dev/null || true
  fi

  echo "[2/2] git-lfs で対象ファイルを取得 ..."
  # 辞書(geneformer/*.pkl)も pull: モデル同様、トークナイズに必須
  if [ "$DOWNLOAD_ALL" = true ]; then
    git -C "$GENEFORMER_DIR" lfs pull
  else
    git -C "$GENEFORMER_DIR" lfs pull \
      --include="$MODEL_DIR/*" \
      --include="geneformer/*.pkl"
  fi

  echo "git-lfs 取得完了。"
  return 0
}

# ------------------------------------------------------------ 直接 HTTPS 取得
https_download() {
  check_prereq curl || { echo "エラー: curl が必要です。" >&2; exit 1; }
  check_prereq git || { echo "エラー: git が必要です(コード取得のため)。" >&2; exit 1; }

  echo "git-lfs を使用せず直接 HTTPS でダウンロードします。"

  # コード本体(setup.py, pyproject.toml, geneformer/*.py 等)も取得する。
  # git-lfs が使えない環境でも、リポジトリ全体を shallow クローンしてコードを確保し、
  # 重い LFS バイナリだけを curl で取得する。
  if [ ! -d "$GENEFORMER_DIR/.git" ]; then
    echo "  [コード] リポジトリを shallow クローン (git-lfs skip) ..."
    GIT_LFS_SKIP_SMUDGE=1 git clone --depth 1 "$GF_REPO_URL" "$GENEFORMER_DIR"
  fi

  local base_url="https://${HF_MIRROR}/ctheodoris/Geneformer/resolve/main"

  # モデルディレクトリのファイル
  local model_files=(
    "model.safetensors:400000000"
    "config.json:500"
    "generation_config.json:80"
    "training_args.bin:1000"
  )
  for entry in "${model_files[@]}"; do
    local f="${entry%%:*}" min="${entry##*:}"
    local dest="$GENEFORMER_DIR/$MODEL_DIR/$f"
    if is_real "$dest" "$min"; then
      echo "  スキップ (既存・実データ): $MODEL_DIR/$f"
      continue
    fi
    echo "  取得: $MODEL_DIR/$f"
    mkdir -p "$(dirname "$dest")"
    curl -L --fail --progress-bar \
      "$base_url/$MODEL_DIR/$f" -o "$dest"
  done

  # パッケージ辞書 (全モデル共通で必要)
  local dict_files=(
    "token_dictionary_gc104M.pkl:400000"
    "gene_median_dictionary_gc104M.pkl:1000000"
    "ensembl_mapping_dict_gc104M.pkl:3000000"
    "gene_name_id_dict_gc104M.pkl:1000000"
  )
  for entry in "${dict_files[@]}"; do
    local f="${entry%%:*}" min="${entry##*:}"
    local dest="$GENEFORMER_DIR/geneformer/$f"
    if is_real "$dest" "$min"; then
      echo "  スキップ (既存・実データ): geneformer/$f"
      continue
    fi
    echo "  取得: geneformer/$f"
    mkdir -p "$(dirname "$dest")"
    curl -L --fail --progress-bar \
      "$base_url/geneformer/$f" -o "$dest"
  done

  ensure_pyproject_toml

  echo "HTTPS 取得完了。"
}

# ---------------------------------------------------------------- 実行
# is_real: 実データファイルか(LFS ポインタの先頭 'version https' でない + 最小サイズ)
is_real() {
  local f="${1:?is_real: file path required}" min_bytes="${2:-1}"
  [ -s "$f" ] && [ "$(stat -c %s "$f" 2>/dev/null || stat -f %z "$f")" -ge "$min_bytes" ] \
    && [[ "$(head -c 8 "$f")" != "version " ]]
}

# verify_all: 必要な全ファイルが実データとして揃っているか。
# 戻り値 0=全部OK, 1=不足あり。不足ファイルは $MISSING にスペース区切りで格納。
verify_all() {
  MISSING=""
  local file min_bytes
  # 相対パス + 実測サイズ(このリポジトリで確認済み)で検証。サイズ未満は LFS ポインタ/欠損とみなす。
  for entry in \
    "$MODEL_DIR/model.safetensors:400000000" \
    "$MODEL_DIR/config.json:500" \
    "$MODEL_DIR/generation_config.json:80" \
    "$MODEL_DIR/training_args.bin:1000" \
    "geneformer/token_dictionary_gc104M.pkl:400000" \
    "geneformer/gene_median_dictionary_gc104M.pkl:1000000" \
    "geneformer/ensembl_mapping_dict_gc104M.pkl:3000000" \
    "geneformer/gene_name_id_dict_gc104M.pkl:1000000" \
    "setup.py:500" \
    "pyproject.toml:50"; do
    file="${entry%%:*}"; min_bytes="${entry##*:}"
    if ! is_real "$GENEFORMER_DIR/$file" "$min_bytes"; then
      echo "  [MISSING] $file"
      MISSING="$MISSING $file"
    fi
  done
  if [ -n "$MISSING" ]; then
    return 1
  fi
  return 0
}

# ensure_pyproject_toml: uv での `-e` インストールに必要な pyproject.toml を生成する。
# upstream ctheodoris/Geneformer は pyproject.toml を git 管理しておらず、
# setup.py のみを持つ。そのため verify_all が必須としている pyproject.toml は
# clone/https いずれの経路でも生成されず、素の環境では必ず不足する。
# ここで最小構成の pyproject.toml を用意する(存在すれば触らない)。
ensure_pyproject_toml() {
  local dest="$GENEFORMER_DIR/pyproject.toml"
  if [ -f "$dest" ]; then
    return 0
  fi
  echo "  [pyproject.toml] 最小構成を生成 (upstream は未管理のため) ..."
  cat > "$dest" <<'EOF'
[build-system]
requires = ["setuptools>=61"]
build-backend = "setuptools.build_meta"
EOF
}

echo ""
echo "===== ダウンロード後の整合性チェック ====="
if verify_all; then
  echo "OK: 必要な全ファイルが実データとして揃っています。"
else
  echo ""
  echo "⚠ 不足/不正ファイル:$MISSING"
  echo "  不足があれば HTTPS で取得します。"
  if [ "$FORCE_HTTPS" = true ]; then
    https_download
  else
    FORCE_HTTPS=true
    https_download
  fi
  if verify_all; then
    echo "OK: 再取得後、必要な全ファイルが揃いました。"
  else
    echo ""
    echo "ERROR: 再取得後も不足しています。"
    echo "  不足ファイル:$MISSING"
    echo "  ネットワーク/権限を確認し、./download.sh --force-https を再実行してください。"
    exit 1
  fi
fi

echo ""
echo "===== 完了 ====="
echo "確認: $GENEFORMER_DIR/$MODEL_DIR/"
ls -la "$GENEFORMER_DIR/$MODEL_DIR/" 2>/dev/null | grep -v '^total' || true
echo ""
echo "モデルファイルのサイズ (不正なら不完全):"
[ -f "$GENEFORMER_DIR/$MODEL_DIR/model.safetensors" ] && \
  echo "  model.safetensors = $(du -h "$GENEFORMER_DIR/$MODEL_DIR/model.safetensors" | cut -f1)"