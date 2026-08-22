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

  echo "git-lfs を使用せず直接 HTTPS でダウンロードします。"
  local base_url="https://${HF_MIRROR}/ctheodoris/Geneformer/resolve/main"

  # モデルディレクトリのファイル
  local model_files=(
    "model.safetensors"
    "config.json"
    "generation_config.json"
    "training_args.bin"
  )
  for f in "${model_files[@]}"; do
    local dest="$GENEFORMER_DIR/$MODEL_DIR/$f"
    if [ -s "$dest" ]; then
      echo "  スキップ (既存): $MODEL_DIR/$f"
      continue
    fi
    echo "  取得: $MODEL_DIR/$f"
    mkdir -p "$(dirname "$dest")"
    curl -L --fail --progress-bar \
      "$base_url/$MODEL_DIR/$f" -o "$dest"
  done

  # パッケージ辞書 (全モデル共通で必要)
  local dict_files=(
    "token_dictionary_gc104M.pkl"
    "gene_median_dictionary_gc104M.pkl"
    "ensembl_mapping_dict_gc104M.pkl"
    "gene_name_id_dict_gc104M.pkl"
  )
  for f in "${dict_files[@]}"; do
    local dest="$GENEFORMER_DIR/geneformer/$f"
    if [ -s "$dest" ]; then
      echo "  スキップ (既存): geneformer/$f"
      continue
    fi
    echo "  取得: geneformer/$f"
    mkdir -p "$(dirname "$dest")"
    curl -L --fail --progress-bar \
      "$base_url/geneformer/$f" -o "$dest"
  done

  echo "HTTPS 取得完了。"
}

# ---------------------------------------------------------------- 実行
# 辞書とコードは常に必要。git-lfs が無い/失敗した場合は HTTPS で取得。
# is_real_pkl: LFS ポインタ(先頭 'version https')は実データ扱いしない
is_real_pkl() {
  local f="$1"
  [ -s "$f" ] && [[ "$(head -c 8 "$f")" != "version " ]]
}

if [ "$FORCE_HTTPS" = true ]; then
  https_download
elif clone_and_lfs_pull; then
  if is_real_pkl "$GENEFORMER_DIR/geneformer/gene_median_dictionary_gc104M.pkl" \
      && is_real_pkl "$GENEFORMER_DIR/geneformer/token_dictionary_gc104M.pkl"; then
    :
  else
    echo "注意: 辞書が LFS ポインタのままです。HTTPS で取得し直します。"
    https_download
  fi
else
  https_download
fi

echo ""
echo "===== 完了 ====="
echo "確認: $GENEFORMER_DIR/$MODEL_DIR/"
ls -la "$GENEFORMER_DIR/$MODEL_DIR/" 2>/dev/null | grep -v '^total' || true
echo ""
echo "モデルファイルのサイズ (不正なら不完全):"
[ -f "$GENEFORMER_DIR/$MODEL_DIR/model.safetensors" ] && \
  echo "  model.safetensors = $(du -h "$GENEFORMER_DIR/$MODEL_DIR/model.safetensors" | cut -f1)"