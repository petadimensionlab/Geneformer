#!/usr/bin/env bash
#
# download.sh — Geneformer に必要なファイルをダウンロードする。
#
# 方針: git / git-lfs を一切使わず、全ファイルを huggingface.co から
# curl で1個ずつ直接取得する。git-lfs や git が無い環境でも確実に動く。
#
# 取得内容:
#   - コード: setup.py, MANIFEST.in, requirements.txt, geneformer/*.py, mtl/*.py
#   - モデル : Geneformer-V2-104M/ (model.safetensors 等)
#   - 辞書   : geneformer/*.pkl (token/median/ensembl/gene_name)
#   - 生成   : pyproject.toml (upstream は未管理のため最小構成を自動生成)
#
# 使い方:
#   ./download.sh                    # geneformer_hf/ に展開
#   ./download.sh --model V1-10M     # 対象モデルを指定
#   ./download.sh --all              # 全モデル
#
# 環境変数:
#   GENEFORMER_DIR   出力先ディレクトリ (既定: <script_dir>/geneformer_hf)
#   HF_MIRROR        Hugging Face ミラー (例: hf-mirror.com)。未指定なら huggingface.co
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
GENEFORMER_DIR="${GENEFORMER_DIR:-$SCRIPT_DIR/geneformer_hf}"
HF_MIRROR="${HF_MIRROR:-huggingface.co}"
HF_REPO="ctheodoris/Geneformer"
BASE_URL="https://${HF_MIRROR}/${HF_REPO}/resolve/main"

MODEL="${MODEL:-V2-104M}"
DOWNLOAD_ALL=false

usage() {
  sed -n '2,34p' "$0"
  exit 0
}

# ---- 引数パース ----
while [[ $# -gt 0 ]]; do
  case "$1" in
    --model)  MODEL="$2"; shift 2 ;;
    --all)    DOWNLOAD_ALL=true; shift ;;
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
[ "$DOWNLOAD_ALL" = true ] && MODEL_DIR="__ALL__"

check_prereq() {
  if ! command -v "$1" >/dev/null 2>&1; then
    echo "エラー: '$1' が見つかりません。インストールしてください。$( [ "$1" = curl ] && echo "" )" >&2
    return 1
  fi
}
check_prereq curl || exit 1

echo "===== Geneformer ダウンロード (curl 直接取得) ====="
echo "出力先   : $GENEFORMER_DIR"
echo "対象モデル: $MODEL_DIR"
echo ""

mkdir -p "$GENEFORMER_DIR"

# is_real: 実データファイルか(LFS ポインタでない + 最小サイズ以上)
is_real() {
  local f="${1:?is_real: file path required}" min_bytes="${2:-1}"
  [ -s "$f" ] && [ "$(stat -c %s "$f" 2>/dev/null || stat -f %z "$f")" -ge "$min_bytes" ] \
    && [[ "$(head -c 8 "$f")" != "version " ]]
}

# fetch: curl で1ファイル取得(既に実データならスキップ、冪等)
fetch() {
  local path="$1" min="${2:-1}"
  local dest="$GENEFORMER_DIR/$path"
  if is_real "$dest" "$min"; then
    echo "  スキップ (既存・実データ): $path"
    return 0
  fi
  echo "  取得: $path"
  mkdir -p "$(dirname "$dest")"
  curl -L --fail --progress-bar "$BASE_URL/$path" -o "$dest"
}

# パッケージコード・設定 (小さく、LFS でないファイル)
code_files=(
  "setup.py:500"
  "MANIFEST.in:50"
  "requirements.txt:50"
  "geneformer/__init__.py:500"
  "geneformer/classifier.py:500"
  "geneformer/classifier_utils.py:500"
  "geneformer/collator_for_classification.py:500"
  "geneformer/emb_extractor.py:500"
  "geneformer/evaluation_utils.py:500"
  "geneformer/in_silico_perturber.py:500"
  "geneformer/in_silico_perturber_stats.py:500"
  "geneformer/mtl/__init__.py:100"
  "geneformer/mtl/collators.py:100"
  "geneformer/mtl/data.py:100"
  "geneformer/mtl/eval_utils.py:100"
  "geneformer/mtl/model.py:100"
  "geneformer/mtl/train.py:100"
  "geneformer/mtl/utils.py:100"
  "geneformer/mtl_classifier.py:100"
  "geneformer/perturber_utils.py:500"
  "geneformer/pretrainer.py:500"
  "geneformer/tokenizer.py:500"
)

# モデル(大きい LFS): デフォルトは対象モデルのみ、--all なら全モデル
model_dirs=()
if [ "$DOWNLOAD_ALL" = true ]; then
  model_dirs=("Geneformer-V1-10M" "Geneformer-V2-104M" "Geneformer-V2-104M_CLcancer" "Geneformer-V2-316M")
else
  model_dirs=("$MODEL_DIR")
fi

# パッケージ辞書 (LFS)。V2 系は gc104M、V1 系は gc30M を使う。
dict_paths=()
for md in "${model_dirs[@]}"; do
  case "$md" in
    Geneformer-V1-10M)
      dict_paths+=(
        "geneformer/gene_dictionaries_30m/token_dictionary_gc30M.pkl:300000"
        "geneformer/gene_dictionaries_30m/gene_median_dictionary_gc30M.pkl:500000"
        "geneformer/gene_dictionaries_30m/ensembl_mapping_dict_gc30M.pkl:1000000"
        "geneformer/gene_dictionaries_30m/gene_name_id_dict_gc30M.pkl:500000"
      ) ;;
    *)
      dict_paths+=(
        "geneformer/token_dictionary_gc104M.pkl:400000"
        "geneformer/gene_median_dictionary_gc104M.pkl:1000000"
        "geneformer/ensembl_mapping_dict_gc104M.pkl:3000000"
        "geneformer/gene_name_id_dict_gc104M.pkl:1000000"
      ) ;;
  esac
done

echo ""
echo "[1/3] コード・設定ファイル"
for entry in "${code_files[@]}"; do
  fetch "${entry%%:*}" "${entry##*:}"
done

echo ""
echo "[2/3] モデル重み ($MODEL_DIR)"
for md in "${model_dirs[@]}"; do
  for entry in \
    "model.safetensors:400000000" \
    "config.json:500" \
    "generation_config.json:80" \
    "training_args.bin:1000"; do
    fetch "$md/${entry%%:*}" "${entry##*:}"
  done
done

echo ""
echo "[3/3] パッケージ辞書 (tokenize/inference に必須)"
# 重複を避ける
seen=""
for entry in "${dict_paths[@]}"; do
  p="${entry%%:*}"; m="${entry##*:}"
  case " $seen " in *" $p "*) continue ;; esac
  seen="$seen $p"
  fetch "$p" "$m"
done

# pyproject.toml (uv の -e install に必要。upstream 未管理のため最小構成を生成)
if [ ! -f "$GENEFORMER_DIR/pyproject.toml" ]; then
  echo "  [pyproject.toml] 最小構成を生成 ..."
  cat > "$GENEFORMER_DIR/pyproject.toml" <<'EOF'
[build-system]
requires = ["setuptools>=61"]
build-backend = "setuptools.build_meta"
EOF
fi

# 整合性チェック
verify_all() {
  MISSING=""
  local file min
  for entry in \
    "setup.py:500" \
    "pyproject.toml:50" \
    "geneformer/__init__.py:500" \
    "geneformer/tokenizer.py:500"; do
    file="${entry%%:*}"; min="${entry##*:}"
    if ! is_real "$GENEFORMER_DIR/$file" "$min"; then
      echo "  [MISSING] $file"
      MISSING="$MISSING $file"
    fi
  done
  for md in "${model_dirs[@]}"; do
    for entry in "$md/model.safetensors:400000000" "$md/config.json:500" \
      "$md/generation_config.json:80" "$md/training_args.bin:1000"; do
      file="${entry%%:*}"; min="${entry##*:}"
      if ! is_real "$GENEFORMER_DIR/$file" "$min"; then
        echo "  [MISSING] $file"
        MISSING="$MISSING $file"
      fi
    done
  done
  for p in "${dict_paths[@]}"; do
    file="${p%%:*}"; min="${p##*:}"
    if ! is_real "$GENEFORMER_DIR/$file" "$min"; then
      echo "  [MISSING] $file"
      MISSING="$MISSING $file"
    fi
  done
  [ -z "$MISSING" ]
}

echo ""
echo "===== ダウンロード後の整合性チェック ====="
if verify_all; then
  echo "OK: 必要な全ファイルが実データとして揃っています。"
else
  echo ""
  echo "ERROR: 不足しています: $MISSING"
  echo "  ネットワーク/権限を確認し、./download.sh を再実行してください。"
  exit 1
fi

echo ""
echo "===== 完了 ====="
echo "確認: $GENEFORMER_DIR/$MODEL_DIR/"
ls -la "$GENEFORMER_DIR/$MODEL_DIR/" 2>/dev/null | grep -v '^total' || true
echo "model.safetensors = $(du -h "$GENEFORMER_DIR/$MODEL_DIR/model.safetensors" 2>/dev/null | cut -f1)"