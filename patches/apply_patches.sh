#!/usr/bin/env bash
#
# apply_patches.sh — Geneformer ダウンロード後のパッチ適用
#
# 使い方:
#   ./download.sh
#   ./patches/apply_patches.sh
#
# 環境変数:
#   GENEFORMER_DIR   出力先ディレクトリ (既定: <repo_root>/geneformer_hf)
#
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"
GENEFORMER_DIR="${GENEFORMER_DIR:-$REPO_ROOT/geneformer_hf}"

if [ ! -d "$GENEFORMER_DIR" ]; then
  echo "エラー: $GENEFORMER_DIR が存在しません。先に ./download.sh を実行してください。"
  exit 1
fi

echo "===== Geneformer パッチ適用 ====="
echo "対象: $GENEFORMER_DIR"
echo ""

# 1) device.py をコピー (常に上書き)
echo "[1/4] device.py を配置 ..."
cp -f "$SCRIPT_DIR/device.py" "$GENEFORMER_DIR/geneformer/device.py"
echo "  OK: $GENEFORMER_DIR/geneformer/device.py"

# 2) git patch を適用 (既に当たっている場合はスキップ)
echo "[2/4] コードベースパッチを適用 ..."
PATCH_FILE="$SCRIPT_DIR/geneformer_multibackend.patch"

if git -C "$GENEFORMER_DIR" apply --check "$PATCH_FILE" 2>/dev/null; then
  git -C "$GENEFORMER_DIR" apply "$PATCH_FILE"
  echo "  OK: パッチを適用しました。"
else
  echo "  スキップ: パッチは既に適用済みです (または競合しています)。"
fi

# 3) bf16 パッチ (GF_DTYPE=bf16 のときに必須)
#    これが未適用だと、埋め込み抽出 / ISP の統計出力で
#    TypeError: Got unsupported ScalarType BFloat16 で落ちます。
echo "[3/4] bf16 パッチを配置 ..."
if [ -d "$SCRIPT_DIR/bf16" ]; then
  cp -f "$SCRIPT_DIR/bf16/emb_extractor.py"   "$GENEFORMER_DIR/geneformer/emb_extractor.py"
  cp -f "$SCRIPT_DIR/bf16/perturber_utils.py" "$GENEFORMER_DIR/geneformer/perturber_utils.py"
  echo "  OK: emb_extractor.py / perturber_utils.py"
else
  echo "  スキップ: patches/bf16/ がありません。"
fi

# 4) チェックポイントパッチ (長時間 fine-tune の中断対策)
#    上流は save_strategy="epoch" / save_total_limit=1 をハードコードしており、
#    エポック途中で落ちると完了済みの全ステップを失います。
echo "[4/4] チェックポイントパッチを配置 ..."
if [ -d "$SCRIPT_DIR/checkpoints" ]; then
  cp -f "$SCRIPT_DIR/checkpoints/classifier.py" "$GENEFORMER_DIR/geneformer/classifier.py"
  echo "  OK: classifier.py (ステップ単位チェックポイント + 自動再開)"
else
  echo "  スキップ: patches/checkpoints/ がありません。"
fi

echo ""
echo "===== 完了 ====="
echo ""
echo "利用可能なバックエンド:"
python3 -c "
import sys
sys.path.insert(0, '$GENEFORMER_DIR')
from geneformer.device import get_device
print(f'  検出されたデバイス: {get_device()}')
" 2>/dev/null || echo "  (device.py のインポート確認は .venv 有効化後に試してください)"