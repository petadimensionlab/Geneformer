# Geneformer DirectML バックエンドセットアップ

このドキュメントでは、NVIDIA GPU を搭載していない Windows / WSL2 環境で **DirectML** バックエンドを使用して Geneformer の GPU アクセラレーションを有効にする方法を説明します。DirectML は AMD、Intel、Qualcomm の GPU を DirectX 12 API 経由でサポートします。

## 概要

Geneformer の `device.py` モジュールは、利用可能な最適な計算バックエンドを自動選択します:

1. **DirectML (dml)** — Windows / WSL2、任意の DirectX 12 GPU
2. **CUDA (cuda)** — NVIDIA GPU または AMD ROCm
3. **MPS (mps)** — Apple Silicon
4. **CPU (cpu)** — フォールバック

`GENEFORMER_DEVICE` 環境変数で特定のバックエンドを強制指定できます。

## 前提条件

- **Windows 10/11** + WSL2 (Ubuntu 22.04+ 推奨)
- **DirectX 12 対応 GPU** (AMD、Intel、NVIDIA、Qualcomm)
- Windows ホストに **WSL GPU ドライバ** がインストール済みであること

## インストール手順

### 1. Python 仮想環境の作成

```bash
cd /path/to/Geneformer
python3 -m venv .venv
source .venv/bin/activate
```

### 2. PyTorch (CPU版) と torch-directml のインストール

```bash
# CPU版 PyTorch をインストール (torch-directml のベースとして必要)
pip install torch==2.4.1 torchvision==0.19.1 --index-url https://download.pytorch.org/whl/cpu

# torch-directml のインストール
pip install torch-directml
```

> **注意:** torch-directml 0.2.5.dev240914 は `torch==2.4.1` および `torchvision==0.19.1` と互換性があります。最新の互換バージョンは [PyPI](https://pypi.org/project/torch-directml/) で確認してください。

### 3. Geneformer コードのダウンロードとパッチ適用

```bash
# コードベースとモデルのダウンロード
./download.sh

# マルチバックエンドパッチ + DML対応 device.py の適用
./patches/apply_patches.sh
```

### 4. 動作確認

```python
from geneformer.device import get_device, get_device_obj

print(f"検出されたデバイス: {get_device()}")      # dml と表示されること
print(f"デバイスオブジェクト: {get_device_obj()}") # privateuseone:0 と表示されること
```

## 使い方

### DirectML の自動検出による実行

デフォルトで Geneformer が DirectML を検出して使用します:

```bash
source .venv/bin/activate
python analysis/04_baseline.py
```

### 特定のバックエンドを強制指定

```bash
# CPU を強制
GENEFORMER_DEVICE=cpu python analysis/04_baseline.py

# DirectML を強制
GENEFORMER_DEVICE=dml python analysis/04_baseline.py
```

## 仕組み

### デバイス検出の優先順位 (`geneformer/device.py`)

1. `get_device()` が `torch-directml` の利用可能性を確認 → `"dml"`
2. 次に `torch.cuda.is_available()` を確認 (NVIDIA CUDA / AMD ROCm)
3. 次に Apple MPS を確認
4. 最後に `"cpu"` にフォールバック

### パッチの概要 (`patches/geneformer_multibackend.patch`)

このパッチは9つのファイルにハードコードされた `"cuda"` 文字列を `get_device_obj()`、`move_to_device()`、`empty_cache()`、`manual_seed_all()` の呼び出しに置き換えます。

主な修正内容:

- `device="cuda"` → `device=get_device_obj()`
- `.to("cuda")` → `.to(get_device_obj())`
- `torch.cuda.empty_cache()` → `empty_cache()`
- `torch.device("cuda" if ... else "cpu")` → `get_device_obj()`
- `df["col"][bool_index]` → `df["col"].iloc[bool_index]` (`tokenizer.py` の FutureWarning 修正)

## トラブルシューティング

| 症状 | 原因 | 解決策 |
|------|------|--------|
| `RuntimeError: Found no NVIDIA driver` | CUDA なし、DirectML 未インストール | `torch-directml` をインストール、または `GENEFORMER_DEVICE=cpu` を設定 |
| `ImportError: undefined symbol` | torch-directml と PyTorch のバージョン不一致 | 互換性のあるバージョンをインストール (インストール手順2を参照) |
| `DML device not found` | WSL GPU ドライバ未インストール | Windows ホストに [WSL 用 GPU ドライバ](https://learn.microsoft.com/en-us/windows/wsl/tutorials/gpu-compute) をインストール |
| `privateuseone:0` が遅い | 初回の JIT コンパイル | ウォームアップとして小さな forward pass を先に実行 |
| `Failed to load CPU gemm_4bit_forward` | bitsandbytes CPU フォールバック警告 | 無害な警告です。純粋な PyTorch にフォールバックします |

## バージョン互換性

| コンポーネント | バージョン |
|---------------|-----------|
| Python | 3.10 – 3.12 |
| PyTorch | 2.4.1 (CPU) |
| torch-directml | 0.2.5.dev240914 |
| torchvision | 0.19.1 (CPU) |
| OS | WSL2 上の Ubuntu 22.04+ (Windows 10/11) |