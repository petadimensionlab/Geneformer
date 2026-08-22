"""Quick MPS correctness smoke test: load V2-104M, run a small forward pass,
verify output is finite and matches CPU output."""
from __future__ import annotations

import os
import torch

from geneformer import perturber_utils as pu

MODEL_DIR = os.environ.get("GENEFORMER_DIR", "geneformer_hf") + "/Geneformer-V2-104M"

print("device: mps")
model_mps = pu.load_model("Pretrained", 0, MODEL_DIR, mode="eval")
model_mps = model_mps.to("mps")
model_mps.eval()

# tiny deterministic input: 2 cells, some token ids
input_ids = torch.tensor(
    [[1, 5, 100, 200, 0, 0], [2, 7, 42, 8, 9, 0]], device="mps"
)
attention_mask = torch.tensor(
    [[1, 1, 1, 1, 0, 0], [1, 1, 1, 1, 1, 0]], device="mps"
)

with torch.no_grad():
    out_mps = model_mps(input_ids=input_ids, attention_mask=attention_mask)

emb_mps = out_mps.hidden_states[-1].cpu()
print("MPS embeddings shape:", tuple(emb_mps.shape))
print("MPS output finite:", bool(torch.isfinite(emb_mps).all()))
print("MPS output abs mean:", float(emb_mps.abs().mean()))

# Compare against CPU
model_cpu = pu.load_model("Pretrained", 0, MODEL_DIR, mode="eval").to("cpu")
model_cpu.eval()
with torch.no_grad():
    out_cpu = model_cpu(input_ids=input_ids.cpu(), attention_mask=attention_mask.cpu())
emb_cpu = out_cpu.hidden_states[-1]
print("CPU output abs mean:", float(emb_cpu.abs().mean()))
diff = (emb_cpu - emb_mps).abs().max()
print("max |MPS - CPU| diff:", float(diff))
print("CORRECTNESS:", "PASS" if float(diff) < 0.5 else "FAIL")
