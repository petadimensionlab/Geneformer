# Fix: exact logit ties crash evaluation (`vote()` returned "tie")

`evaluation_utils.vote()` returns the arg-max class index, but on an **exact
float tie** between the top classes it returned the string `"tie"`. That string
then reached sklearn as part of `y_pred`:

```
File "geneformer/classifier.py", line 1348, in evaluate_model
  conf_mat, macro_f1, acc, roc_metrics = eu.get_metrics(...)
File "geneformer/evaluation_utils.py", line 157, in get_metrics
  conf_mat = confusion_matrix(y_true, y_pred, labels=list(labels))
ValueError: Mix of label input types (string and number);
  Got [ 0  1 ... 24] and ['0' ... '24' 'tie']
```

So a **single tied cell** aborts the whole fine-tune run at the very end — the
checkpoint is already written, but the training script exits non-zero and the
evaluation/report tables are never produced. It hit the V2-316M AD_spleen run
(1 epoch, ~2 h) and did not hit the V2-104M run, which had no exact ties.

Fix: break the tie deterministically (lowest class index) and keep the value
numeric. Non-tie behaviour is unchanged.

```bash
cp patches/eval_tie/evaluation_utils.py geneformer_hf/geneformer/evaluation_utils.py
```

Consequence to be aware of: a tied cell is silently counted as a prediction of
the lowest-index tied class instead of being flagged. Exact float ties are rare
(they need two class logits to be bit-identical), so the effect on accuracy /
macro-F1 is negligible — but the count is no longer visible in the output.
