# Training chain (stages 7-10): benchmark strip -> Node2Vec -> UPNA-DDAL decoder -> evaluation

Run in order, from this directory:

```
python 01_strip_benchmark.py          # -> ../../data/zenodo/network/G_train.pickle
python 02_train_node2vec.py --seed 4707   # -> ../../data/zenodo/models/node2vec_s4707_rebuilt.npz
python 03_train_decoder.py --seed 4707 --embedding ../../data/zenodo/models/node2vec_s4707_rebuilt.npz
python 04_evaluate_ensemble.py         # scores shipped IxIDN+/- pairs with the FROZEN ensemble
```

Repeat `02`/`03` for each of the 5 manuscript seeds (4707, 11, 202, 3407, 8191) to retrain the
full ensemble from scratch. `04_evaluate_ensemble.py` always scores against the shipped frozen
weights in `../../data/zenodo/models/` (the reported numbers come from those exact files, never
from a local retrain); point it at your own retrained `*_rebuilt.*` files only if
you are checking retraining fidelity, not reproducing the reported numbers.

`03_train_decoder.py --pos-weight` defaults to the literal value recorded in
`FROZEN_MECHDIST_MODEL.json`; see `../../data/zenodo/models/MANIFEST.md` for an unresolved
discrepancy between that recorded default and the frozen model's own prose disclosure.

`01_strip_benchmark.py` and `04_evaluate_ensemble.py` read the IxIDN+/IxIDN- edge tables from
`../../data/02S_build_IxIDN/` and `../../data/02S_build_IxIDN_neg/`.
