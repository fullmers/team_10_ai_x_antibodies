# Run log — HER2 CDRH3 design (HER2 novel)

All work inside `~/amar`. Run folder: `~/amar/data/docs` (the folder holding `Data/` and `Settings/`).
Date: 2026-11-14.

## Environment / versions

| Item | Version |
| --- | --- |
| Host | Linux 6.11.0-1016-nvidia, x86_64 |
| GPU | NVIDIA H100 80GB HBM3 (81559 MiB) |
| Docker | 29.9.2, build a82e528 |
| Design image | `quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09` (image id `49d8ccd71fb4`, built 3 days ago) |
| Design weight | `her2-novel-2026-11-09` (RFD4-Proteina v0.3.0) |
| Folding image | `quay.io/apheris/oslo-hackathon:apheris-boltz-v1.3.0-boltz2` (image id `5ff19499c91b`) — Boltz-2, apheris-boltz v1.3.0 |
| torch (in design image) | 2.14.1+cu130 |
| Ready-made script | `data/docs/Scripts/run-her2.sh` (unmodified) |

## Commands

### 1. Inspection (read-only)

```bash
cd /home/hackathon/amar && ls -la && ls -la data/docs/
cd /home/hackathon/amar/data/docs && cat index.md log.md 12-her2.md 18-vm.md \
  Scripts/run-her2.sh Settings/12-her2.json Data/12-her2/her2_cdrh3.yaml
grep -o '<title>[^<]*</title>' Tutorial.html
docker --version; nvidia-smi --query-gpu=name,memory.total --format=csv; docker images | grep -i hackathon
docker run --rm --entrypoint sh quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09 \
  -c 'python -c "import torch;print(torch.__version__)"'
```

### 2. Design — 20 designs, HER2 novel

Run via the ready-made script, overriding only the batch size (`N`) and the number
refolded (`TOPN`). `IMG` is left at its default, which is already `her2-novel-2026-11-09`.
HER2 novel is used with `sequenceTsrK` unset, per the models table in `12-her2.md`;
the script passes no `--sequence-tsr-k`, so this is correct.

```bash
cd /home/hackathon/amar/data/docs
N=20 TOPN=5 OUT=runs/her2 bash Scripts/run-her2.sh design
```

Effective container command:

```bash
docker run --rm --gpus all --shm-size 16g --user $(id -u):$(id -g) -e HOME=/tmp \
  -v /etc/passwd:/etc/passwd:ro -v $PWD:/work -w /work \
  quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09 \
  module-run design --input Data/12-her2/her2_cdrh3.yaml --output runs/her2/design \
  --num-replicates 20 --seed 1 --recycle-strategy high/front_loaded
```

Checkpoint actually loaded (from `per_example_metrics.csv`, column `ckpt_name`):
`rfd4-proteina-her2-novel-2026-11-09-ema.ckpt` — confirms the requested weight.

Result: exit 0, 20/20 designs written. ~7 min wall clock.

Note: a first attempt failed (exit 1) because I redirected stdout to a scratchpad
file whose directory did not exist. No container ran; nothing was written. Rerun as above.

### 3. Refold — prepare, predict (Boltz-2), score

```bash
cd /home/hackathon/amar/data/docs
N=20 TOPN=5 OUT=runs/her2 bash Scripts/run-her2.sh prepare
N=20 TOPN=5 OUT=runs/her2 bash Scripts/run-her2.sh predict
N=20 TOPN=5 OUT=runs/her2 bash Scripts/run-her2.sh score
```

Effective container commands:

```bash
# prepare  (rfd4-proteina her2-novel-2026-11-09)
module-run refold-prepare --input runs/her2/design --output runs/her2/prepared --top-n 5

# predict  (apheris-boltz-v1.3.0-boltz2)
python -m apheris_boltz.cli predict \
  --input runs/her2/prepared/holo/queries.json --output runs/her2/boltz2 \
  --template-paths runs/her2/prepared/holo/template.cif \
  --write-detailed-scores true --num-diffusion-samples 1 --seeds 42 --query-failures skip

# score  (rfd4-proteina her2-novel-2026-11-09)
module-run refold-score --input runs/her2/design --output runs/her2/scored-boltz2 \
  --predictions runs/her2/boltz2 --model boltz2
```

Result: exit 0. 5/5 predicted, 5 rows in `refold_metrics.csv`.

### 4. Post-hoc analysis (not part of the tutorial)

CDRH3 sequences extracted from chain B residues 99-108 of each `*_generated.cif` with a
throwaway script, `/tmp/claude-1001/-home-hackathon-amar/db28f190-a17b-4689-803e-fc6c38dc5455/scratchpad/cdrh3.py`.

## Results

20 designs, 20 distinct CDRH3 loops (parent trastuzumab loop: `WGGDGFYAMD`).

Top 5 refolded with Boltz-2, ranked by `min_ipae` (lower is better):

| design | CDRH3 | min_ipae | iptm | binder_plddt | plddt |
| --- | --- | --- | --- | --- | --- |
| her2_cdrh3_3 | WASSSEYYFD | 0.440 | 0.841 | 97.46 | 95.64 |
| her2_cdrh3_10 | TTPAVAYAHD | 0.445 | 0.805 | 97.85 | 96.18 |
| her2_cdrh3_14 | WHSLHEYYFD | 0.456 | 0.796 | 97.42 | 95.68 |
| her2_cdrh3_12 | DLSLSVYAFD | 0.473 | 0.770 | 97.36 | 95.50 |
| her2_cdrh3_4 | YWSLHYYAMD | 0.488 | 0.734 | 96.02 | 93.73 |

`pass` is empty for every row, as `18-vm.md` says it will be for this example
(`pass_note`: "missing refolded_target_aligned_binder_rmsd.rmsd; fewer than 3 target
atoms fixed in absolute coordinates").

## Deviations from the documented example

- `--num-replicates 20` instead of 100 (requested test batch size).
- `--top-n 5` instead of the script's default 2 and the Settings file's 20 — a readable
  sample of 20 designs. My choice, not from the docs.
- `sequenceTsrK` left unset, per the models table for HER2 novel. `Settings/12-her2.json`
  sets it to 2, but that file configures HER2 *balanced*, not novel.
- Novelty check now run; see the novelty section below. Originally noted as not yet
  known to be novel against the screen.


## Main run — attempt 1 (FAILED): 300 designs, HER2 novel, seed 2 (`runs/her2-main`)

Started 2026-11-14T11:34:28+00:00, finished 2026-11-14T11:34:28+00:00. Exit code: 2.

N=300, TOPN=60, OUT=runs/her2-main, weight `her2-novel-2026-11-09`, seed 2.
`Scripts/run-her2.sh` hardcodes `--seed 1` and has no seed override, so the design
step was run as the container command directly, identical to the script apart from
`--num-replicates 300 --seed 2`. prepare/predict/score use the script unmodified.

### Design

```bash
docker run --rm --gpus all --shm-size 16g --user $(id -u):$(id -g) -e HOME=/tmp \
  -v /etc/passwd:/etc/passwd:ro -v $PWD:/work -w /work \
  quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09 \
  module-run design --input Data/12-her2/her2_cdrh3.yaml --output runs/her2-main/design \
  --num-replicates 300 --seed 2 --recycle-strategy high/front_loaded
```

No designs written. The design module caps `--num-replicates` at 200 per input
(`--num-replicates INTEGER  Designs generated per input, up to 200`), so N=300 is
rejected before any sampling starts. Awaiting a decision on how to get 300.

<details><summary>Design stdout/stderr (tail)</summary>

```
error[usage]: num_replicates must be <= 200
```

</details>

## Novelty check — 20 test designs (`runs/her2`, seed 1)

Dataset located. Both provided tables were already present locally at
`data/hackathon-data/her2-trastuzumab/`, and match the bucket
(`aws s3 ls s3://apheris-oslo-hackathon/data/her2-trastuzumab/`):

| File | Rows |
| --- | --- |
| `tz_her2_facs_3classes_524346.csv` (Chinery et al. FACS) | 524,346 |
| `trastuzumab_her2_spr_all_768.csv` (SPR) | 768 |
| Unique CDRH3 across both | 525,063 (525,042 of length 10) |

`scripts/novelty_check.py` reports, per loop: mutations and identity % vs
trastuzumab `WGGDGFYAMD`, minimum Hamming distance to any provided sequence,
and pass/fail. It evaluates **both** readings of the filter in every run:

- `--rule docs` (default) — as written in `12-her2.md`: identity **<= 70%**
  AND min distance **>= 2**.
- `--rule strict` — identity **< 70%** AND min distance **> 2**.

```bash
python3 scripts/novelty_check.py --loops scripts/designs20.txt --rule docs
python3 scripts/novelty_check.py --loops scripts/designs20.txt --rule strict --csv scripts/novelty20.csv
```

Result: **18/20 pass the documented rule, 16/20 pass the strict rule.**
Full table in `scripts/novelty20.csv`. The two rules disagree on
`her2_cdrh3_2` (identity exactly 70%, dist exactly 2) and `her2_cdrh3_16` (dist exactly 2).
Failing under both: `her2_cdrh3_9` (dist 1) and `her2_cdrh3_13` (identity 80%, dist 1).

Validation: the FACS table's own `edit_distance` column agrees with our Hamming
distance to `WGGDGFYAMD` on 524,346/524,346 rows, and an independent brute-force
scan reproduced the min distances for the borderline designs.

## Main run — attempt 2: 200 designs, HER2 novel, seed 2 (`runs/her2-main`)

Started 2026-11-14T11:46:22+00:00, finished 2026-11-14T12:47:48+00:00. Exit code: 0.

N=200, TOPN=60, OUT=runs/her2-main, weight `her2-novel-2026-11-09`, seed 2.
`Scripts/run-her2.sh` hardcodes `--seed 1` and has no seed override, so the design
step was run as the container command directly, identical to the script apart from
`--num-replicates 200 --seed 2`. prepare/predict/score use the script unmodified.

### Design

```bash
docker run --rm --gpus all --shm-size 16g --user $(id -u):$(id -g) -e HOME=/tmp \
  -v /etc/passwd:/etc/passwd:ro -v $PWD:/work -w /work \
  quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-her2-novel-2026-11-09 \
  module-run design --input Data/12-her2/her2_cdrh3.yaml --output runs/her2-main/design \
  --num-replicates 200 --seed 2 --recycle-strategy high/front_loaded
```

Designs written: 200 `*_generated.cif` under `runs/her2-main/design/generation/inference_samples/`,
200 rows in `generation/inference_metrics/per_example_metrics.csv`.
Resolved config confirms `seed=2, num_replicates=200, recycle_strategy=high/front_loaded`.

<details><summary>Design stdout/stderr (tail)</summary>

```
/usr/local/lib/python3.12/site-packages/torch/jit/_script.py:1491: FutureWarning: `torch.jit.script` is deprecated. Please switch to `torch.compile` or `torch.export`.
  warnings.warn(
cuEquivariance is available and will be used.                                   
Environment variable CCD_MIRROR_PATH not set. Will not be able to use function  
requiring this variable. To set it you may:                                     
  (1) add the line 'export VAR_NAME=path/to/variable' to your .bashrc or .zshrc 
file                                                                            
  (2) set it in your current shell with 'export VAR_NAME=path/to/variable'      
  (3) write it to a .env file in the root of the atomworks.io repository        
Environment variable PDB_MIRROR_PATH not set. Will not be able to use function  
requiring this variable. To set it you may:                                     
  (1) add the line 'export VAR_NAME=path/to/variable' to your .bashrc or .zshrc 
file                                                                            
  (2) set it in your current shell with 'export VAR_NAME=path/to/variable'      
  (3) write it to a .env file in the root of the atomworks.io repository        
generated new fontManager                                                       
[rank: 0] Prevalidating input example: her2_cdrh3                               
[rank: 0]                                                                       
--- condition_spec: her2_cdrh3 ---                                              
Sequence: 326 residues                                                          
                                                                                
Hotspot: 4 atoms                                                                
                                                                                
Coordinate groups: 3 groups                                                     
                                                                                
InternalCoordinate: 2105 atoms                                                  
                                                                                
REMINDER: Conditions not listed above (including coordinates and sequence) are  
NOT provided to the model and will be freely generated.                         
[rank: 0] Dataset inference-dataset: 1 examples (her2_cdrh3_0)                  
```

</details>

### Prepare

```bash
N=200 TOPN=60 OUT=runs/her2-main bash Scripts/run-her2.sh prepare
```

Started 2026-11-14T12:48:16+00:00, finished 2026-11-14T12:48:21+00:00. Exit code: 0.

<details><summary>Output (tail)</summary>

```
```

</details>

### Predict

```bash
N=200 TOPN=60 OUT=runs/her2-main bash Scripts/run-her2.sh predict
```

Started 2026-11-14T12:48:21+00:00, finished 2026-11-14T12:55:43+00:00. Exit code: 0.

<details><summary>Output (tail)</summary>

```
Assigned templates to 240 entities across 60 queries                            
Read 60 queries from /work/runs/her2-main/prepared/holo/queries.json            
Reading config                                                                  
Inferred Boltz v2 from weights file /app/weights/boltz2_conf.ckpt               
Creating inference module                                                       
Seed set to 42                                                                  
Creating trainer                                                                
Using bfloat16 Automatic Mixed Precision (AMP)                                  
GPU available: True (cuda), used: True                                          
TPU available: False, using: 0 TPU cores                                        
💡 Tip: For seamless cloud logging and experiment tracking, try installing      
[litlogger](https://pypi.org/project/litlogger/) to enable LitLogger, which logs
metrics and artifacts automatically to the Lightning Experiments platform.      
Loading weights                                                                 
Seed set to 42                                                                  
Run prediction                                                                  
Seed set to 42                                                                  
💡 Tip: For seamless cloud uploads and versioning, try installing               
[litmodels](https://pypi.org/project/litmodels/) to enable LitModelCheckpoint,  
which syncs automatically with the Lightning model registry.                    
LOCAL_RANK: 0 - CUDA_VISIBLE_DEVICES: [0]                                       
/app/.pixi/envs/boltz-image/lib/python3.11/site-packages/pytorch_lightning/utilities/_pytree.py:21: `isinstance(treespec, LeafSpec)` is deprecated, use `isinstance(treespec, TreeSpec) and treespec.is_leaf()` instead.
/app/.pixi/envs/boltz-image/lib/python3.11/site-packages/pytorch_lightning/trainer/connectors/data_connector.py:434: The 'predict_dataloader' does not have many workers which may be a bottleneck. Consider increasing the value of the `num_workers` argument` to `num_workers=15` in the `DataLoader` to improve performance.
/app/.pixi/envs/boltz-image/lib/python3.11/site-packages/pytorch_lightning/utilities/data.py:123: Your `IterableDataset` has `__len__` defined. In combination with multi-process data loading (when num_workers > 1), `__len__` could be inaccurate if each worker is not configured independently to avoid having duplicate data.
Read 4 entities from `_entity` table.                                           
Read 4 polymer types from `_entity_poly` table.                                 
Read 4 entities with a total of 326 residues from `_entity_poly_seq` block.     
Read 4 chains from `_struct_asym` table.                                        
Read 2509 atoms from `_atom_site` table.                                        
Identified 326 resolved residue(s) from `_atom_site` table.                     
Parsed molecule with 2513 atoms (2509 resolved), 2578 bonds, 326 residues, and 4
chains.                                                                         
/app/.pixi/envs/boltz-image/lib/python3.11/site-packages/pytorch_lightning/loops/prediction_loop.py:307: Couldn't infer the batch indices fetched from your dataloader: `DataLoader`
/app/.pixi/envs/boltz-image/lib/python3.11/site-packages/torch/jit/_script.py:1488: DeprecationWarning: `torch.jit.script` is deprecated. Please switch to `torch.compile` or `torch.export`.
  warnings.warn(
/app/.pixi/envs/boltz-image/lib/python3.11/site-packages/cuequivariance_ops_torch/triangle_attention.py:623: UserWarning: Generic kernel expects bias to be float32 so it's going to be cast to torch.float32. Check if you can change your code for maximum performance.
  warn_eager_only(
/app/.pixi/envs/boltz-image/lib/python3.11/site-packages/boltz-community/src/boltz/data/feature/featurizerv2.py:1742: DeprecationWarning: __array_wrap__ must accept context and return_scalar arguments (positionally) in the future. (Deprecated NumPy 2.0)
  deletion = np.pi / 2 * np.arctan(deletion / 3)
Predicting ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━ 60/60 0:06:49 • 0:00:00 0.15it/s 
```

</details>

### Score

```bash
N=200 TOPN=60 OUT=runs/her2-main bash Scripts/run-her2.sh score
```

Started 2026-11-14T12:55:43+00:00, finished 2026-11-14T12:55:51+00:00. Exit code: 0.

<details><summary>Output (tail)</summary>

```
generation/inference_samples/her2_cdrh3_170/her2_cdrh3_170_generated.cif  her2_cdrh3_170  boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_171/her2_cdrh3_171_generated.cif  her2_cdrh3_171  boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_175/her2_cdrh3_175_generated.cif  her2_cdrh3_175  boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_179/her2_cdrh3_179_generated.cif  her2_cdrh3_179  boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_18/her2_cdrh3_18_generated.cif    her2_cdrh3_18   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_181/her2_cdrh3_181_generated.cif  her2_cdrh3_181  boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_184/her2_cdrh3_184_generated.cif  her2_cdrh3_184  boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_192/her2_cdrh3_192_generated.cif  her2_cdrh3_192  boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_193/her2_cdrh3_193_generated.cif  her2_cdrh3_193  boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_195/her2_cdrh3_195_generated.cif  her2_cdrh3_195  boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_197/her2_cdrh3_197_generated.cif  her2_cdrh3_197  boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_22/her2_cdrh3_22_generated.cif    her2_cdrh3_22   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_25/her2_cdrh3_25_generated.cif    her2_cdrh3_25   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_27/her2_cdrh3_27_generated.cif    her2_cdrh3_27   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_28/her2_cdrh3_28_generated.cif    her2_cdrh3_28   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_30/her2_cdrh3_30_generated.cif    her2_cdrh3_30   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_31/her2_cdrh3_31_generated.cif    her2_cdrh3_31   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_35/her2_cdrh3_35_generated.cif    her2_cdrh3_35   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_37/her2_cdrh3_37_generated.cif    her2_cdrh3_37   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_39/her2_cdrh3_39_generated.cif    her2_cdrh3_39   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_43/her2_cdrh3_43_generated.cif    her2_cdrh3_43   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_44/her2_cdrh3_44_generated.cif    her2_cdrh3_44   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_5/her2_cdrh3_5_generated.cif      her2_cdrh3_5    boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_57/her2_cdrh3_57_generated.cif    her2_cdrh3_57   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_60/her2_cdrh3_60_generated.cif    her2_cdrh3_60   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_61/her2_cdrh3_61_generated.cif    her2_cdrh3_61   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_64/her2_cdrh3_64_generated.cif    her2_cdrh3_64   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_65/her2_cdrh3_65_generated.cif    her2_cdrh3_65   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_66/her2_cdrh3_66_generated.cif    her2_cdrh3_66   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_7/her2_cdrh3_7_generated.cif      her2_cdrh3_7    boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_70/her2_cdrh3_70_generated.cif    her2_cdrh3_70   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_71/her2_cdrh3_71_generated.cif    her2_cdrh3_71   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_80/her2_cdrh3_80_generated.cif    her2_cdrh3_80   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_82/her2_cdrh3_82_generated.cif    her2_cdrh3_82   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_85/her2_cdrh3_85_generated.cif    her2_cdrh3_85   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_9/her2_cdrh3_9_generated.cif      her2_cdrh3_9    boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_91/her2_cdrh3_91_generated.cif    her2_cdrh3_91   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_94/her2_cdrh3_94_generated.cif    her2_cdrh3_94   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_95/her2_cdrh3_95_generated.cif    her2_cdrh3_95   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
generation/inference_samples/her2_cdrh3_98/her2_cdrh3_98_generated.cif    her2_cdrh3_98   boltz2  predicted  42    1       B             A;C;D          antibody         "missing refolded_target_aligne
```

</details>

### Results — main run

60/60 refolded and scored (`status=predicted`), 60 rows in `scored-boltz2/refold_metrics.csv`.
As in the 20-design run, `pass` is empty for every row (`pass_note`: missing
`refolded_target_aligned_binder_rmsd.rmsd`), so ranking is by `min_ipae`.

`min_ipae`: best 0.308, median 0.380, worst 0.480 — better across the board than the
20-design seed-1 run (best 0.440).

**Duplicate loops.** The 200 designs contain only **164 unique CDRH3s**; 21 loops recur,
one up to 5 times. The 60 refolded are only **52 unique loops**. refold-prepare ranks by
score and does not deduplicate, so duplicate sequences consume refold slots and appear
repeatedly in the rankings with identical metrics.

Top unique loops by `min_ipae`, joined with the novelty filter:

| design | CDRH3 | min_ipae | iptm | binder_plddt | identity | min_dist | docs | strict |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| her2_cdrh3_140 | SLSLSNYAFD | 0.308 | 0.926 | 98.59 | 30% | 3 | pass | pass |
| her2_cdrh3_11 | WGSPGFYAMD | 0.309 | 0.929 | 98.33 | 80% | 1 | FAIL | FAIL |
| her2_cdrh3_70 | WGPPGFYAFD | 0.314 | 0.923 | 97.78 | 70% | 0 | FAIL | FAIL |
| her2_cdrh3_98 | TRAAVGYAHD | 0.323 | 0.903 | 98.32 | 30% | 4 | pass | pass |
| her2_cdrh3_31 | WSSAVGYAHD | 0.324 | 0.926 | 98.16 | 40% | 3 | pass | pass |
| her2_cdrh3_39 | SLYVNNYSYD | 0.329 | 0.902 | 98.13 | 20% | 4 | pass | pass |
| her2_cdrh3_1 | SLSSSVYAFD | 0.336 | 0.860 | 98.01 | 30% | 2 | pass | FAIL |
| her2_cdrh3_130 | WHSLSEYYFD | 0.340 | 0.905 | 97.40 | 30% | 2 | pass | FAIL |
| her2_cdrh3_14 | WPDGVYYAHD | 0.340 | 0.917 | 97.07 | 40% | 2 | pass | FAIL |
| her2_cdrh3_27 | DLFFPYYAFD | 0.345 | 0.892 | 95.99 | 30% | 3 | pass | pass |

`her2_cdrh3_70` (`WGPPGFYAFD`, min_dist 0) is **already in the provided screen data** —
not a novel design at all, despite ranking 3rd on interface confidence.

### Novelty over all 200 designs

```bash
python3 scripts/novelty_check.py --loops designs200.txt --rule docs --csv novelty200.csv
```

**172/200 pass the documented rule; 142/200 pass the strict rule.**
Of the 52 unique refolded loops: 47 pass docs, 35 pass strict.
CDRH3s extracted from chain B residues 99-108 of each `*_generated.cif`.

## Refold rerun — 2026-11-14T12:56:38..13:04:14 (redundant, reproducibility check)

The prepare/predict/score pipeline above had already been run (outputs on disk at
12:55:51). I relaunched the identical pipeline at 12:56:38 without checking for
existing outputs. `Scripts/run-her2.sh` begins each step with `rm -rf`, so the
original `prepared/` and `boltz2/` were deleted and regenerated; `scored-boltz2/`
was copied aside first and then overwritten.

Outcome: the rerun reproduces the earlier scores exactly.
`runs/her2-main/scored-boltz2/refold_metrics.csv` is byte-identical to the
preserved copy (md5 `5d1677de1a1762bd6c4ff6a885ffa4b0`, 60 rows + header) —
same image, same `--seeds 42`, same inputs.

```bash
TOPN=60 OUT=runs/her2-main bash Scripts/run-her2.sh prepare   # 12:56:38 -> 12:56:43
TOPN=60 OUT=runs/her2-main bash Scripts/run-her2.sh predict   # 12:56:43 -> 13:04:06, 60/60
TOPN=60 OUT=runs/her2-main bash Scripts/run-her2.sh score     # 13:04:06 -> 13:04:14
```

## CDRH3 extraction — all 200 main-run designs

`scripts/extract_cdrh3.py` (new) reads the `_atom_site` loop of each
`*_generated.cif` and emits chain B (auth), author residues 99-108, as
`<loop>  <design_id>` — the input format `novelty_check.py --loops` takes.

Validation: re-deriving the 20 pilot loops from `runs/her2/design` reproduces
`scripts/designs20.txt` exactly.

```bash
python3 -I scripts/extract_cdrh3.py data/docs/runs/her2-main/design > scripts/designs200.txt
```

200 loops written, **164 distinct** — the sampler repeats some loops across
replicates (e.g. `her2_cdrh3_7` and `her2_cdrh3_140` are both `SLSLSNYAFD`,
and refold to identical metrics).

## Novelty check — 200 main-run designs (`runs/her2-main`, seed 2)

```bash
python3 -I scripts/novelty_check.py --loops scripts/designs200.txt --rule docs   --csv scripts/novelty200_docs.csv
python3 -I scripts/novelty_check.py --loops scripts/designs200.txt --rule strict --csv scripts/novelty200_strict.csv
```

Result: **172/200 pass the documented rule (86%), 142/200 pass the strict rule (71%).**
Each pass took ~33 s against the 525,042 length-10 provided sequences.
Full tables in `scripts/novelty200_docs.csv` and `scripts/novelty200_strict.csv`
(both carry `pass_docs` and `pass_strict` columns; the files differ only in the
`verdict` column).

Restricted to the 60 designs that were refolded: **52 pass docs, 37 pass strict.**

### Top 10 refolded by `min_ipae`, with novelty

| design | CDRH3 | min_ipae | iptm | binder_plddt | min_dist | identity % | docs/strict |
| --- | --- | --- | --- | --- | --- | --- | --- |
| her2_cdrh3_140 | SLSLSNYAFD | 0.308 | 0.926 | 98.59 | 3 | 30 | pass/pass |
| her2_cdrh3_7 | SLSLSNYAFD | 0.308 | 0.926 | 98.59 | 3 | 30 | pass/pass |
| her2_cdrh3_11 | WGSPGFYAMD | 0.309 | 0.929 | 98.33 | 1 | 80 | fail/fail |
| her2_cdrh3_28 | WGSPGFYAMD | 0.309 | 0.929 | 98.33 | 1 | 80 | fail/fail |
| her2_cdrh3_70 | WGPPGFYAFD | 0.314 | 0.923 | 97.78 | 0 | 70 | fail/fail |
| her2_cdrh3_98 | TRAAVGYAHD | 0.323 | 0.903 | 98.32 | 4 | 30 | pass/pass |
| her2_cdrh3_31 | WSSAVGYAHD | 0.324 | 0.926 | 98.16 | 3 | 40 | pass/pass |
| her2_cdrh3_39 | SLYVNNYSYD | 0.329 | 0.902 | 98.13 | 4 | 20 | pass/pass |
| her2_cdrh3_1 | SLSSSVYAFD | 0.336 | 0.860 | 98.01 | 2 | 30 | pass/fail |
| her2_cdrh3_130 | WHSLSEYYFD | 0.340 | 0.905 | 97.40 | 2 | 30 | pass/fail |

Note `her2_cdrh3_70` (`WGPPGFYAFD`, min_dist 0) is already present in the provided
screen data, and the two 80%-identity hits are near-parent loops — the structurally
best-scoring designs skew toward the trastuzumab loop, so novelty has to be applied
as a filter rather than read off the ranking. The best design passing both readings
is **`her2_cdrh3_140` / `her2_cdrh3_7` — `SLSLSNYAFD`, min_ipae 0.308, iptm 0.926,
3 mutations from anything in the screen, 30% identity to trastuzumab.**

`pass` is empty for all 60 rows, as for the pilot (`pass_note`: missing
`refolded_target_aligned_binder_rmsd.rmsd`).

## Stability proxy — 37 strict-passing designs (`scripts/stability_proxy.csv`)

Loop-level glycine count, proline count, longest G/S run, plus whole-binder refold
pLDDT, each against the parent `WGGDGFYAMD` (3 Gly, 0 Pro, G/S run 2). **Proxy only,
not a Tm prediction; uncalibrated.** No design has more glycines than the parent
(33 have none, 4 have one). 4 flagged worse on flexible-run length only
(`DLSSSDHAFD`, `DLSSSDFAFD`, `DLSSSDYAFD`, `YASSSDYYFD`: `SSS` run of 3 vs 2). 3 carry
a proline (reported, not scored). The four picks are not worse than the parent.
The parent was not refolded, so pLDDT is compared with the design median (97.79)
instead; binder pLDDT is global, not loop-local, and spans only 95.3-98.6.
Script: `scripts/stability_proxy.py`.

## De novo VHH vs HER2 domain IV (`runs/denovo`, antibody-2026-10-30, seed 4)

Goal: a new VHH to HER2 domain IV (not a CDRH3 redesign of trastuzumab), same epitope
as trastuzumab. Weight/image: `...apheris-rfd4-proteina-v0.3.0-antibody-2026-10-30`.

### Input format (from `15-antibody-design.md`, `03-antibodies.md`, and the shipped example)

The shipped VHH example (`Data/03-vhh/pmvk_caplacizumab.cif.zst`) carries the whole
design spec inside the mmCIF, with no separate YAML. Chain H is an 86-residue
caplacizumab framework in which each CDR is one `HETATM UNK` placeholder with NaN
coordinates; its `comp_id` is the allowed length range (H23 `13-13`, H38 `9-10`,
H76 `12-21`). Per-atom columns carry `condition_sequence`, `condition_hotspot`,
`condition_cath-label` (`2.60.40.x` on the antibody) and `condition_CoordinateGroup`.
Only backbone atoms (N, CA, C, O) are in a coordinate group; side chains are -1, and
hotspots are side-chain atoms. The example's antigen (PMVK) was not usable for HER2,
so a hybrid had to be built: `scripts/build_denovo_input.py` grafts caplacizumab
chain H (3 placeholders) with HER2 chains C (489-580, 92 residues) and D (591-607,
17 residues) from `her2_fv_1n8z.cif`, with the same 4 hotspots as the CDRH3 run
(`D/593/NZ`, `C/573/CD1`, `D/591/CG2`, `C/572/CD`). 1,484 atoms.

### Input build: four builds, three failures, and the coordinate-group finding

Validation = a 1-2 design test run (`module-run design --num-replicates 2`).

| # | build | result |
| --- | --- | --- |
| 1 | HER2 on coordinate group 0, all non-hotspot atoms incl. side chains; HER2 rows first | FAIL 13:41 |
| 2 | HER2 on group 0, backbone atoms only (as in the example) | FAIL 13:48 |
| 3 | as 2, with antibody rows first (as in the example) | FAIL 13:49 |
| 4 | **relative groups: framework = 1, HER2 backbone = 2, no group 0** | **PASS 13:50** |

Builds 1-3 all failed with the same message:
`Transforms failed at stage MaterializeInternalCoords: Absolute coordinate targets
must exist exactly on coordinate group 0.` Neither the side-chain/backbone split nor
row order was the cause.

**Finding:** with a two-chain HER2 antigen, pinning HER2 as absolute coordinates on
group 0 (the PMVK example's convention) is rejected. The convention that validates is
the one the existing `her2_cdrh3.yaml` run uses: two relative internal-geometry groups
(antibody framework 1, HER2 backbone 2) and no absolute group. Cause not diagnosed
beyond that. **Consequence:** the VHH-HER2 pose is not pinned, so hotspots are the only
steer toward the epitope. `scripts/build_denovo_input.py --antigen-mode relative`
builds the working input; `absolute` reproduces the failing one. The failed input is
kept as `runs/denovo/input_absolute_FAILED/`. (`runs/` is gitignored; the input is
regenerable from the script.) The 2-design test gave two real VHHs (123 and 118
residues) within 3.3-3.7 A of the hotspot atoms.

### Run

```bash
docker run --rm --gpus all --shm-size 16g --user $(id -u):$(id -g) -e HOME=/tmp \
  -v /etc/passwd:/etc/passwd:ro -v $PWD:/work -w /work \
  quay.io/apheris/oslo-hackathon:apheris-rfd4-proteina-v0.3.0-antibody-2026-10-30 \
  module-run design --input runs/denovo/input --output runs/denovo/design \
  --num-replicates 150 --seed 4 --recycle-strategy high/front_loaded
IMG=antibody-2026-10-30 TOPN=30 OUT=runs/denovo bash Scripts/run-her2.sh prepare   # unmodified script
IMG=antibody-2026-10-30 TOPN=30 OUT=runs/denovo bash Scripts/run-her2.sh predict
IMG=antibody-2026-10-30 TOPN=30 OUT=runs/denovo bash Scripts/run-her2.sh score
```

Timeline (2026-11-14): design 13:51:36 -> 14:18:59 (150/150, exit 0, ~27 min, 8.6-13
s/design); prepare 14:18:59 -> 14:19:01; predict 14:19:01 -> 14:23:05 (30/30); score
-> 14:23:08. Design step wrapped in `timeout` to a 15:30 hard stop; not needed.
Variable CDR lengths did NOT break the shared-template check in refold-prepare.

### Refold results (30 designs, Boltz-2, `min_ipae` in A)

min_ipae min 1.64 / median 12.4 / max 21.6; iptm 0.35 / 0.45 / 0.89; binder pLDDT
85.7 / 91.7 / 97.6. iptm > 0.7: 4/30; > 0.8: 2/30. `min_ipae` < 2.0 A (antibody pass
line): **2/30** (`rep117` 1.638, `rep106` 1.701), then `rep48` 3.39, `rep33` 5.54,
`rep136` 5.97. `pass` is empty/false for all 30; most fail `oracle_ppi.min_ipae_coi`,
and the top two have empty `pass` (missing `refolded_target_aligned_binder_rmsd`).
These are not comparable with the CDRH3-redesign run (0.31-0.48 A): that was one loop
inside a known complex, this is de novo with an unpinned pose. Only the top 30 of 150
were refolded, so the pass rate of the other 120 is unknown.

### Check 1 — epitope (`scripts/epitope_check.py`, `scripts/epitope_check_top5.csv`)

For the top 5 by min_ipae, in the refolded Boltz-2 poses: minimum distance from any
VHH heavy atom to each hotspot atom, the VHH footprint on HER2 (residues within
4.5 A), and overlap with trastuzumab's own epitope (19 HER2 residues within 4.5 A of
the Fv in 1N8Z). HER2 residue mapping (Boltz renumbers 1..N) uses the query's
`residue_ids` and asserts residue names match 1N8Z; no mismatch occurred.

| design | min_ipae | iptm | D/593/NZ | C/573/CD1 | D/591/CG2 | C/572/CD | hotspots <=4.5 A | footprint (res) | in Tz epitope | % of footprint in Tz | % of Tz covered |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| rep117 | 1.64 | 0.89 | 3.4 | 4.6 | 4.2 | 3.8 | 3/4 | 14 | 11 | 79% | 58% |
| rep106 | 1.70 | 0.89 | 4.5 | 3.8 | 5.0 | 3.7 | 3/4 | 15 | 11 | 73% | 58% |
| rep48 | 3.39 | 0.78 | 4.0 | 4.9 | 6.2 | 3.3 | 2/4 | 9 | 9 | 100% | 47% |
| rep33 | 5.54 | 0.71 | 8.0 | 7.6 | 15.9 | 5.3 | 0/4 | 17 | 9 | 53% | 47% |
| rep136 | 5.97 | 0.66 | 8.1 | 7.8 | 11.6 | 3.4 | 1/4 | 21 | 8 | 38% | 42% |

The top three sit on the trastuzumab epitope; `rep33` and `rep136` have drifted partly
off it (0-1 hotspots in contact, 38-53% of the footprint inside the epitope). Caveats:
single Boltz-2 sample per design; the input is only a domain-IV fragment (HER2
489-607), so "elsewhere" can only mean elsewhere on that fragment; being on the
epitope in a prediction is not evidence of binding.

### Check 2 — novelty vs SAbDab2 ML subset (`scripts/sabdab_novelty.py`, `scripts/novelty_sabdab2.csv`)

**Reference set is a curated subset: the SAbDab2 ML dataset (Zenodo 20083995, CC-BY-4.0,
md5 `0dbb4cc499e9eb77f14008b232f2c38c`), 15,641 structures / 8,641 PDB entries / 5,301
unique antibodies. It is NOT full SAbDab2 (21,237 distinct instances), and not OAS or
any other collection.** (An earlier attempt to download "SAbDab" fetched a 2 KB HTML
page from a stale URL; I had not opened it before treating it as a download.)

Identity on CDRs, not the full sequence: the 83 fixed caplacizumab framework residues
are 65-71% of each design, so full-length identity to caplacizumab-like VHHs is ~65-70%
by construction. Identity = LCS / max(len), which can only over-report similarity.
CDR3 = residues between Cys104 and the J-region Trp (validated against the database:
`C`+CDRH3+`W` in 14,831/15,542 = 95.4%). CDR1 = placeholder[3:11] (validated 11,856/
14,011 = 84.6%). CDR2 = placeholder[1:-1] (inferred; FR3 spacing 37 in 15,129 of the
database, as predicted). So the CDR3 result is firmer than the CDR1+2+3 result.

Fraction of the 30 with max identity < 70%: **CDR3 27/30; CDR1+2+3 29/30.**
Failing CDR3: `rep133` 0.714 (FV 8bw0), `rep15` 0.769 (FAB 9g6s), `rep127` 0.706 (SD-H
6b73); failing CDR1+2+3: `rep15` 0.724. Top 5 by min_ipae all pass (CDR3 0.62-0.67).

**The threshold is weak here.** Control: shuffling the CDR3 middle (ends kept) and
re-searching gives median max identity 0.61, the same as the designs' median 0.61.
So the designs' nearest neighbours are at chance level: "no design is a near copy of a
known CDR3 (none above 0.77)" is supported, "the designs are distinct from known
antibodies in a way a random loop would not be" is not. Nothing here establishes
novelty against antibodies outside the subset.

### Not done / open

- Epitope check covered the top 5 only. Pose-based evidence for the other 25 is absent.
- 120 of 150 designs were never refolded.
- `pass` unavailable for all 30 (missing rmsd evidence for the top two).
- Boltz-2 iptm here is not comparable with the Chai-1 calibration in
  `scripts/calibration.csv`.
- `git push` failed earlier (`could not read Username for 'https://github.com'`): no
  credentials on this VM, so the first commit (`5a43142`) exists locally only.
