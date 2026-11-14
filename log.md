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
