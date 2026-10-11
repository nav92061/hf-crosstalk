#!/bin/bash
# Usage: bash run_lctsc.sh <limit|0> [prev_csv]
set -u
mkdir -p out
nvidia-smi --query-gpu=name,driver_version,memory.total --format=csv | tee out/gpu_info.txt
python -c "import torch, importlib.metadata as m; print('torch', torch.__version__, 'cuda_available', torch.cuda.is_available(), 'TotalSegmentator', m.version('TotalSegmentator'), 'nnunetv2', m.version('nnunetv2'))" | tee -a out/gpu_info.txt
python -c "import torch, sys; sys.exit(0 if torch.cuda.is_available() else 3)" || { echo "CUDA not available to torch - aborting"; exit 3; }
pip freeze > out/pip_freeze.txt
EXTRA=""
if [ -n "${2:-}" ]; then EXTRA="--done $2"; fi
python lctsc_external_accuracy_job.py lctsc_case_list.csv out/external_accuracy_cases_raw.csv --limit "$1" $EXTRA
