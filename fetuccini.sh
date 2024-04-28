#!/bin/bash
# darwen_run.sh
source /usr/local/miniconda3/etc/profile.d/conda.sh
cd darwen
conda activate py310
python Fetuccini.py >> output.txt