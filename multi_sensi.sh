#!/bin/bash
# multi_sensi.sh
source /usr/local/miniconda3/etc/profile.d/conda.sh
cd darwen
conda activate py310
python multi_sensi.py >> output_sensi.txt