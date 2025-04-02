#!/bin/bash

python crd_convert.py
python trj_analysis.py
python wham_analysis.py
python mbar_analysis_umbrella_1d.py
python mbar_analysis_umbrella_block.py
python avecrd_analysis.py
python kmeans_clustering.py
python hb_analysis_count_atom.py
python hb_analysis_count_snap.py
python rmsd_analysis.py
python drms_analysis.py
python rg_analysis.py
python msd_analysis.py
python diffusion_analysis.py
