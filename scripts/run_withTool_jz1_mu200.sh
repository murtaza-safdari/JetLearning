#!/bin/bash
python JetLearning/scripts/Run_withTool.py \
    --submitDir "JetLearningOut/jz1_mu200_test" \
    -w \
    --inputFileList \
    --inputFiles "/u/at/acukierm/nfs/Voronoi_xAOD/mc15_14TeV.147911.JZ1W.e2403_s2638_s2206_r7701.filelist.txt" \
    --nevents 10 \
    --files_per_job 400 \
    --driver "direct"
exit 0
