#!/bin/bash
python JetLearning/scripts/Run_withTool.py \
    --submitDir "JetLearningOut/JZ2_LC_ClusterInfo" \
    -w \
    --doLC \
    --inputFileList \
    --inputFiles mc15_13TeV.361022.JZ2W.e3668_s2576_s2132_r6765_r6282.filelist.txt \
    --files_per_job 1 \
    --driver "lsf"
exit 0


