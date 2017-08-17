#!/bin/bash
python JetLearning/scripts/Run_withTool.py \
    --submitDir "JetLearningOut/jz1_mu200_test" \
    -w \
    --doClusterInfo \
    --doTracks \
    --doJetReclustering \
    --inputFileList \
    --inputFiles user.cdelitzs.147911.Pythia8_AU2CT10_jetjet_JZ1W.recon.DAOD_Trimming.e2403_s3142_s3143_r9589_v0_StreamAOD.filelist.txt \
    --nevents 5 \
    --files_per_job 400 \
    --driver "direct"
exit 0
