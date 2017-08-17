#!/bin/bash
python JetLearning/scripts/Run_withTool.py \
    --submitDir "JetLearningOut/HL_JZ1_EM_PileupJets" \
    -w \
    --doClusterInfo \
    --doTracks \
    --doJetReclustering \
    --inputFileList \
    --inputFiles user.cdelitzs.147911.Pythia8_AU2CT10_jetjet_JZ1W.recon.DAOD_Trimming.e2403_s3142_s3143_r9589_v0_StreamAOD.filelist.txt \
    --files_per_job 1 \
    --driver "lsf"
exit 0


# DQ2:
#    --inputDQ2 \
#    --inputFiles "mc15_13TeV.361021.Pythia8EvtGen_A14NNPDF23LO_jetjet_JZ1W.merge.AOD.e3569_s2576_s2132_r6765_r6282/" \
# Local:
#    --inputFileList \
#    --inputFiles mc15_13TeV.361021.JZ1W.e3569_s2576_s2132_r6765_r6282.filelist.txt \
# No PU:
#    --inputFileList \
#    --inputFiles "/u/at/acukierm/nfs/Voronoi_xAOD/mc15_13TeV.361021.JZ1W.e3569_s2576_s2132_r8256.filelist.txt" \
