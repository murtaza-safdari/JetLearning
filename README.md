# JetLearning
Simple analysis for writing out ntuples for Jet-related performance analyses, including Machine Learning.

## Dependencies
This package makes use of [UChicago](https://github.com/UCATLAS)'s [xAODAnaHelpers](https://github.com/UCATLAS/xAODAnaHelpers) package.

## Installing
To install,
```bash
mkdir myRootCore && cd $_
setupATLAS
rcSetup Base,2.4.33
git clone https://github.com/UCATLAS/xAODAnaHelpers
git clone https://github.com/AvivCukierman/JetLearning
mv Voronoi_xAOD MyAnalysis
rc checkout_pkg atlasoff/AsgExternal/Asg_FastJet/tags
rc checkout_pkg atlasoff/AsgExternal/Asg_FastJetContrib/tags
rc find_packages
rc compile
```

## Running
Two sample run scripts (`run_withTool_test.sh`, for running on a file stored locally) and (`run_withTool_jz1_lsf.sh`, for running on a file stored on the grid; uses DQ2 which is currently deprecated) are included in the scripts folder. The `--driver` option allows you to run on batch, by specifying `--driver lsf`. I'll document this better later but just ask me for help if you need it for now.
