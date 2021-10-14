# BIPSPI v2: Enhancing partner-specific binding site prediction with better data

[BIPSPI](http://bipspi.cnb.csic.es/xgbPredApp/) (xgBoost Interface Prediction of Specific Partner Interactions) is a partner specific 
binding site predictor that employs the XGBoost algorithm. BIPSPI v2 has been trained on
newer and larger datasets compiled from the PDB, severely improving its performance.
<br>BIPSPI v2 can be employed to predict partner-specific binding sites given two atomic models, 
two sequences, or one sequence and one structure.

BIPSPI v2 is distributed as a Docker image and as a GitHub repository. Installation is only required in
the latter case. Complete guide can be found in:

- [DOCKER image](docs/docker_help.md)
- [GitHub repository](docs/repo_help.md)


### New features


- Homo-complexes and hetero-complex specific models, boosting performance especially when predicting homo-complexes.
- Sequence vs structure mode. Version 1 could only be executed with two sequences or two structures. Now one sequence and one structure could be used as input as well. 
- Atomatic correlated mutations (optinal) 



