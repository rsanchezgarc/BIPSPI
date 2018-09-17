################################################################################
# BIPSPI: xgBoost Interface Prediction of Specific Partner Interactions        #
################################################################################

ACADEMIC USE ONLY. This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY

xgBoost based Interface Prediction of Specific Partner Interactions (BIPSPI) is a new method for 
the prediction of partner-specific protein interfaces from pdb files or input sequences. 
BIPSPI employs Extreme Gradient Boosting (XGBoost) models trained on the residue pairs of the protein complexes 
compiled in Protein-Protein Docking Benchmark version 5 and an scoring function that converts pair prediction 
to interface residue predictions. contact: rsanchez@cnb.csic.es; jsegura@cnb.csic.es

CONTENT:
  1) Installation
  2) Use
    2.1) Train model
    3.2) Predict
    
-------------------------
- 1. Installation       -
-------------------------

BIPSPI make use of several bioinformatics tool that are distributed within its docker. No need for installation
if this docker is used. You only have to compile an uniref90 sequence database for psiblast and, optionally,
a uniclust30 database for hhblits if you want to use correlated mutations. Path to these databases must be
indicated in ./configFiles/configFile.cfg


By using BIPSPI you are accepting the Terms and Conditions of the licenses of the following packages:

- PSAIA 1.0 (http://bioinfo.zesoi.fer.hr/index.php/en/10-category-en-gb/tools-en/19-psaia-en)
- DSSP  (https://swift.cmbi.umcn.nl/gv/dssp/index.html)
- AL2CO  (http://prodata.swmed.edu/al2co/al2co.php)
      Warning. Default code has too small buffers for input names, code was modified from char[500]  to char[1024] and compiled
    AL2CO dependencies:
      - cd-hit (http://weizhongli-lab.org/cd-hit/)
      - clustalw (http://www.clustal.org/)
- qhull (http://www.qhull.org/)
- psiblast 2.2.31+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
- SPIDER2 (http://sparks-lab.org/yueyang/server/SPIDER2/)
- hhblits (OPTIONAL, needed if correlated mutations want to be used)
          (https://github.com/soedinglab/hh-suite)
- ccmpred (OPTIONAL, needed if correlated mutations want to be used)
          (https://github.com/soedinglab/CCMpred)
- Anaconda 5.0.1 (https://anaconda.org/)

- Python packages (as reported by Anaconda):
    name: xgbpred
    channels:
      - bioconda
      - conda-forge
      - anaconda
      - defaults
    dependencies:
      - enum34=1.1.6=py27h99a27e9_1
      - freetype=2.8=hab7d2ae_1
      - funcsigs=1.0.2=py27h83f16ab_0
      - joblib=0.11=py27_0
      - jpeg=9b=h024ee3a_2
      - libgcc-ng=7.2.0=hdf63c60_3
      - libgfortran=3.0.0=1
      - libpng=1.6.34=hb9fc6fc_0
      - libstdcxx-ng=7.2.0=hdf63c60_3
      - libtiff=4.0.9=he85c1e1_1
      - llvmlite=0.19.0=py27_0
      - msgpack-python=0.5.6=py27h6bb024c_0
      - numba=0.34.0=np112py27_0
      - olefile=0.45.1=py27_0
      - openblas=0.2.19=0
      - pillow=4.2.1=py27h7cd2321_0
      - pip=9.0.1=py27_1
      - python=2.7.13=0
      - python-dateutil=2.7.3=py27_0
      - pytz=2018.4=py27_0
      - readline=6.2=2
      - reportlab=3.4.0=py27_0
      - setuptools=39.1.0=py27_0
      - simplejson=3.11.1=py27_0
      - singledispatch=3.4.0.3=py27h9bcb476_0
      - six=1.10.0=py27_0
      - sqlite=3.13.0=0
      - tk=8.5.18=0
      - wget=1.18=0
      - wheel=0.31.1=py27_0
      - xz=5.2.4=h14c3975_4
      - zlib=1.2.11=ha838bed_2
      - biopython=1.70=np112py27_0
      - mmtf-python=1.0.2=py27_0
      - blas=1.1=openblas
      - numpy=1.12.1=py27_blas_openblas_200
      - pandas=0.21.0=py27_0
      - scikit-learn=0.19.1=py27_blas_openblas_200
      - scipy=0.19.1=py27_blas_openblas_202
      - xgboost=0.6a2=py27_2
      - asn1crypto=0.24.0=py27_0
      - ca-certificates=2018.03.07=0
      - certifi=2018.4.16=py27_0
      - cffi=1.11.5=py27h9745a5d_0
      - chardet=3.0.4=py27hfa10054_1
      - cryptography=2.2.2=py27h14c3975_0
      - idna=2.6=py27h5722d68_1
      - ipaddress=1.0.22=py27_0
      - libffi=3.2.1=hd88cf55_4
      - openssl=1.0.2o=h20670df_0
      - pycparser=2.18=py27hefa08c5_1
      - pyopenssl=18.0.0=py27_0
      - pysocks=1.6.8=py27_0
      - requests=2.18.4=py27hc5b0589_1
      - urllib3=1.22=py27ha55213b_0
      - pip:
        - bz2file==0.98
        - gputil==1.3.0
    prefix: /services/xgbpred/app/miniconda2/envs/xgbpred


Otherwise, you should install all them manually and edit
./configFiles/configFile.cfg file consequently to point to installation location
  
  
To compile uniref90 blastDb you can use the following comands
  wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
  gunzip -v uniref90.fasta.gz
  makeblastdb -in uniref90.fasta -dbtype prot -out uniref90.fasta -hash_index
  
-------------------------
- 2. Use                -
-------------------------

*** 2.1 Train Model ***
  
In order to train a model you need a set of protein complexes with the format of Docking Benchmark v5 stored in a 
directory. For each complex, 4 pdb files must be provided, 2 for ligand (bound and unbound state) and
2 for receptor (bound and unbound). If just bound pdb files available, you must symlink them in order to
have four different files. filenames follow this setting: prefix_X_Y.pdb, where prefix is an id for the complex (a pdb id or any other
unique string), X is l or r (ligand or receptor) and Y is u or b (bound or unbound).

For example:
~/path/to/trainPdbs/
    1A2K_l_b.pdb
    1A2K_l_u.pdb
    1A2K_r_b.pdb
    1A2K_r_u.pdb        
    1ACB_l_b.pdb
    1ACB_l_u.pdb
    1ACB_r_b.pdb
    1ACB_r_u.pdb 
    
Then, edit the following fields in ./configFile/configFile.cfg

  ncpu: int. number of cpu's to run in parallel (subprocess for features computing and threads for model training)
  modelType: "mixed" or "seq". type of model you want to train, sequence-only (seq) or sequence and structure (mixed)
  N_KFOLD: int. Type of cross validation. -1 for leave-one-complex out, possitive values for k= N_KFOLD cross-validation
  psiBlastNThrs: int. number of threads to use in psiblast
  minNumResiduesPartner: Minimum number of amino acids of a partner
  maxNumResiduesPartner: Maximum number of amino acids of a partner
  
  
  pdbsIndir: path where pdb files used to train benchmark are stored (can be removed after training)
  computedFeatsRootDir: directory where features files will be stored as subdirectories (can be removed after training)
  codifiedDataRootDir: str. Directory where ready to train joblib pickle files will be stored (can be removed after training)
  resultsRootDir: str. Directory where cross validation results will be stored
  savedModelsPath: str. Directory where xgBoost models will be saved.

  psiBlastDB: path where psiblast uniref90 database is placed
  
Next, load anaconda environment
  source activate xgbpred

Finally execute python script
  python generateBIPSPIModel.py
NOTE: tmux or screen are recommended when training the model.
  e.g. screen -dmSL trainSession python generateBIPSPIModel
  
  
*** 2.2 Predict ***  

  In order to obtain predictions you need a set of pdb or fasta files stored in a directory. 
  For each complex, 2 files must be provided, one for ligand and other for the receptor partner.
  have four different files. filenames follow this rule: prefix_X_u.Y, where prefix is an id for the complex (a pdb id or any other
  unique adress), X is l or r (ligand or receptor) and Y is .pdb or .fasta .
  
  For example:
  ~/path/to/predictSequences/
      1ACB_l_u.fasta
      1ACB_r_u.fasta
      seq1_l_u.fasta
      seq1_r_u.fasta
  or
  
  ~/path/to/predictPDBs/
      1ACB_l_u.pdb
      1ACB_r_u.pdb
      c1_l_u.pdb
      c1_r_u.pdb  
If files are pdbs, sequence-based and structural features are used, otherwise, sequence-based features.

Then, edit the following fields in ./configFile/configFile.cfg

  ncpu: int. number of cpu's to run in parallel (subprocess for features computing and threads for model training/prediction)
  savedModelsPath: str. Directory where xgBoost models are loaded. Already trained models are located at 
                        ~/xgbModels
  psiBlastNThrs: int. number of threads to use in psiblast
  psiBlastDB: path where psiblast uniref90 database is placed
  minNumResiduesPartner: Minimum number of amino acids of a partner
  maxNumResiduesPartner: Maximum number of amino acids of a partner
#The following filds are just used in  training and thus, ignored
  modelType: Ignored
  N_KFOLD: Ignored
  pdbsIndir: Ignored
  computedFeatsRootDir: Ignored
  codifiedDataRootDir: Ignored
  resultsRootDir: str. Ignored

Next, load anaconda environment
  source activate xgbpred

Finally execute python script
  python predictComplexes.py path/where/inputFiles/areLocated path/where/predictions/are/stored/path/to/results
NOTE: tmux or screen are recommended when predicting several complexes
  e.g. screen -dmSL trainSession python predictComplexes.py path/where/inputFiles/areLocated path/where/predictions/are/stored/path/to/results
   
   
For each complex, 3 results file are generated. in path/where/predictions/are/stored/path/to/results/preds
  -prefix.tab.res: predition of Residue-Residue Contacts. Has the following columns
        chainIdL structResIdL resNameL chainIdR structResIdR resNameR categ prediction
        categ colum is ignored.
        predictions go from 0 to 1, 1 contact, 0 no contact.
        
  -prefix.tab.res.lig: predition of ligand binding site. Has the following columns
        chainId resId categ prediction
        categ colum is ignored.
        predictions go from 0 to +infinite, 0 no binding site        

  -prefix.tab.res.rec: predition of receptor binding site. Has the following columns
        chainId resId categ prediction
        categ colum is ignored.
        predictions go from 0 to +infinite, 0 no binding site  
