################################################################################
# BIPSPI: xgBoost Interface Prediction of Specific Partner Interactions        #
################################################################################

ACADEMIC USE ONLY. This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY

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
  channels:
    - https://conda.binstar.org/travis
    - kalefranz
    - samoturk
    - conda-forge
    - bioconda
    - anaconda
    - defaults
  dependencies:
    - alabaster=0.7.10=py27_0
    - astroid=1.4.9=py27_0
    - babel=2.3.4=py27_0
    - backports=1.0=py27_0
    - backports_abc=0.5=py27_0
    - bleach=1.5.0=py27_0
    - cairo=1.14.8=0
    - chardet=2.3.0=py27_0
    - configparser=3.5.0=py27_0
    - cycler=0.10.0=py27_0
    - cython=0.26=py27_0
    - dbus=1.10.10=0
    - decorator=4.0.11=py27_0
    - docutils=0.13.1=py27_0
    - entrypoints=0.2.2=py27_1
    - enum34=1.1.6=py27_0
    - expat=2.1.0=0
    - fontconfig=2.12.1=3
    - freetype=2.5.5=2
    - funcsigs=1.0.2=py27_0
    - functools32=3.2.3.2=py27_0
    - get_terminal_size=1.0.0=py27_0
    - glib=2.50.2=1
    - gst-plugins-base=1.8.0=0
    - gstreamer=1.8.0=0
    - html5lib=0.999=py27_0
    - icu=54.1=0
    - imagesize=0.7.1=py27_0
    - ipykernel=4.5.2=py27_0
    - isort=4.2.5=py27_0
    - jbig=2.1=0
    - jedi=0.9.0=py27_1
    - jinja2=2.9.5=py27_0
    - joblib=0.11=py27_0
    - jpeg=8d=2
    - jsonschema=2.5.1=py27_0
    - lazy-object-proxy=1.2.2=py27_0
    - lcms=1.19=0
    - libffi=3.2.1=1
    - libgcc=5.2.0=0
    - libgfortran=3.0.0=1
    - libiconv=1.14=0
    - libpng=1.6.27=0
    - libsodium=1.0.10=0
    - libtiff=4.0.6=2
    - libxcb=1.12=1
    - libxml2=2.9.4=0
    - llvmlite=0.19.0=py27_0
    - markupsafe=0.23=py27_2
    - matplotlib=2.0.0=np112py27_0
    - mistune=0.7.3=py27_0
    - mkl=2017.0.1=0
    - msgpack-python=0.4.8=py27_0
    - nbconvert=5.1.1=py27_0
    - nbformat=4.3.0=py27_0
    - numba=0.34.0=np112py27_0
    - numpydoc=0.6.0=py27_0
    - openssl=1.0.2k=0
    - pandocfilters=1.4.1=py27_0
    - path.py=10.1=py27_0
    - pathlib2=2.2.0=py27_0
    - patsy=0.4.1=py27_0
    - pcre=8.39=1
    - pep8=1.7.0=py27_0
    - pexpect=4.2.1=py27_0
    - pickleshare=0.7.4=py27_0
    - pil=1.1.7=py27_2
    - pillow=3.2.0=py27_0
    - pip=9.0.1=py27_1
    - pixman=0.34.0=0
    - prompt_toolkit=1.0.9=py27_0
    - psutil=5.2.0=py27_0
    - ptyprocess=0.5.1=py27_0
    - pycairo=1.10.0=py27_0
    - pyflakes=1.5.0=py27_0
    - pygments=2.2.0=py27_0
    - pylint=1.6.4=py27_1
    - pyparsing=2.1.4=py27_0
    - pyqt=5.6.0=py27_2
    - python=2.7.13=0
    - python-dateutil=2.6.0=py27_0
    - pytz=2016.10=py27_0
    - pyzmq=16.0.2=py27_0
    - qt=5.6.2=2
    - qtawesome=0.4.4=py27_0
    - qtconsole=4.2.1=py27_1
    - qtpy=1.2.1=py27_0
    - readline=6.2=2
    - reportlab=3.4.0=py27_0
    - requests=2.13.0=py27_0
    - rope=0.9.4=py27_1
    - scandir=1.5=py27_0
    - setuptools=27.2.0=py27_0
    - simplegeneric=0.8.1=py27_1
    - simplejson=3.11.1=py27_0
    - singledispatch=3.4.0.3=py27_0
    - sip=4.18=py27_0
    - six=1.10.0=py27_0
    - snowballstemmer=1.2.1=py27_0
    - sphinx=1.5.1=py27_0
    - spyder=3.1.3=py27_0
    - sqlalchemy=1.1.13=py27_0
    - sqlite=3.13.0=0
    - ssl_match_hostname=3.4.0.2=py27_1
    - statsmodels=0.8.0=np112py27_0
    - subprocess32=3.2.7=py27_0
    - testpath=0.3=py27_0
    - tk=8.5.18=0
    - tornado=4.4.2=py27_0
    - traitlets=4.3.2=py27_0
    - wcwidth=0.1.7=py27_0
    - wget=1.18=0
    - wheel=0.29.0=py27_0
    - wrapt=1.10.8=py27_0
    - xlrd=1.1.0=py27_0
    - xz=5.2.2=0
    - zeromq=4.1.5=0
    - zlib=1.2.8=3
    - biopython=1.70=np112py27_0
    - mmtf-python=1.0.2=py27_0
    - amqp=2.1.4=py27_0
    - anyjson=0.3.3=py27_1
    - billiard=3.5.0.2=py27_0
    - blas=1.1=openblas
    - h5py=2.7.0=np112py27_1
    - hdf5=1.8.18=0
    - kombu=4.1.0=py27_0
    - numpy=1.12.1=py27_blas_openblas_200
    - openblas=0.2.19=2
    - pandas=0.21.0=py27_0
    - scikit-learn=0.19.1=py27_blas_openblas_200
    - scipy=0.19.1=py27_blas_openblas_202
    - vine=1.1.4=py27_0
    - xgboost=0.6a2=py27_2
    - uwsgi=2.0.2=py27_0
    - httpd=2.2.29=6
    - freeglut=3.0.0=4
    - glew=2.0.0=0
    - pmw=2.0.0=py27_1
    - pymol=1.8.6.0=py27_0
    - pip:
      - backports-abc==0.5
      - backports.shutil-get-terminal-size==1.0.0
      - backports.ssl-match-hostname==3.4.0.2
      - boto==2.46.1
      - bs4==0.0.1
      - bz2file==0.98
      - gensim==2.0.0
      - gputil==1.3.0
      - httplib2==0.9.2
      - mock==2.0.0
      - olefile==0.44
      - pbr==1.10.0
      - prompt-toolkit==1.0.7
      - protobuf==3.2.0
      - smart-open==1.5.2

Otherwise, you should install all them manually and edit
./configFiles/configFile.cfg file consequently to point to installation location
  
  
-------------------------
- 2. Use                -
-------------------------

*** 2.1 Train Model ***
  
In order to train a model you need a set of protein complexes with the format of Docking Benchmark v5 stored in a 
directory. For each complex, 4 pdb files must be provided, 2 for ligand (bound and unbound state) and
2 for receptor (bound and unbound). If just bound pdb files available, you must symlink them in order to
have four different files. filenames are prefix_X_Y.pdb, where prefix is an id for the complex (a pdb id or any other
unique adress), X is l or r (ligand or receptor) and Y is u or b (bound or unbound).

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
  have four different files. filenames are prefix_X_u.Y, where prefix is an id for the complex (a pdb id or any other
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
