# Repository usage guide

### Content:
 
 1. [INSTALLATION](#installation)
 
 2. [USAGE](#usage)  

    2.1. [Train model](#Train-model) 
    
    2.2. [Predict](#Predict)
  


### Installation


1) Install dependencies. The following bioformatics packages are required:
    - PSAIA 1.0 (http://bioinfo.zesoi.fer.hr/index.php/en/10-category-en-gb/tools-en/19-psaia-en)
    - DSSP  (https://swift.cmbi.umcn.nl/gv/dssp/index.html)
    - AL2CO  (http://prodata.swmed.edu/al2co/al2co.php)
          Warning. Default code has too small buffers for input names, code was modified from char[500]  to char[1024] and
          recompiled
    
    - cd-hit 4.6 (http://weizhongli-lab.org/cd-hit/) (AL2CO dependency)
    - clustalw 1.83 (http://www.clustal.org/) (AL2CO dependency)
    - psiblast 2.2.31+ (https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Web&PAGE_TYPE=BlastDocs&DOC_TYPE=Download)
    - SPIDER2 (http://sparks-lab.org/yueyang/server/SPIDER2/)

    You need to install all them and edit the [dependencies.cfg](../configFiles/cmdTool/dependencies.cfg) (`./configFiles/cmdTool/dependencies.cfg`) file 
    accordingly to point to installation locations. If other config file is intented to be used, it can be set using the `--configDependencies` flag.

2) Prepare BLAST database
    1. Create a directory that will contain the database
     ```
    mkdir /blast/database/directory/
    cd /blast/database/directory/
    
    ```   
    2. Donwload the raw sequences and decompress them.
    ```
    wget ftp://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref90/uniref90.fasta.gz
    gunzip -v uniref90.fasta.gz
    ```
    3. Compile database
    ```
     makeblastdb -in uniref90.fasta -dbtype prot -out uniref90.fasta -hash_index
    ```
    4. Edit  [dependencies.cfg](../configFiles/cmdTool/dependencies.cfg) (`./configFiles/cmdTool/dependencies.cfg`) file 
    to point to the BLAST database. For this example, the name is `/blast/database/directory/uniref90.fasta`


3) Install Python dependencies using conda.

```
conda env create -f  bipspiV2_environ.yml
```


At this point you should be able to execute BIPSPI v2.
In order to validate if your installation works, you could execute BIPSPI with its testing complexes:
```
python generateBIPSPIModel.py --N_KFOLD 4 --scopeFamiliesFname ./docs/scopes_example.tab --wdir /tmp/test_bipspi --pdbsIndir ./docs/trainingPDBsExample --ncpu 2
```

### Usage            


 ##### Train Model
  
In order to train a model you need a set of protein complexes with the format of the Protein Docking Benchmark v5 stored in a 
directory. For each complex, 4 pdb files must be provided, 2 for ligand (bound and unbound state) and
2 for receptor (bound and unbound). If only bound pdb files available, you must symlink them in order to
have four different files, pretending that the unbound file is the same that the bound one.
Filenames follow the pattern prefix_X_Y.pdb, where prefix is an id for the complex (a pdb id or any other
unique id), X is 'l' or 'r' (ligand or receptor) and Y is 'b' or 'u' (bound or unbound). Finally, the prefix must
be written in capital letters if the complex needs to be evaluated, while written in lowercase if the complex is
a train-only one.

You can find an example in [trainingPDBsExample](docs/trainingPDBsExample):



<details>
  <summary>For example:</summary>

  
```
$ls -l docs/trainingPDBsExample
total 5020
-rw-r--r-- 1 user user 129438 Oct 13 23:01 1A2K_l_b.pdb
-rw-r--r-- 1 user user 129643 Oct 13 23:01 1A2K_l_u.pdb
-rw-r--r-- 1 user user 161190 Oct 13 23:01 1A2K_r_b.pdb
-rw-r--r-- 1 user user 156345 Oct 13 23:01 1A2K_r_u.pdb
-rw-r--r-- 1 user user  42282 Oct 13 23:02 1ACB_l_b.pdb
-rw-r--r-- 1 user user  45429 Oct 13 23:02 1ACB_l_u.pdb
-rw-r--r-- 1 user user 143289 Oct 13 23:02 1ACB_r_b.pdb
-rw-r--r-- 1 user user 142125 Oct 13 23:02 1ACB_r_u.pdb

...

-rw-rw-r-- 1 user user 285289 Jul 21  2020 3qdg_l_b.pdb
lrwxrwxrwx 1 user user     12 Oct 14 01:18 3qdg_l_u.pdb -> 3qdg_l_b.pdb
-rw-rw-r-- 1 user user 252402 Jul 21  2020 3qdg_r_b.pdb
lrwxrwxrwx 1 user user     12 Oct 14 01:18 3qdg_r_u.pdb -> 3qdg_r_b.pdb

```

</details>
<br>

Default training mode is K-fold cross-validation. The number of folds is set  in [dependencies.cfg](../configFiles/cmdTool/dependencies.cfg) (`./configFiles/cmdTool/dependencies.cfg`) file parameter
N_KFOLD. If -1, leave one out is performed. In the last step, only predictions for capital letter prefixes will be produced.
If all complexes are wanted to be predicted, edit `SKIP_LOWER_PREDICTION=False` in [trainAndTest.py](../trainAndTest/trainAndTest.py)

Alternatively, you can provide a list of complexes to be considered in each fold. In that case N_KFOLD will point
towards that file.


Additionally, for cross-validation, you generally require a file with information about the families of the protein chains to ensure that the cross-validation partitions are independend.
In our study, we employ SCOPe familes, but any grouping strategy should be valid (e.g sequence identity).
An example of such file is found in [scopes_example.tab](../docs/scopes_example.tab) and the scheme is as follows:
```
COMPLEX_ID  CHAINS_PARTNER_1  CHAINS_PARTNER_2 SCOPES_PARTNER_1 SCOPES_PARTNER_2

###################################      DESCRIPTION        ##########################
COMPLEX_ID: lowercase if train-only complex and upper-case if test-complex
CHAINS_PARTNER_1 and CHAINS_PARTNER_2: The chains to consider in each complex. No id chains is denoted as '+'. If several chains, they are concatenated with no delimiter
SCOPES_PARTNER_1 and SCOPES_PARTNER_1:  The family/cluster ids. Familes of different chains are separeted by ';'. Within the same chain, familes are separated by ":". 

```

 
The config file for the training experiment  [configFile.cfg](../configFiles/cmdTool/configFile.cfg) (`./configFiles/cmdTool/dependencies.cfg`)  
needs to be eddited to point towards the location of the train complexes and the scopes families file.
Another config file located in a non estandard location could be used with  the flag `--configFile path/to/configFile.cfg)`

The following parameters are the most important config parameters for training:

<details>
  <summary>PARAMETERS:</summary>

  
```
#ncpu: int. number of cpu's to run in parallel (subprocess for features computing and threads for model training/prediction)
ncpu 1
#modelType: "struct", "mixed" or "seq". type of model you want to train, sequence-only (seq) or sequence and structure (struct) or one partner seq and the other struct (mixed)
modelType struct
#checkHomoInteractionInTraining Corrects contact maps to consider homologous residues as contact if one does. Set it to True if using a homo-dataset for training. For heter-datasets impact is minor
checkHomoInteractionInTraining True
#N_KFOLD: int. Type of cross validation. -1 for leave-one-complex out, positive values for k= N_KFOLD cross-validation. If a path to file, used as a train test list file
N_KFOLD 10
#N_KFOLD_path ./docs/customFolds.txt
#scopeFamiliesFname_path: str. (optional) path to scope families to guaranty independency
scopeFamiliesFname_path ./docs/scopes_example.tab

#pdbsIndir: str. path where pdb files used to train benchmark are stored (can be removed after training)
pdbsIndir_path ./docs/trainingPDBsExample

#wdir_path: the working directory for the project
wdir_path ~/tmp/prueba_bipspi
#computedFeatsRootDir_path: str. directory where features files will be stored as subdirectories (can be removed after training)
computedFeatsRootDir_path %(wdir)s/computedFeatures
#codifiedDataRootDir_path: str. Directory where ready to train joblib pickle files will be stored (can be removed after training)
codifiedDataRootDir_path %(wdir)s/codifiedInput
#resultsRootDir_path: str. Directory where cross validation results will be stored
resultsRootDir_path %(wdir)s/results
#savedModelsPath_path: str. Directory where xgBoost models will be saved.
#Already trained models are stored in:
#savedModelsPath_path ~/Tesis/rriPredMethod/pyCode/webApp/rriPredWeb/media/xgbModels
savedModelsPath_path %(wdir)s/modelsComputed

#Minimun size and maximun size of partners for each complex
minNumResiduesPartner 12
maxNumResiduesPartner 2999

#temporal directory
tmp_path ~/tmp

```

</details>
<br>

Next, load anaconda environment
  ```
  conda activate bipspiV2
  ```

Finally execute Python script:
```
python generateBIPSPIModel.py  #Will run according config.cfg file
```
NOTE: tmux or screen are recommended when training the model to prevent network issues
  ``` 
  screen -dmSL trainSession python generateBIPSPIModel.py
  
  or
  
  screen -dmSL trainSession python generateBIPSPIModel.py --configFile /path/to/myConfigFile
  ```
  
  
BISPI v2 allows to train using PDB-PDB, PDB-Sequence or Sequence sequence. Edit the config file `modelType` field to select
the model you want to train.

Lastly, if you don't want to edit the config  file, you can modify any of the options just by adding
to the command the `--OPTION VALUE` flag. E.g.:
```
python generateBIPSPIModel.py --N_FOLDS 5 #Perform a 5-fold cross-validation independently of what the config file says
python generateBIPSPIModel.py --wdir /path/to/results  --modelType mixed # Save results to /path/to/results and train a mixed model independently of what config file says
```


 ##### Predict  

  In order to obtain predictions you need a set of pdb or fasta files stored in a directory. 
  For each complex, 2 files must be provided, one for the ligand and other for the receptor partner.
  filenames are prefix_X_u.Y, where prefix is an id for the complex (a pdb id or any other
  unique identifier), X is l or r (ligand or receptor) and Y is .pdb or .fasta .
  
  For example:
  ```
  $ls ~/path/to/predictSequences/
      1ACB_l_u.fasta
      1ACB_r_u.fasta
      seq1_l_u.fasta
      seq1_r_u.fasta
  ```
  or
  ```
  $ls ~/path/to/predictPDBs/
      1ACB_l_u.pdb
      1ACB_r_u.pdb
      c1_l_u.pdb
      c1_r_u.pdb
```
or
```
  $ls ~/path/to/predictMixed/
      1ACB_l_u.pdb
      1ACB_r_u.fasta
      c1_l_u.pdb
      c1_r_u.fasta
```
If files are pdbs or pdbs and fasta, sequence-based and structural features are used, 
otherwise, sequence-based features only.

Then, edit the following fields in [configFile_pred.cfg](../configFile/configFile_pred.cfg) (`./configFile/configFile_pred.cfg`)


The following parameters are the most important config parameters for prediction:

<details>
  <summary>PARAMETERS:</summary>

```
PREDICT_ONLY True
#ncpu: int. number of cpu's to run in parallel (subprocess for features computing and threads for model training/prediction)
ncpu 1
#checkHomoInteractionInTraining Set it to True if using a homo-dataset for training
checkHomoInteractionInTraining True

#Already trained models are stored in:  #TODO: change this location
savedModelsPath ../webApp/rriPredWeb/media/xgbModels

#Minimun size and maximun size of partners for each complex
minNumResiduesPartner 12
maxNumResiduesPartner 2999

#temporal directory
tmp ~/tmp
```
</details>
<br>

If you don't want to edit the config file, you can always provide the options using the command line
flags `--OPTION value`.


Next, load anaconda environment
 ```
  conda activate bipspiV2
  ```

Finally execute python script

```
python predictComplexes.py -i ./docs/trainingPDBsExample -o path/where/predictions/are/stored/
```

You can also provide additional arguments using the `-OPTION value` instead of editing the config file
```
python predictComplexes.py -i ./docs/trainingPDBsExample --wdir  path/where/predictions/are/stored/ --ncpu 32 # Use 32 cpus.

```

Additionally, instead a direectory with files, you can provide a csv file with pdbIds and chain Ids to be used as input. 
BIPSPI will automatically download then and execute the struct mode. Only one chain per pdbId is superted at this moment.
See [predict_from_pdbIds.csv](../docs/predict_from_pdbIds.csv) for an example

```
python predictComplexes.py -f ./docs/predict_from_pdbIds.csv --wdir path/where/predictions/are/stored/

```



For each complex, 3 results file are generated in ` path/where/predictions/are/stored/results/(seq|struct|mixed)_2`
  ```
  -prefix.tab.res.gz: predition of Residue-Residue Contacts. Has the following columns
        chainIdL resIdL resNameL chainIdR resIdR resNameR categ prediction
        categ colum is ignored.
        predictions go from 0 to 1, 1 contact, 0 no contact.
        
  -prefix.tab.res.lig.gz: predition of ligand binding site. Has the following columns
        chainId resId categ prediction
        categ colum is ignored.
        predictions go from 0 to +infinite, 0 no binding site        

  -prefix.tab.res.rec.gz: predition of receptor binding site. Has the following columns
        chainId resId categ prediction
        categ colum is ignored.
        predictions go from 0 to +infinite, 0 no binding site  
```
