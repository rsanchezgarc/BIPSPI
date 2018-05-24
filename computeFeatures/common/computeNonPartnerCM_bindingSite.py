import os, sys
import pandas as pd
from Config import Configuration
from subprocess import Popen, PIPE, check_call
from ast import literal_eval as make_tuple

from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist

cd_hit_path="/home/rsanchez/Tesis/rriPredMethod/dependencies/bioinformaticTools/cdhit-master/cd-hit"

computedFeatsRootDir= Configuration().computedFeatsRootDir
sequencesPath_root= os.path.join(computedFeatsRootDir, "seqStep","extractedSeqs")
sequencesPath= os.path.join(sequencesPath_root, "seqsData")
resIdsMap= os.path.join(sequencesPath_root, "seqToStructMap")
cMapPath= os.path.join(computedFeatsRootDir, "common","contactMaps")
newCmPath= os.path.join(computedFeatsRootDir, "common","contactMapsBinding")

aligsFastaName= "/tmp/aligsFastaName.fa"
cdhit_out= "/tmp/cdhit.out"

scoreMat= matlist.blosum62


def computeCD_hit():
  cmd= "awk 1 %s/*.fasta  > %s"%(sequencesPath, aligsFastaName)
  print(cmd)
  check_call(cmd, shell=True)
  proc= Popen([cd_hit_path, "-i", aligsFastaName, "-o", cdhit_out])
  outCdhit= proc.communicate()
  #print(outCdhit)
  os.remove(aligsFastaName)
def load_seqs():
  seqsDict={}
  for fname in os.listdir(sequencesPath):
    if fname.endswith(".fasta"):
      pdbId, chainType, chainId, __= fname.split("_")
      smallName= "%s_%s_%s"%(pdbId, chainType, chainId)
#      print(smallName)
      with open(os.path.join(sequencesPath, fname)) as seqF:
        seq= seqF.readlines()[1]
#        print(seq)
      with open(os.path.join(resIdsMap, smallName+"_u.seqStruMap")) as seqMapF:
        seqMapF.readline()
        idsMap= [ make_tuple(elem.split(";")[-1]) for elem in seqMapF.readlines()]
        idsMap= [str(elem[1])+(elem[2].strip()) for elem in idsMap]
      seqsDict[smallName]= (seq, chainId, idsMap)
      
  return seqsDict      
    
def parseCD_hit():
  clustersList= []
  with open(cdhit_out+".clstr") as f:
    newCluster=[]
    insideCluster=False
    for line in f:
      if line.startswith(">"):
        newCluster=[]
        if insideCluster:
          clustersList.append(newCluster)
        insideCluster=True
        continue
      if insideCluster:
        newCluster.append( line.split("...")[0].split(">")[-1])
  return clustersList

def loadCM(prefix):
#  print(cMapPath, prefix+".cMap.tab")
  fname=os.path.join(cMapPath, prefix+".cMap.tab")
  df= pd.read_table(fname,sep='\s+', header='infer', comment="#", 
                        dtype= {"chainIdL":str, "chainIdR":str, "structResIdL":str, "structResIdR":str})
  df_ligand=   df[["chainIdL", "structResIdL", "categ"]].drop_duplicates()
  df_receptor= df[["chainIdR", "structResIdR", "categ"]].drop_duplicates()

#  print(df_ligand.head(), list(df_ligand["categ"]).count(1))
#  print(df_receptor.head())
#  if prefix.startswith("1A2K"):
#    raw_input()
  return df_ligand, df_receptor
  
def getDfSeqId(mapId, cmap, chain):
  cmap=cmap.loc[cmap.iloc[:,0]==chain,:].drop_duplicates()
  seqPositions=[]
#  print(cmap.tail(), mapId[-5:])
  for i in range(cmap.shape[0]):
    seqPositions.append( mapId.index(cmap.iloc[i,1]))
  cmapFin= cmap.copy()
  cmapFin["seqId"]= seqPositions
  return cmapFin


def processOneClus(listOfClusMembers):
  names, seqs, chains, mapIds= zip(* listOfClusMembers)
#  print(mapIds)
  contactMapsDict= {}
  for name, seq, chain, mapId in listOfClusMembers:
    prefix=name.split("_")[0]
    if prefix not in contactMapsDict:
      contactMapsDict[prefix]= loadCM(prefix)
#      print(prefix, contactMapsDict[prefix][0].head(), contactMapsDict[prefix][1].head())

  for name0, seq0, chain0, mapId0 in listOfClusMembers:
    prefix0=name0.split("_")[0]
    cmap0= contactMapsDict[prefix0][0] if "_l_" in name0 else contactMapsDict[prefix0][1]
    posCmap0= getDfSeqId(mapId0, cmap0, chain0)
    posCmap0= posCmap0[(posCmap0["categ"]==1) & (posCmap0[posCmap0.columns[0]]== chain0) ]
    for name1, seq1, chain1, mapId1 in listOfClusMembers:
      if name0== name1: continue
#      print(name1)
      prefix1=name1.split("_")[0]
      cmap1= contactMapsDict[prefix1][0] if "_l_" in name1 else contactMapsDict[prefix1][1]
      posCmap1= getDfSeqId(mapId1, cmap1, chain1)
      posCmap1= posCmap1[(posCmap1["categ"]==1) & (posCmap1[posCmap1.columns[0]]== chain1) ]
      alignments = pairwise2.align.localds(seq0, seq1, scoreMat, -11, -0.5)
      ali0= alignments[0][0]
      ali1= alignments[0][1]
      index0= -1
      index1=-1    
      for letter1, letter0 in zip(ali1, ali0):
        if letter1 !="-":
          index1+=1
        if letter0 !="-":
          index0+=1
        if letter1 !="-" and letter0 !="-":
          if sum(posCmap1["seqId"]==index0)>0:
            print("update cmap0", index0, mapId0[index0], name0, chain0)
            cmap0.loc[(cmap0[cmap0.columns[0]]==chain0) & (cmap0[cmap0.columns[1]]==mapId0[index0]), "categ"]= 1
#            print(cmap0.loc[(cmap0[cmap0.columns[0]]==chain0) & (cmap0[cmap0.columns[1]]==mapId0[index0]), :])
#            print(cmap0.head())
#            raw_input("press enter")
    cmap0_Chain= cmap0[cmap0[cmap0.columns[0]]==chain0]
    cmap0_Chain.columns= ["chainId", "structResId", "categ"]
    cmap0_Chain.to_csv(os.path.join(newCmPath,name0+"_u.binding.Cmap"), sep=" ", index=False)
#  raw_input("enter for next cluster")

def copyCmap(extendedPrefix):
  prefix, chainType, chainId =extendedPrefix.split("_")
  if chainType=="l":
    raw_cmap=loadCM(prefix)[0]
  else:
    raw_cmap=loadCM(prefix)[1]
  cmap= raw_cmap[raw_cmap[raw_cmap.columns[0]]==chainId]
  cmap.columns= ["chainId", "structResId", "categ"]
  cmap.to_csv(os.path.join(newCmPath,extendedPrefix+"_u.binding.Cmap"), sep=" ", index=False)

  return 
def processClusters(clusList, seqsDict):
  nClus=0
  for clus in clusList:
    print(clus)
    
#    if (not '1A2K_l_A' in clus): continue
    if len(clus)>1:
      nClus+=1
      processOneClus( [(smallName,) + seqsDict[smallName] for smallName in clus])
    else:
      copyCmap( clus[0] )
  print(nClus)
  
def keepDiferentPdbsClust(clustersList):
  diffClusList=[]
  for clus in clustersList:
    if len(clus)<2: continue
    firstElem=clus[0][:4]
    isValid=False
    for elem in clus[1:]:
      if elem[:4]!= firstElem:
        isValid=True
    if isValid:
      diffClusList.append(clus)
  return diffClusList
if __name__== "__main__":
  computeCD_hit()
  clustersList= parseCD_hit()
  print(keepDiferentPdbsClust(clustersList))
  raise ValueError("Stopped after cd-hit")
  seqsDict= load_seqs()
  processClusters(clustersList, seqsDict)
  
