
import os, sys
from myPDBParser import myPDBParser
from utils import tryToRemove, myMakeDir
from Bio.PDB import Structure, PDBIO
from subprocess import check_output
import random
import getopt
import glob
import socket

DIST_CUT= 0.8 #in nm
MAX_ELEMS_PER_CLUS= 10
if socket.gethostname()=="servet":
  GMX_PATH="gmx"
else:
  GMX_PATH="/mnt/big2/jsegura/app/gromacs-5.0.7/bin/gmx"
  
def mergePDBs(filesPattern, mergedFileName):
  fileNames=[]
  if not os.path.isfile(mergedFileName):
    check_output(["touch", mergedFileName])
    for i, fileName in enumerate(glob.iglob( os.path.expanduser(filesPattern))):
      print(fileName)
      fileNames.append(fileName)
      check_output('echo "MODEL      %d" >> %s'%(i, mergedFileName), shell=True)
      if fileName.endswith(".gz"):
        check_output('zcat %s >> %s'%(fileName, mergedFileName), shell=True)
      else:
        check_output('cat %s >> %s' % (fileName, mergedFileName), shell=True)
    check_output('echo "ENDMDL" >> %s'%(mergedFileName), shell=True)
    check_output('echo "END" >> %s'%(mergedFileName), shell=True)      
  else:
    return list(glob.iglob( os.path.expanduser(filesPattern)))
  return fileNames

def processOneCluster(fileNames, centroid, members, outPath, maxPerClus):
  print(members)
  raw_input("enter")
  fnamesList=[ fileNames[elem] for elem in members if elem!=centroid]
  fnamesList= [fileNames[centroid]]+ list(random.sample(fnamesList, min(len(fnamesList),max(1,maxPerClus-1))) )
  print(fnamesList)
  for fname in fnamesList:
    baseName= os.path.basename(fname)
    print(fname,  os.path.join(outPath, baseName) )
    symLinkName=  os.path.join(outPath, baseName)
    if not os.path.isfile(symLinkName):
      os.symlink(fname, symLinkName)
    if "_l_" in baseName :
      oriChainType="_l_"
      newChainType="_r_"
    else:
      oriChainType="_r_"
      newChainType="_l_"
    oriPath, oriBase= os.path.split(fname) 
    fname= os.path.join(oriPath, oriBase.replace(oriChainType, newChainType) )
    symLinkName=  os.path.join(outPath, baseName.replace(oriChainType, newChainType) )
    if not os.path.isfile(symLinkName):
      print(fname, symLinkName)
      os.symlink(fname, symLinkName)
      
def clusterize(mergedFileName, prefix, fileNames, outPath, maxPerClus):
  logsClusters= os.path.join( os.path.expanduser("~/tmp"), prefix+".clusters.log")
  distCut= DIST_CUT
  gmx= GMX_PATH
  cwd= os.getcwd()
  os.chdir(os.path.expanduser("~/tmp"))
  cmd= ("echo 5 | %(gmx)s cluster -f %(mergedFileName)s -s %(mergedFileName)s "+
               "-cutoff %(distCut)f -g %(logsClusters)s -nofit -method gromos")%locals()
               
  if not os.path.isfile(logsClusters):
    print(cmd)
    try:
      check_output(cmd, shell=True)
      os.chdir(cwd)
    except (Exception, KeyboardInterrupt):
      tryToRemove( logsClusters) 
  with open(logsClusters) as f:
    for line in f:
      if line.startswith("cl."):
        break
    members=None
    for line in f:
      lineArray= line.split("|")
      print(lineArray)
      if not  lineArray[0].strip().isdigit():
        members+= [ int(elem) for elem in lineArray[-1].split() ]
      else:
        if members:
          processOneCluster(fileNames, centroid, members, outPath, maxPerClus)
        centroid= lineArray[2].split(".")[0]
        if centroid[-1].isdigit():
          centroid= centroid[:-1]
        centroid= int(centroid)
        members= [ int(elem) for elem in lineArray[-1].split() ]
  processOneCluster(fileNames, centroid, members, outPath, maxPerClus)
          
def getClusterRepresentatives(prefix, filesPattern, outPath, maxPerClus= MAX_ELEMS_PER_CLUS):
  outPath= myMakeDir (os.path.expanduser( outPath))
  if maxPerClus is None:
    path, base = os.path.split(os.path.expanduser(filesPattern))
    nPos= check_output("ls %s/*_T* | wc -l "%(path), shell=True)
    maxPerClus= int(nPos)/2

  mergedFileName= os.path.join( os.path.expanduser("~/tmp"), prefix+".merged.pdb")
  try:
    fileNames= mergePDBs(filesPattern, mergedFileName)
  except (Exception, KeyboardInterrupt) as e:
    print("exception:", e)
    tryToRemove( mergedFileName)
    raise
  print(fileNames)
  clusterize(mergedFileName, prefix, fileNames, outPath, maxPerClus)

if __name__== "__main__":

  options, remainder= getopt.getopt(sys.argv[1:], 'i:o:p:',
                                    ['input=','output=','prefix=',])

  #    print("options:",cmds[1:], options, remainder)
  #    raw_input()
  inputPattern=None
  outPath= None
  prefix=None
  for opt, arg in options:
    if opt in   ('-i', '--input'):
      inputPattern = os.path.abspath(os.path.expanduser(arg))
    elif opt in ('-o', '--output'):
      outPath = os.path.abspath(os.path.expanduser(arg))
    elif opt in ('-p', '--prefix'):
      prefix= arg
  print(options)
  assert inputPattern and outPath and prefix, "Example:\npython -m pythonTools.clusterizePdbs -p 1A2K_l -o \"~/Tesis/dockScoring/data/develData/sampledDecoys/1A2K*pdb\" -i  ~/Tesis/dockScoring/data/develData/decoys/1A2K/1A2K_A_l_.pdb.gz"
  getClusterRepresentatives(prefix, inputPattern, outPath)

