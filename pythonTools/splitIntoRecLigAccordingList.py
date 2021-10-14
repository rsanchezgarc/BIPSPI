import sys, os
from extractPDBChain import ChainSplitter


#pdbListFile="/home/rsanchez/Tesis/rriPredMethod/data/joanDimers/117_dimers_list.tsv"
#pdbsPath="/home/rsanchez/Tesis/rriPredMethod/data/joanDimers/pdbFiles/rawPDBs"
#outPath="/home/rsanchez/Tesis/rriPredMethod/data/joanDimers/pdbFiles/pdbsForTrain"
chainL_idx=1
chainR_idx=2

def splitAll(pdbListFile, inPath, outPath):
  chainSelector = ChainSplitter(outPath)
  with open(pdbListFile) as f:
    for line in f:
      lineArray= line.split()
      pdbId= lineArray[0]
      print(pdbId)
      chainL= lineArray[chainL_idx]
      chainR= lineArray[chainR_idx]
      ligOutName= chainSelector.make_pdb(os.path.join(inPath, pdbId+".pdb"), chainL, rejectInsteadAccept=False)
      os.rename(ligOutName, os.path.join(outPath, pdbId+"-%s%s_l_u.pdb"%(chainL,chainR)))
      recOutName= chainSelector.make_pdb(os.path.join(inPath, pdbId+".pdb"), chainR, rejectInsteadAccept=False)
      os.rename(recOutName, os.path.join(outPath, pdbId+"-%s%s_r_u.pdb")%(chainL,chainR))
  
if __name__=="__main__":
  if len(sys.argv)!=4:
    raise ValueError("error, incorrect number of arguments")
  pdbListFile= os.path.expanduser( sys.argv[1])
  pdbsPath= os.path.expanduser( sys.argv[2])
  outPath= os.path.expanduser( sys.argv[3])
  splitAll(pdbListFile, pdbsPath, outPath)
  
