import sys, os
from subprocess import call
from extractPDBChain import ChainSplitter

    
def spiltInFile(fname, inPath, outPath):
  splitter=ChainSplitter(outPath)
  with open(fname) as f:
    for line in f:
      pdbId, chainL, chainR= line.split()
      pdbOutName= splitter.make_pdb(os.path.join(inPath, pdbId+".pdb"), chainL, rejectInsteadAccept=False)
      os.rename(pdbOutName, os.path.join(outPath, pdbId+"-"+chainL+chainR+"_l_u.pdb"))
      pdbOutName= splitter.make_pdb(os.path.join(inPath, pdbId+".pdb"), chainR, rejectInsteadAccept=False)

      os.rename(pdbOutName, os.path.join(outPath, pdbId+"-"+chainL+chainR+"_r_u.pdb"))
      
if __name__=="__main__":

  if len(sys.argv)==4:
    pdbListFile= os.path.abspath(os.path.expanduser(sys.argv[1]))
    inPath= os.path.abspath(os.path.expanduser(sys.argv[2]))
    outPath= os.path.abspath(os.path.expanduser(sys.argv[3]))
  print(  pdbListFile, inPath, outPath)
  spiltInFile(pdbListFile, inPath, outPath)
  
