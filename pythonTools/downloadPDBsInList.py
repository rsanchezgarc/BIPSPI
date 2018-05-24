import sys, os
from subprocess import call


pdbListFile="/home/rsanchez/Tesis/rriPredMethod/data/joanDimers/117_dimers_list.tsv"
outPath="/home/rsanchez/Tesis/rriPredMethod/data/joanDimers/pdbFiles/rawPDBs"
USE_BIO_UNIT=False

def downloadPDB(pdbId, pdbOutPath, useBioUnit):
## descargar pdb: wget ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/1i1q.pdb2.gz  o ya descomprimido
## wget -qO- ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/1i1q.pdb2.gz |zcat  > 1i1q.pdb

  outName= os.path.join(pdbOutPath,pdbId+'.pdb')
  if not os.path.isfile(outName):
    if useBioUnit:
      cmd= 'wget -qO- ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/%s.pdb1.gz |zcat  > %s'%(pdbId.lower(), outName)
    else:
      cmd= 'wget -qO- http://www.pdb.org/pdb/files/%s.pdb  | cat > %s'%(pdbId.upper(), outName)
    print(cmd)
    call(cmd, shell= True)
    
def downloadInFile(fname, outPath, useBioUnit=True):
  with open(fname) as f:
    for line in f:
      pdbId= line.split()[0]
      print(pdbId)
      downloadPDB(pdbId, outPath, useBioUnit)
      
if __name__=="__main__":

  if len(sys.argv)==3:
    pdbListFile= os.path.abspath(os.path.expanduser(sys.argv[1]))
    outPath= os.path.abspath(os.path.expanduser(sys.argv[2]))
  print(  pdbListFile, outPath)
  downloadInFile(pdbListFile, outPath, USE_BIO_UNIT)
  
