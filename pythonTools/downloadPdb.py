import os, sys
import requests
from subprocess import Popen, PIPE
from extractPDBChain import ChainSplitter

def downloadPDB(pdbid, pdbOutPath, chainId= None, bioUnit=None):

  if not len(pdbid)==4 :
      raise ValueError("bad format pdbid '%s'"%(pdbid))
      
  pdbid= pdbid.lower()    
  print("Trying to download pdbId %s"%(pdbid))
  r= requests.get("https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/%s"%(pdbid))
  if r.status_code!=200:
    raise ValueError("pdb id %s does not exist according to https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/%s"%(pdbid,pdbid))
  outName= os.path.join(pdbOutPath,pdbid+'.pdb')
  if not os.path.isfile(outName):
    if bioUnit:
      cmd= 'wget -qO- ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/%s.pdb%d.gz --timeout=30 | zcat  > %s'%(
                                                                    pdbid, bioUnit, outName)
    else:
      cmd= 'wget -qO- http://www.rcsb.org/pdb/files/%s.pdb.gz --timeout=30 | zcat  > %s'%(pdbid, outName)
    print(cmd)
    p= Popen(cmd, stdin=PIPE, stdout= PIPE, stderr=PIPE, shell=True)
    out= p.communicate()
#    print(out)
    if len(out[1])>0:
      try:
        os.remove(outName)
      except OSError:
        pass
      raise ValueError("Error downloading pdb_id: %s biounit: %s. Does biounit exist?"%(pdbid, bioUnit))    
    if not chainId is None:
      print("extracting chain '%s'"%(chainId))
      try:
        splitter = ChainSplitter( os.path.split(outName)[0])
        outName_extract= splitter.make_pdb(outName, chainId, rejectInsteadAccept=False)
        os.rename(outName_extract, outName)
      except ValueError:
        try:
          os.remove(outName)
        except OSError:
          pass
        raise
  else:
    print("Already computed pdb_id %s"%(pdbid))
  return outName

if __name__=="__main__":
  usage="python downloadPdb.py pathWherePdbWillBeSaved PDBID [bioUnitNumber] \nPDBID: fourLeters[:ChainId]  e.g. 1A2K; 1A2K:B"
  if len(sys.argv)>=2:
    pdbOutPath= os.path.abspath(os.path.expanduser(sys.argv[1]))
    pdbid= sys.argv[2]
  else:
    print(usage)
    raise ValueError("Bad arguments: %s"%(str(sys.argv)))
  bioUnit= None

  if len(sys.argv)==4:
    bioUnit= int(sys.argv[3])
  elif len(sys.argv)>4:
    print(usage)
    raise ValueError("Bad arguments: %s"%(str(sys.argv)))
  try:    
    downloadPDB(pdbid, pdbOutPath, bioUnit)
  except ValueError as e:
    print(str(e))
