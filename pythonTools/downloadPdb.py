import os, sys
import requests
from subprocess import Popen, PIPE


from Bio.PDB.mmtf import MMTFParser
from Bio.PDB.PDBIO import PDBIO

try:
  from myPDBParser import myPDBParser as PDBParser
except ImportError:
  sys.path.append(os.path.split(__file__)[0])
  from myPDBParser import myPDBParser as PDBParser

from extractPDBChain import ChainSplitter, NoChainInPdb

try:
  from bigExceptions import NoValidPDBFile
except ImportError:
  class NoValidPDBFile(Exception): pass
  
try:
  from urllib2 import HTTPError
except ImportError:
  from urllib.error import HTTPError

MAX_NUMBER_OF_CHAINS=48

def downloadUsingMmtf(pdbId, fnameOut, maxNumberOfChains= MAX_NUMBER_OF_CHAINS):
  print("downloadUsingMmtf")
  try:
    parser= MMTFParser()
    struct= parser.get_structure_from_url(pdbId)
    if not 0 in struct:
      return False
    if len(struct[0])>maxNumberOfChains:
      raise NoValidPDBFile("The maximun number of allowed chains is %d (%d) for %s"%(maxNumberOfChains, len(struct[0]), pdbId))
    writter=PDBIO()
    writter.set_structure(struct)
    writter.save(fnameOut)
    return True  
  except (Exception, ValueError, HTTPError) as e:
    print(e)
    if isinstance(e, NoValidPDBFile):
      raise e   
    return False

def checkIfSuccess(fname, maxNumberOfChains=MAX_NUMBER_OF_CHAINS):
  try:
    parser= myPDBParser()
    struct= parser.get_structure("pdb", fname)
    if not 0 in struct:
      return False
    if len(struct[0])>maxNumberOfChains:
      raise NoValidPDBFile("The maximun number of allowed chains is %d (%d) for %s"%(maxNumberOfChains, len(struct[0]), pdbId))
    return True
  except (Exception, ValueError) as e:
    print(e)
    if isinstance(e, NoValidPDBFile):
      raise e   
    return False
    
def downloadPDB(pdbId, pdbOutPath, chainId= None, bioUnit=None, removeOtherChains=False, checkObsolete=True):

  if not len(pdbId)==4 :
      raise ValueError("bad format pdbId '%s'"%(pdbId))
      
  pdbId= pdbId.lower()    
  outName= os.path.join(pdbOutPath,pdbId+'.pdb')
  print("Trying to download pdbId %s to %s"%(pdbId, outName))
  success=False
  if not os.path.isfile(outName):
    if not bioUnit: 
      success= downloadUsingMmtf(pdbId, outName)
    if not success or bioUnit:
      if checkObsolete:
        r= requests.get("https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/%s"%(pdbId))
        if r.status_code!=200:
          raise ValueError("pdb id %s does not exist according to https://www.ebi.ac.uk/pdbe/api/pdb/entry/summary/%s. It may be obsolete or cif"%(pdbId,pdbId))
      if bioUnit:
        cmd= 'wget -qO- ftp://ftp.wwpdb.org/pub/pdb/data/biounit/coordinates/all/%s.pdb%d.gz --timeout=30 | zcat  > %s'%(
                                                                      pdbId, bioUnit, outName)
      else:
        cmd= 'wget -qO- http://www.rcsb.org/pdb/files/%s.pdb.gz --timeout=30 | zcat  > %s'%(pdbId, outName)
      print(cmd)
      p= Popen(cmd, stdin=PIPE, stdout= PIPE, stderr=PIPE, shell=True)
      out= p.communicate()
      if len(out[1])>0 or not checkIfSuccess(outName):
        print(out)
        try:
          os.remove(outName)
        except OSError:
          pass
        msg= "Error downloading pdb_id: %s biounit: %s. Does biounit exists?"%(pdbId, bioUnit) if bioUnit else \
              "Error downloading pdb_id: %s"%pdbId
        raise NoValidPDBFile( msg )
  else:
    print("Already downloaded pdb_id %s"%(pdbId))
  if not chainId is None:
    print("extracting chain '%s'"%(chainId))
    try:
      splitter = ChainSplitter( os.path.split(outName)[0])
      outName_extract= splitter.make_pdb(outName, chainId, rejectInsteadAccept=False)
#      print(outName_extract, outName)
      if removeOtherChains:
        os.remove(outName)
    except NoChainInPdb:
      raise        
    except (ValueError, OSError):
      try:
        os.remove(outName_extract)
      except (OSError, UnboundLocalError):
        pass
      raise
    return outName_extract
  else:
    return outName    


def test():
  downloadPDB("1atn", "/tmp/downloadTest", chainId= 'A', bioUnit=0, removeOtherChains=False)
  
if __name__=="__main__":
#  test(); sys.exit()
  usage="python downloadPdb.py pathWherePdbWillBeSaved pdbId [bioUnitNumber] \npdbId: fourLeters[:ChainId]  e.g. 1A2K; 1A2K:B"
  if len(sys.argv)>=2:
    pdbOutPath= os.path.abspath(os.path.expanduser(sys.argv[1]))
    pdbId= sys.argv[2]
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
    downloadPDB(pdbId, pdbOutPath, bioUnit)
  except ValueError as e:
    print(str(e))
