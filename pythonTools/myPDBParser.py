from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser
import re
import gzip
import requests
try:
    from StringIO import StringIO ## for Python 2
except ImportError:
    from io import StringIO ## for Python 3

class myPDBParser(PDBParser, MMCIFParser):
  def __init__(self, structure_builder=None, QUIET=False, removeHeteroDuplicated=True):
    self.removeHeteroDuplicated= removeHeteroDuplicated
    PDBParser.__init__(self, structure_builder= structure_builder, QUIET=QUIET)
    MMCIFParser.__init__(self, structure_builder= structure_builder, QUIET=QUIET)

  def get_structure(self, *args):
    if len(args)==2:
      pdbId, fileName= args
    elif len(args)==1:
      fileName= args[0]
      pdbId, fileName= str(fileName), fileName
    else:
      raise ValueError("Error, input should be (id, fileName) or (fileName))")

    if re.match("http(s?)://",fileName):
      r= requests.get(fileName)
      if r.ok:
        fileName= StringIO(r.text)
      else:
        raise Exception("Error downloading pdb")

    try:
      if not isinstance(fileName, str) or not fileName.endswith(".gz"):
        structure= PDBParser.get_structure( self,pdbId, fileName)
      else:
        with gzip.open(fileName) as f:
          structure= PDBParser.get_structure( self, pdbId, f)
    except Exception as e:
      print(e)
      structure= MMCIFParser.get_structure(self, pdbId, fileName)
    if self.removeHeteroDuplicated:
      structure= self.filterOutDuplicated(structure)
    return structure

  def filterOutDuplicated(self, structure):
    seenFlag= {}
    for model in structure:
      seenFlag[model]={}
      for chain in model:
        seenFlag[model][chain]=set([])
        for res in chain:
          resId_tuple= res.get_id()
          resNum= ("%d%s"%(resId_tuple[1], resId_tuple[-1])).strip()
          if resNum not in seenFlag[model][chain]:
            seenFlag[model][chain].add(resNum)
          else:
            #remove duplicated
            chain.detach_child(resId_tuple)
    return structure

def test():
  pdbFile="/home/rsanchez/Tesis/rriPredMethod/pyCode/webApp/rriPredWeb/media/inputExamples/1a79_l_u.pdb"
  mmCifFile="/home/rsanchez/Tesis/rriPredMethod/pyCode/webApp/rriPredWeb/media/inputExamples/1acb.cif"
  parser= myPDBParser(QUIET=True)
  struct= parser.get_structure("pdb", pdbFile)
  print(struct[0].get_list())
  struct= parser.get_structure("mmCif", mmCifFile)
  print(struct[0].get_list())
  return None

if __name__=="__main__":
  test()
