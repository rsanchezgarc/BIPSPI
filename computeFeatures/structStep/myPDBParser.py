from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.MMCIFParser import MMCIFParser 


class myPDBParser(PDBParser, MMCIFParser):
  def __init__(self, structure_builder=None, QUIET=False):
    PDBParser.__init__(self, structure_builder= structure_builder, QUIET=QUIET)
    MMCIFParser.__init__(self, structure_builder= structure_builder, QUIET=QUIET)
  
  def get_structure(self, structure_id, filename):
    try:
      return PDBParser.get_structure( self, structure_id, filename)
    except Exception:
      return MMCIFParser.get_structure(self, structure_id, filename)
    
    
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
