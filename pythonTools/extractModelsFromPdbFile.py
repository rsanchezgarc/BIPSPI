import sys, os
import re
import StringIO
from Bio.PDB import PDBParser, PDBIO
from Bio.PDB.Model import Model
from Bio.PDB.Structure import Structure


class StructWriter():
  def __init__(self, structId="subset"):
    self.structId= structId
    self.pdbParser = PDBParser(QUIET=True)
    self.structure = Structure(structId)

  def addModel(self, pdb_as_str, modelId):
    pdbLikeFile = StringIO.StringIO()
    pdbLikeFile.write(pdb_as_str)
    pdbLikeFile.flush()
    pdbLikeFile.seek(0, 0)
    # print( "--->", pdbLikeFile.getvalue())
    new_struct = self.pdbParser.get_structure(str(modelId), pdbLikeFile)
    print( new_struct.child_list )
    model= new_struct[0]
    model.detach_parent()
    model.id= int(modelId)
    model.serial_num= model.id
    model.get_full_id()
    self.structure.add(model)
    print( "Current struct", self.structure.child_list )

  def saveStruct(self, fname, desiredOrder):
    io= PDBIO(use_model_flag=True)
    if  desiredOrder is not None:
      children= self.structure.child_list

    self.structure = Structure(self.structId)
    for modelId in desiredOrder:
      child= [model for model in children if model.id==modelId][0]
      child.detach_parent()
      self.structure.add(child)

    io.set_structure(self.structure)
    io.save(fname) #,  preserve_atom_numbering=True)

  def __len__(self):
    return len(self.structure.child_list )

def extractModelsFromPdbFile(fnameDecoys, modelIds, outName):


  newStructWriter=  StructWriter()

  currentModelId=None
  pdb_str=""
  with open(fnameDecoys) as f:
    for i, line in enumerate(f):
      # print(line)
      matchObj=re.match("^MODEL\s+(\d+)", line)
      if matchObj:
        if currentModelId in modelIds:
          newStructWriter.addModel(pdb_str, currentModelId)
        currentModelId=matchObj.group(1)
        pdb_str=""
        if len(newStructWriter)>=len(modelIds): break
      else:
        pdb_str+=line

  if currentModelId in modelIds:
    newStructWriter.addModel(pdb_str, currentModelId)

  newStructWriter.saveStruct(outName, desiredOrder=map(int, modelIds))

if __name__=="__main__":
  fnameDecoys = os.path.expanduser(sys.argv[1])
  modelIds = sys.argv[2].split(",")
  outName= os.path.expanduser(sys.argv[3])

  assert len(modelIds) > 0, "Error, modelIds must be provided as argument"
  extractModelsFromPdbFile(fnameDecoys, modelIds, outName)