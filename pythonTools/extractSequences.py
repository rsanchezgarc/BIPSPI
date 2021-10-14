from __future__ import absolute_import
import os,sys
from Bio.PDB.PDBParser import PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_one

SMALL_CHAINS_LIMIT=9

def extractSeqsOnePDB(struct, outNameFasta):
  chainToSeq={}
  for chain in struct[0]:
    chainId= chain.get_id()
    if chainId==" ":
      chainId="+"
    nResStandard= sum([ 1 for res in chain if is_aa(res,standard=True)])
    resList= [ res for res in sorted(chain.child_list, key=lambda x: x.get_id()[1:]) if is_aa(res, standard=False)]
    nResAll= len(resList)
    if nResStandard<= int(0.5 * nResAll): continue #skip if most residues are not standard
    if len(resList)> SMALL_CHAINS_LIMIT: #Too small chains will not be considered
      sequence=[]
      for i,res in enumerate(resList):
        try:
          letter= three_to_one(res.resname)
        except KeyError:
          print("Exception", res)
          letter= "X"
          if i== (nResAll-1): break #This case is for TCGR....TLRX where X is GDP or other molecule
        sequence.append(letter)
      sequence= "".join(sequence)
      chainToSeq[chain.get_id()]=(sequence, [ res for res in resList])
      if outNameFasta is not None:
        outNameFasta= outNameFasta%chainId
        with open(outNameFasta,"w") as f:
          f.write(">"+ os.path.basename(outNameFasta)+"\n"+sequence)
  return chainToSeq

def extractSeqsOnePDBFromFile(fnamePdb, prefix, outPath):
  if outPath is not None:
    outNameFasta= os.path.join(outPath,prefix+"_%(chainId)s.fasta")
  else:
    outNameFasta=None
  struct= PDBParser().get_structure( prefix, fnamePdb)
  extractSeqsOnePDB(struct, outNameFasta)

def extractAll(pdbsPath, seqsOutPath):

  badPrefixes= set([])
  for fname in os.listdir(pdbsPath):
    print(fname)
    prefix= fname.split(".")[0].split("_")
    if prefix[1]=="l" or prefix[1]=="r":
      prefix= "_".join(prefix[:2])
    else:
      prefix= prefix[0]
    try:
      extractSeqsOnePDB(os.path.join(pdbsPath, fname), prefix, seqsOutPath)
    except KeyError:
      badPrefixes.add(prefix)

  for fname in  pdbsPath:
    prefix= fname.split(".")[0].split("_")
    if prefix[1]=="l" or prefix[1]=="r":
      prefix= "_".join(prefix[:2])
    else:
      prefix= prefix[0]
    if prefix in badPrefixes:
      os.remove(os.path.join(pdbsPath,fname))

if __name__=="__main__":

#  pdbsPath="/home/rsanchez/Tesis/rriPredMethod/data/partnerSpecificity/trimers/rawPDBs"
#  seqsOutPath="/home/rsanchez/Tesis/rriPredMethod/data/partnerSpecificity/trimers/seqs"
  if len(sys.argv)==3:
    pdbsPath= sys.argv[1]
    seqsOutPath= sys.argv[2]
    extractAll(pdbsPath, seqsOutPath)
  else:
    print("usage:\npython extractSequences.py pdbsDirPath seqsOutPath")



