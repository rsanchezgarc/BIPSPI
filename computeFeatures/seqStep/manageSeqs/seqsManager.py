from __future__ import absolute_import
import os,sys
from computeFeatures.structStep.myPDBParser import myPDBParser as PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_one, d1_to_index
from ..SeqFeatComputer import SeqFeatComputer, FeatureComputerException
from utils import myMakeDir, tryToRemove #utils is at the root of the package

SMALL_CHAINS_LIMIT=9
class SeqsManager(SeqFeatComputer):
  '''
  Extends SeqFeatComputer class. Extracts sequences from pdbFiles to fasta files. Then it allows to easily enumerate
  sequences (all letters) and fasta files and also allows for mapping between seqIndices and structIndices and 
  vice versa
  '''
  def __init__(self, rFname, lFname, computedFeatsRootDir= None):
    '''
      @param rFname: str. path to receptor pdb or fasta file
      @param lFname: str. path to ligand pdb or fasta file
      @param computedFeatsRootDir: str. path where features will be stored. If None they will be stored
                                        at default path (assigned in ../Config.py)
    '''

    SeqFeatComputer.__init__(self, rFname, lFname, computedFeatsRootDir)
    self.outPath=        myMakeDir(self.computedFeatsRootDir, "extractedSeqs")
    self.fastaOutDir=    myMakeDir(self.outPath,"seqsData")
    self.seqToStructDir= myMakeDir(self.outPath,"seqToStructMap")
    
    self.parser= PDBParser( QUIET=True)
    self.seqsDict= {}
    self.seqToStruct= {}
    self.structToSeq= {}
    self.seqToStructFnames= {}

  def computeOneFile(self, fileName, chainType):
    '''
      Gets the seq to struct mapping for a given pdb file
      @param fileName: str. fname to pdb file
      @param chainType: str. "l" for ligand and "r" for receptor
    '''
    if self.checkIfIsFasta(fileName):
      self.computeOneFileFromFasta( fileName, chainType)
    else:
      self.computeOneFileFromPDB( fileName, chainType)


  def addResiduesToSeqToStructMap(self, chainType, chainId, seqStr, resIds):
  
    assert len(seqStr)== len(resIds)
    fastaFname= self.seqsDict[chainType][chainId][1]
    self.seqsDict[chainType][chainId]= (seqStr, fastaFname)
    f=open(fastaFname,"w")
    f.write(">"+ os.path.split(fastaFname)[-1]+"\n"+seqStr)
    f.close()
#    print(self.seqToStruct)
#    raw_input("press enter to continue")
    for key in sorted(self.seqToStruct):
      if key[:2]== (chainType, chainId):
        del self.seqToStruct[key]
    listForFile=[]
    for i, resId in enumerate(resIds):
      key_seqStruct=  (chainType, chainId, i)
      flag=" "
      if resId=="-" : continue
      if not resId[-1].isdigit():
        flag= resId[-1]
        resId= resId[:-1]
      else:
        resId= int(resId)
      self.seqToStruct[key_seqStruct]= (" ",resId,flag)
      key_structSeq=  (chainType, chainId, (" ",resId,flag))      
      if key_structSeq in self.structToSeq:
        self.structToSeq[key_structSeq]= i
        listForFile.append( "%d;%s;%s"%(i,seqStr, str((" ",resId,flag))) )

    outName, prefixAndChainType= self.seqToStructFnames[(chainType, chainId)]        
    f=open(outName,"w")
    f.write(">"+ prefixAndChainType+"_"+chainId+"\n"+"\n".join(listForFile))
    f.close()

#    print(self.seqToStruct)
#    raw_input("press enter to continue")
    
  def computeOneFileFromPDB(self, fileName, chainType):
    '''
      Gets the seq to struct mapping for a given pdb file
      @param fileName: str. fname to pdb file
      @param chainType: str. "l" for ligand and "r" for receptor
    '''
    self.seqsDict[chainType]={}
    
    if not ( fileName.endswith("_r_u.pdb") or fileName.endswith("_l_u.pdb")):
      prefixAndChainType = (os.path.split(fileName)[-1]).split(".pdb")[0] +"_"+chainType
    else:
      prefixAndChainType = (os.path.split(fileName)[-1]).split("_u.pdb")[0]
##    print(fileName)
    struct= self.parser.get_structure( prefixAndChainType, fileName)
    for chain in struct[0]:
      chainId= chain.get_id()
      if chainId==" ":
        chainId="*"
      nResStandard= sum([ 1 for res in chain if is_aa(res,standard=True)])
      resList= [ res for res in sorted(chain.child_list, key=lambda x: x.get_id()[1:]) if is_aa(res, standard=False)] #New version feature
      nResAll= len(resList)
#      print(chainId, len(resList))
      if nResStandard<= int(0.5 * nResAll): continue #skip if most residues are not standard
      if len(resList)> SMALL_CHAINS_LIMIT: #Too small chains will not be considered
        sequence=[]
        resIds=[]
        for i,res in enumerate(resList):
          try:
            letter= three_to_one(res.resname)
          except KeyError: # New version feature
            print("Exception", res)
            letter= "X"
            if i== (nResAll-1): break #This case is for TCGR....TLRX where X is GDP or other molecule 
          resId= res.get_full_id()[3]
          sequence.append(letter)
##          print(sequence[-1])
          resIds.append( "%d;%s;%s"%(i,letter,resId) )
          self.seqToStruct[(chainType,chainId, i)]= resId
          self.structToSeq[(chainType,chainId, resId)]= i
        sequence= "".join(sequence)
        outNameFasta= os.path.join(self.fastaOutDir,prefixAndChainType+"_"+chainId+"_u.fasta")
        f=open(outNameFasta,"w")
        f.write(">"+ prefixAndChainType+"_"+chainId+"\n"+sequence)
        f.close()

        resIds= "\n".join(resIds)

        outName= os.path.join(self.seqToStructDir,prefixAndChainType+"_"+chainId+"_u.seqStruMap")
        self.seqToStructFnames[(chainType, chainId)]= (outName, prefixAndChainType)        
        f=open(outName,"w")
        f.write(">"+ prefixAndChainType+"_"+chainId+"\n"+resIds)
        f.close()

        self.seqsDict[chainType][chainId]= (sequence,outNameFasta)

  def computeOneFileFromFasta(self, fileName, chainType):
    '''
      Gets the seq to struct mapping for a given fasta file (dummy, used for compatibility)
      @param fileName: str. fname to fasta file
      @param chainType: str. "l" for ligand and "r" for receptor
    '''
    self.seqsDict[chainType]={}
    
    if not ( fileName.endswith("_r_u.fasta") or fileName.endswith("_l_u.fasta")):
      prefixAndChainType = (os.path.split(fileName)[-1]).split(".pdb")[0] +"_"+chainType
    else:
      prefixAndChainType = (os.path.split(fileName)[-1]).split("_u.fasta")[0]
#    print(fileName,prefixAndChainType, chainType)

    seq= self.parseFasta( fileName)
    chainId=None
    if chainType=="l":
      chainId="L"
    elif chainType=="r":
      chainId="R"
    else:
      raise FeatureComputerException("Error, bad chainType %s for computeOneFileFromFasta, must be 'r' or 'l'"%chainType)

    if len(seq)> SMALL_CHAINS_LIMIT: #Too small chains will not be considered
      sequence=[]
      resIds=[]
      for i, resname in enumerate(seq):
        if not resname in d1_to_index:
          resname= "X"
        resId= (' ', i, ' ')
        sequence.append(resname)
##          print(sequence[-1])
        resIds.append( "%d;%s;%s"%(i,resname,resId) )
        self.seqToStruct[(chainType,chainId, i)]= resId
        self.structToSeq[(chainType,chainId, resId)]= i
      sequence= "".join(sequence)
      outNameFasta= os.path.join(self.fastaOutDir,prefixAndChainType+"_"+chainId+"_u.fasta")
      if not os.path.isfile(outNameFasta):
        f=open(outNameFasta,"w")
        f.write(">"+ prefixAndChainType+"_"+chainId+"\n"+sequence)
        f.close()

      resIds= "\n".join(resIds)

      outName= os.path.join(self.seqToStructDir,prefixAndChainType+"_"+chainId+"_u.seqStruMap")
      if not os.path.isfile(outName):
        f=open(outName,"w")
        f.write(">"+ prefixAndChainType+"_"+chainId+"\n"+resIds)
        f.close()

      self.seqsDict[chainType][chainId]= (sequence,outNameFasta)
    else:
      raise FeatureComputerException("Error, %s is to short (10 AA min) "%prefixAndChainType)
    
  def getSeq(self,chainType, chainId):
    '''
      gets the desired seq of a pdb complex that matches chainType and chainId
      @param chainType: str. "l" for ligand and "r" for receptor
      @param chainId: str. chain id of sequence to be extracted
      @return  (seqStr:str, fastaFileName:str). Tuple. 1st element sequence as str and second element
               path to a fasta file where sequence was extracted
    '''
    return self.seqsDict[chainType][chainId]
  
  def enumSeqs(self, chainType):
    '''
      yields all the sequences contained at pdb file.
      @param chainType: str. "l" for ligand and "r" for receptor
      @yields chainType: str chainId:str, (seqStr:str, fastaFileName:str)  
                            chainType and chain id of sequence to be extracted. 

    '''
    for chainId in self.seqsDict[chainType]:
      yield chainType, chainId

  def getSeqsOutDir(self):
    '''
      returns the path where fasta files are saved for each of the chains of a pdb file
      @return fastaOutDir:str. Path to fasta file
    '''
    return self.fastaOutDir

  def getSeqsMapperOutDir(self):
    '''
      returns the path where seq to struct maps have been saved (No needed)
      @return seqToStructDir:str. Path to seq to struct map
    '''
    return self.seqToStructDir

  def seqToStructIndex(self, chainType, chainId, seqIndex, asString= False):
    '''
      gets the struct id that matches to the chainType, chainId, seqIndex asked
      @param chainType: str. "l" for ligand and "r" for receptor
      @param chainId: str. chain id of sequence
      @param seqIndex: int. Position of the residue at the sequence
      @param asString: boolean. If False, the returned value will be a tuple provided by Bio.PDB.Residue.get_full_id()[3]
                                If True it will be a string obtained by concatenating the tuple and using strip()
      @return None if there is no mapping. Otherwise
                Bio.PDB.Residue.get_full_id()[3] if asString== False
                "".join(Bio.PDB.Residue.get_full_id()[3][1:])).strip() if asString== True
    '''
    try:
#      print(">>", self.seqToStruct[(chainType, chainId, seqIndex)])
#      raw_input("press enter to continue")
      if asString:
        valList= [ str(elem) for elem in self.seqToStruct[(chainType, chainId, seqIndex)]]
        valList= "".join(valList[1:]).strip()
        return valList
      else:
        return self.seqToStruct[(chainType, chainId, seqIndex)]
    except KeyError:
      return None

  def structToSeqIndex(self, chainType,chainId,structIndex):
    '''
      gets the seq index that matches to the chainType, chainId, structIndex asked
      @param chainType: str. "l" for ligand and "r" for receptor
      @param chainId: str. chain id of sequence
      @param structIndex: int. resId as the one provided by Bio.PDB.Residue.get_full_id()[3]
      @return seqIndex: integer. The sequential index of residue with resId==structIndex
    '''
    return  self.structToSeq[(chainType,chainId, structIndex)]
    
def testModuleStruct():
  rFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1FQJ_r_u.pdb"
  lFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1FQJ_l_u.pdb"
  
#  rFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1A22.pdb"
#  lFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1ABR.pdb"
  computedFeatsRootDir="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/tmpComputedFeatures"
  comp= SeqsManager(rFname, lFname, computedFeatsRootDir)
  comp.computeFun()
#  print( comp.structToSeqIndex("l","B",(" ", 22," ")))
#  print( comp.seqToStructIndex("l","B",22))
  raise FeatureComputerException("\nDone: Debug Mode. Comment out testModuleStruct() line 128")
  
def testModuleSeq():
  lFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/seq1_l_u.fasta"
  rFname="/home/rsanchez/Tesis/rriPredMethod/data/develData/inputFasta/seq1_r_u.fasta"
  
#  rFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1A22.pdb"
#  lFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1ABR.pdb"
  computedFeatsRootDir="/home/rsanchez/tmp/computedFeats"
  comp= SeqsManager(rFname, lFname, computedFeatsRootDir)
  comp.computeFun()
#  print( comp.structToSeqIndex("l","B",(" ", 22," ")))
#  print( comp.seqToStructIndex("l","B",22))
  raw_input("test done. Press enter to finish")
  raise FeatureComputerException("\nDone: Debug Mode. Comment out testModuleStruct() line 128")
  
if __name__=="__main__":
#  testModuleStruct()
  testModuleSeq()
  
  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3"
  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures"
  
  SeqFeatComputer.computeAllComplexes(SeqsManager,pdbsIndir= pdbsIndir ,computedFeatsRootDir= computedFeatsRootDir)

