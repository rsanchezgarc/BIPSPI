from __future__ import absolute_import
import os
import gzip
from Bio.PDB.Polypeptide import is_aa
from computeFeatures.featuresComputer import FeatureComputerException
from computeFeatures.seqStep.seqToolManager import SeqToolManager
from computeFeatures.seqStep.seqToolManagers.seqExtraction.seqAligner import SeqAligner
from utils import myMakeDir
from ast import literal_eval as make_tuple

class SeqStructMapper(SeqToolManager, SeqAligner):
  
  def __init__(self, computedFeatsRootDir, statusManager=None ):
    '''
      :param computedFeatsRootDir: str. root path where results will be saved
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''   
    SeqToolManager.__init__(self, computedFeatsRootDir, winSize=None, statusManager=statusManager)

    self.outPath=        myMakeDir(self.computedFeatsRootDir, "extractedSeqs")
    self.fastaOutDir=    myMakeDir(self.outPath,"seqsData")
    self.seqToStructDir= myMakeDir(self.outPath,"seqToStructMap")   
    self.seqsDict= {}
    self.seqToStruct= {}
    self.structToSeq= {}
    self.seqToStructFnames= {}
    self.seqToRefSeq= {}
    self.refSeqToSeq= {}    
    self.currentSeq={}
    
  def reset(self):
    self.seqsDict= {}
    self.seqToStruct= {}
    self.structToSeq= {}
    self.seqToStructFnames= {}
    self.seqToRefSeq= {}
    self.refSeqToSeq= {}
    self.currentSeq={}
    
  def setCurrentSeq(self, seq, chainType, chainId):
    '''
    sets a sequence to modify  seqToStructIndex behaviour
    If     self.currentSeq=None. Mapping will be refSeq ->structIndex
    else  currentSeq -> refSeq ->structIndex

    '''
    self.currentSeq[( chainType, chainId) ]= seq
    matchingIndices= self.getMatchingSeqsIndices(seq, self.seqsDict[chainType][chainId][0])
    if not matchingIndices is None:
      self.seqToRefSeq[( chainType, chainId) ]= matchingIndices
      self.refSeqToSeq[( chainType, chainId) ]= { b:a for a,b in matchingIndices.iteritems() }
    
  def computeComplex(self, fnameL, fnameR, structureL, structureR):
    '''
    Applies self.computeOneFile method to both ligand and receptor pdbFiles.
    :param fnameOrStructL: str or Bio.PDB.Structure. fname or structure to pdbfile or fname to fasta of ligand
    :param fnameOrStructR: str or Bio.PDB.Structure. fname or structure to pdbfile or fname to fasta of receptor
    :param structureL: Bio.PDB.Structure.Structure. Structure of ligand protein (contained in fnameL). None if fasta
    :param structureR: Bio.PDB.Structure.Structure. Structure of receptor protein (contained in fnameR). None if fasta
    '''
    self.reset()
    if fnameL:
      self.computeOneFile( fnameL, structureL, "l")
    if fnameR:
      self.computeOneFile( fnameR, structureR, "r")
    
  def computeOneFile(self, fileName, struct, chainType):
    '''
      Gets the seq to struct mapping for a given pdb file
      :param fileName: str. fname to pdb file
      :param struct: Bio.PDB.Structure or None if fileName points to fasta file
      :param chainType: str. "l" for ligand and "r" for receptor
    '''


    extendedPrefix= self.getExtendedPrefix(fileName)
    if struct is None:
      struct= self.getStructFromFasta(fileName, chainType)
      
    self.computeOneFileFromPDB( extendedPrefix, struct, chainType)

  def saveFasta(self, seqStr, fastaFname):
    with open(fastaFname,"w") as f:
      f.write(">"+ os.path.basename(fastaFname)+"\n"+seqStr)
      
  def saveSeqStructMap(self, listOfResIdsSeqIds, seqStruMapFname):
    '''
      listOfResIdsSeqIds= [ (i,aaLetter,resId)  ]
    '''
    with open(seqStruMapFname,"w") as f:
      f.write(">"+ os.path.basename(seqStruMapFname)+"\n")
      for (seqId, resName, resId) in listOfResIdsSeqIds:
        f.write("%d;%s;%s\n"%(seqId, resName, resId) )
    
  def loadSeqStruct(self, seqStruMapFname):
    listResMaps= []
    with open(seqStruMapFname) as f:
      header= f.readline()
      for line in f:
        (seqId, resName, resId_tupleStr) = line.split(";")
        resId= make_tuple(resId_tupleStr)
        seqId= int(seqId)
        resId_tuple= (resId[0], int(resId[1]), resId[2])
        listResMaps.append( (seqId, resName, resId_tuple)  )
    return listResMaps
    
  def addResiduesToSeqToStructMap(self, chainType, chainId, seqStr, resIds):
    '''
      Given an already mapped seq to struct object, modify one chain to potentially add new residues
      Needed if 3dcons is used as it generally report all residues in sequence and not just the included in pdb
    '''
    assert len(seqStr)== len(resIds), "error, there must be the same number of residues as amino acids are in sequence"
    for key in sorted(self.seqToStruct): #Remove first old residues
      if key[:2]== (chainType, chainId):
        del self.seqToStruct[key]
        
    listOfResIdsSeqIds=[]
    for i, resId in enumerate(resIds):
      key_seq2Struct= (chainType, chainId, i)
      flag=" "
      if resId=="-" : continue
      if not resId[-1].isdigit():
        flag= resId[-1]
        resId= resId[:-1]
      else:
        resId= int(resId)
      tuple_resId= (" ",resId,flag)
      self.seqToStruct[key_seq2Struct]= tuple_resId
      listOfResIdsSeqIds.append( (i, seqStr[i], tuple_resId ) )
      key_struct2Seq=  (chainType, chainId, tuple_resId )      
      if key_struct2Seq in self.structToSeq:
        self.structToSeq[key_struct2Seq]= i

    outNameSeqStructMap, prefixAndChainType= self.seqToStructFnames[(chainType, chainId)]
    __, outNameFasta= self.seqsDict[chainType][chainId]
    self.seqsDict[chainType][chainId]= (seqStr, outNameFasta)                  
    self.saveFasta(seqStr, outNameFasta)
    self.saveSeqStructMap( listOfResIdsSeqIds, outNameSeqStructMap)
    
  def computeOneFileFromPDB(self, extendedPrefix, struct, chainType):
    '''
      Gets the seq to struct mapping for a given pdb file
      :param fileNameOrStruct: str. fname to pdb file or Bio.PDB.Structure
      :param chainType: str. "l" for ligand and "r" for receptor
    '''

    self.seqsDict[chainType]={}
    split_extendedPrefix=  self.splitExtendedPrefix(extendedPrefix)
    prefix= split_extendedPrefix[0]
    if len(split_extendedPrefix)==3:
      boundUnbound= split_extendedPrefix[-1]+"_"
    else:
      boundUnbound= ""
    
    for chain in struct[0]:
      chainId= chain.get_id()
      if chainId==" ": chainId= type(self).UNKNOWN_CHAIN_SYMB
      baseNamePrefix= "%s_%s_%s_%s"%(prefix, chainType, chainId, boundUnbound)
      outNameFasta= os.path.join(self.fastaOutDir, baseNamePrefix+".fasta")
      outNameSeqStructMap= os.path.join(self.seqToStructDir, baseNamePrefix+".seqStruMap")
      
      if os.path.isfile(outNameFasta) and os.path.isfile(outNameSeqStructMap):
        print("loading seq to struct map %s"%(outNameSeqStructMap) )
        seqStr= self.parseFasta(outNameFasta, chainType)
        listOfResIdsSeqIds= self.loadSeqStruct(outNameSeqStructMap)
      else:
        nResStandard= sum( ( 1 for res in chain if is_aa(res,standard=True) ) )
        resList= [ res for res in sorted(chain.child_list, key=lambda x: x.get_id()[1:]) if is_aa(res, standard=False)] #New version feature
        nResAll= len(resList)
        if nResStandard<= int(0.5 * nResAll): continue #skip if most residues are not standard
        nBadContinuous=0
        sequence=[]
        listOfResIdsSeqIds=[]
##        print(resList); raw_input("enter")
        for i,res in enumerate(resList):
          __, __, resName= self.fromRes2ChainResIdAndName(res)
          if resName is None: continue
#          print(i, resId, resName)
          resId= res.get_id()
          if resName=="X" and resId[0].strip()!="":
            nBadContinuous+=1
          else:
            nBadContinuous=0
          listOfResIdsSeqIds.append( (i, resName, resId, ) )
          sequence.append(resName)
        if nBadContinuous>0:
          listOfResIdsSeqIds= listOfResIdsSeqIds[:-nBadContinuous]
        seqStr= "".join(sequence)
#        print(seqStr); raw_input("enter")
      if len(seqStr) > self.minNumResiduesPartner: #Too small chains will not be considered      
        self.saveFasta(seqStr, outNameFasta)
        self.saveSeqStructMap( listOfResIdsSeqIds, outNameSeqStructMap)
        for i, resName, resId in listOfResIdsSeqIds:
          self.seqToStruct[(chainType,chainId, i)]= resId
          self.structToSeq[(chainType,chainId, resId)]= i  
        self.seqsDict[chainType][chainId]= (seqStr, outNameFasta)                  
        self.seqToStructFnames[(chainType, chainId)]= (outNameSeqStructMap, baseNamePrefix)


  def getSeq(self,chainType, chainId):
    '''
      gets the desired seq of a pdb complex that matches chainType and chainId
      :param chainType: str. "l" for ligand and "r" for receptor
      :param chainId: str. chain id of sequence to be extracted
      :return  (seqStr:str, fastaFileName:str). Tuple. 1st element sequence as str and second element
               path to a fasta file where sequence was extracted
    '''
    return self.seqsDict[chainType][chainId]
  
  def enumSeqs(self, chainType):
    '''
      yields all the sequences contained at pdb file.
      :param chainType: str. "l" for ligand and "r" for receptor
      @yields chainType: str chainId:str, (seqStr:str, fastaFileName:str)  
                            chainType and chain id of sequence to be extracted. 

    '''
    for chainId in self.seqsDict[chainType]:
      yield chainType, chainId

  def getSeqsOutDir(self):
    '''
      returns the path where fasta files are saved for each of the chains of a pdb file
      :return fastaOutDir:str. Path to fasta file
    '''
    return self.fastaOutDir

  def getSeqsMapperOutDir(self):
    '''
      returns the path where seq to struct maps have been saved (No needed)
      :return seqToStructDir:str. Path to seq to struct map
    '''
    return self.seqToStructDir

  def seqToStructIndex(self, chainType, chainId, seqIndex, asString= False):
    '''
      gets the struct id that matches to the chainType, chainId, seqIndex asked
      :param chainType: str. "l" for ligand and "r" for receptor
      :param chainId: str. chain id of sequence
      :param seqIndex: int. Position of the residue at the sequence
      :param asString: boolean. If False, the returned value will be a tuple provided by Bio.PDB.Residue.get_full_id()[3]
                                If True it will be a string obtained by concatenating the tuple and using strip()
      :return None if there is no mapping. Otherwise
                Bio.PDB.Residue.get_full_id()[3] if asString== False
                "".join(Bio.PDB.Residue.get_full_id()[3][1:])).strip() if asString== True
    '''
    try:
      if (chainType, chainId) in self.seqToRefSeq[(chainType, chainId)]:
        seqIndex= self.seqToRefSeq[(chainType, chainId)][seqIndex]
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
      :param chainType: str. "l" for ligand and "r" for receptor
      :param chainId: str. chain id of sequence
      :param structIndex: int. resId as the one provided by Bio.PDB.Residue.get_full_id()[3]
      :return seqIndex: integer. The sequential index of residue with resId==structIndex
    '''
    try:
      seqIndex=  self.structToSeq[(chainType,chainId, structIndex)]
      if (chainType, chainId) in self.refSeqToSeq[(chainType, chainId)]:
        seqIndex= self.refSeqToSeq[(chainType, chainId)][[seqIndex]]
    except KeyError:
      return None    



