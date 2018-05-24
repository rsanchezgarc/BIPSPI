from __future__ import absolute_import
import sys, os
from Bio.PDB.NeighborSearch import NeighborSearch
from computeFeatures.structStep.myPDBParser import myPDBParser as PDBParser
from Bio.PDB.Polypeptide import is_aa, three_to_one, PPBuilder, CaPPBuilder
import numpy as np

from ..FeaturesComputer import FeaturesComputer, FeatureComputerException
from .boundUnboundMapper import BoundUnboundMapper
from .fixHomoOligomers import HomoOligomerFinder
from utils import myMakeDir, myMakeDirUnique, tryToRemove #utils is at the root of the package
from bigExceptions import NoValidPDBFile, BadNumberOfResidues


JUST_INTERACTING_CHAINS= False #True to ignore chains that are not involved in interaction
CONSIDER_HOMOOLIG_AS_POS= False #True to consider that if a residue is interacting in one chain it is also interacting in all other homooligomers

class ContactMapper(FeaturesComputer):
  '''
    Extends FeaturesComputer class. Extracts res and chainIds for training and predicting and computes contact maps 
    for training for a given complex
  '''
  def __init__(self, rFname, lFname, computedFeatsRootDir= None, boundAvailable= True, res2res_dist= 6.0,
               isForPrediction=False, statusManager=None):
    '''
      @param rFname: str. path to receptor pdb file
      @param lFname: str. path to ligand pdb file      
      @param computedFeatsRootDir: str. path where features will be stored
      @param boundAvailable: bool. True if bound structures are available. False otherwise. Bound structures must be located
                                   at the same path that unbound structures and need to be named as in the following example:
                                    1A2K_l_u.pdb  1A2K_r_b.pdb
      @param res2res_dist: float. max distance between any heavy atoms of 2 amino acids to be considered as interacting
                                  (Amstrongs)
      @param isForPrediction: bool. False to compute contacts between amino acids, True otherwise. Positive contacts will
                                    be tag as 1, negative as -1. If True, all amino acids will have as tag np.nan
      @param statusManager: class that implements .setStatus(msg) to communicate
    '''
    FeaturesComputer.__init__(self, rFname, lFname, computedFeatsRootDir)

    self.prefixR= os.path.split(rFname)[1].split(".")[0].split("_")[0]
    self.prefixL= os.path.split(lFname)[1].split(".")[0].split("_")[0] 
    if self.prefixR == self.prefixL:
      self.prefix= self.prefixR
    else:
      if "<" in self.prefixL:
        raise FeatureComputerException("Error. Ligand pdbFile name %s must not contain '<' or '>' character"%lFname)
      if ">" in self.prefixR:
        raise FeatureComputerException("Error. Receptor pdbFile name %s must not contain '<' or'>' character"%rFname)
      self.prefixR= self.getExtendedPrefix(rFname)
      self.prefixL= self.getExtendedPrefix(lFname)      
        
      self.prefix= self.prefixL+"<->"+self.prefixR

    self.isForPrediction= isForPrediction      
    self.res2res_dist= res2res_dist
    self.boundAvailable=boundAvailable
    self.outPath= myMakeDir(self.computedFeatsRootDir, "common/contactMaps")
    self.outName= os.path.join(self.outPath,self.prefix+".cMap.tab")    
    self.parser= PDBParser(QUIET=True)
#    self.ppb=PPBuilder( radius= 200) # To not worry for broken chains
    self.ppb=CaPPBuilder()
    self.computeFun= self.contactMapOneComplex

  def mapBoundToUnbound(self, structureUnbound, structureBound, skipBoundChainsIds=set([])):
    '''
      Obtains correspondence between unbound structure and bound structure when available. Returns a dictionary
      that maps bound_residue --> equivalent unbound_residue
      
      @param structureUnbound: Bio.PDB.Structure. Structure in bound state
      @param structureBound:   Bio.PDB.Structure. Structure in unbound state
      @param skipBoundChainsIds:   Set of Chars. Set of chain ids that will be skipped for calculations. 
      @return bound2UnboundMapDict: Dict {Bio.PDB.Residue (from bound structure): Bio.PDB.Residue (from unbound structure)}
      
    '''
    bound2UnboundMapDict={}
    pp_list_unbound= self.ppb.build_peptides(structureUnbound, aa_only= False)    
    if structureBound is None: # if there is no bound structure, use just unbound.
      boundToUnboundMap= lambda x: x  #For a given residue will return the same residue
      pp_list_bound= pp_list_unbound
    else:
      pp_list_bound= self.ppb.build_peptides(structureBound, aa_only= False)    
      mapper= BoundUnboundMapper( pp_list_unbound,pp_list_bound) # res_bound->res_unbound mapper object
      mapper.build_correspondence()
      boundToUnboundMap = mapper.mapBoundToUnbound  #For a given bound residue will return its unbound equivalent
    for pp in pp_list_bound:
      for resBound in pp:
        chainBound= resBound.get_full_id()[2] # str chainId
        if chainBound in skipBoundChainsIds: continue
        resUnbound= boundToUnboundMap( resBound)
        if not resUnbound is None: #In case there is no equivalent unbound residue for a given bound residue        
          bound2UnboundMapDict[ resBound ]= resUnbound
    return bound2UnboundMapDict

  def fixHomooligomers(self, structureL, structureR, positiveContacts, chainsInContactL, chainsInContactR):
    '''
      For each interacting pair of residues (resL_1, resR_2), it will add to positiveContacts (res_1L', resR_2) and/or
      (resL_1, resR_2') where resL_1' is an equivalent residue in homooligomers of ligand
      
      @param structureL: Bio.PDB.Structure. Structure of ligand
      @param structureR:   Bio.PDB.Structure. Structure of receptor
      @param positiveContacts:  [(ligandResId, receptorResId)]: ligandResId and receptorResIds are full_ids of Bio.PDB.Residue
      @param chainsInContactL:  [(ligandResId)]: ligandResId and receptorResIds are full_ids of Bio.PDB.Residue
      @param chainsInContactR:  [(receptorResId)]: ligandResId and receptorResIds are full_ids of Bio.PDB.Residue
      @return positiveContacts, chainsInContactL, chainsInContactR. Updated with equivalent residues interactions added
      
    '''
    pp_list_l= self.ppb.build_peptides(structureL, aa_only=False)
    equivalentLmapper= HomoOligomerFinder(pp_list_l, positiveContacts, chainType="l")
    positiveContacts, chainsInContactL= equivalentLmapper.update_interactions()
    pp_list_r= self.ppb.build_peptides(structureR, aa_only=False)
    equivalentRmapper= HomoOligomerFinder(pp_list_r, positiveContacts, chainType="r")
    positiveContacts, chainsInContactR= equivalentRmapper.update_interactions()
    return positiveContacts, chainsInContactL, chainsInContactR
    
  def getPairsOfResiduesInContact(self, structureL, structureR):
    '''
      Computes which amino acids of ligand are in contact with which amino acids of receptor
      
      @param structureL: Bio.PDB.Structure. Structure of ligand (bound state if available)
      @param structureR:   Bio.PDB.Structure. Structure of receptor (bound state if available).
      @return positiveContacts:  Set {(Bio.PDB.Residue.fullResId (from bound structure structureL), Bio.PDB.Residue.fullResId (from bound structure structureR))}
      @return chainsNotContactL: Set { str(chainId structureL)}
      @return chainsNotContactR: Set { str(chainId structureR)}
      
    '''
    try:
      atomListL = [atom for atom in structureL.child_list[0].get_atoms() if not atom.name.startswith("H")]
    except IndexError:
      raise NoValidPDBFile("Problems parsing pdbFile 1")
    try:      
      atomListR = [atom for atom in structureR.child_list[0].get_atoms() if not atom.name.startswith("H")]
    except IndexError:
      raise NoValidPDBFile("Problems parsing pdbFile 2")

    searcher= NeighborSearch(atomListL+atomListR)
    allNeigs= searcher.search_all(self.res2res_dist, level= "R")
    lStructId= structureL.get_id()
    rStructId= structureR.get_id()
    positiveContacts = set([])
    chainsInContactL= set([])
    chainsInContactR= set([])
    for res1,res2 in allNeigs:
      pdbId1, modelId1, chainId1, resId1 = res1.get_full_id()
      pdbId2, modelId2, chainId2, resId2 = res2.get_full_id()
      fullResId1= res1.get_full_id()
      fullResId2= res2.get_full_id()
      if pdbId1 == lStructId and pdbId2 == rStructId:
        positiveContacts.add(( fullResId1, fullResId2 ))
        chainsInContactL.add(fullResId1[2])
        chainsInContactR.add(fullResId2[2])       
      elif pdbId1 == rStructId and pdbId2 == lStructId:
        positiveContacts.add(( fullResId2, fullResId1))
        chainsInContactL.add(fullResId2[2])
        chainsInContactR.add(fullResId1[2])
    if CONSIDER_HOMOOLIG_AS_POS:
      positiveContacts, chainsInContactL, chainsInContactR= self.fixHomooligomers(structureL, structureR, 
                                                                 positiveContacts, chainsInContactL, chainsInContactR)
    allChainsL= set([ elem.get_id() for elem in  structureL[0].get_list()])
    allChainsR= set([ elem.get_id() for elem in  structureR[0].get_list()])
    chainsNotContactL= allChainsL.difference(chainsInContactL)
    chainsNotContactR= allChainsR.difference(chainsInContactR)
    return positiveContacts, chainsNotContactL, chainsNotContactR

  def contactMapOneComplex(self):
    '''
      Computes the contact map of a complex. Initial input for complex codification. Contact map is a file written at
      self.computedFeatsRootDir/common/contactMaps/ with name prefix.cMap.tab where prefix is either the common name of
      ligand and receptor pdb files or the concatenation of ligand and receptor names.
      1A2K_l_u.pdb and 1A2K_r_u.pdb  --> 1A2K.cMap.tab
      1A2K_l_u.pdb and 1A22.pdb  --> 1A2K-1A22.cMap.tab
      
    '''    
    outName= self.outName
    print (outName)
    if os.path.isfile(outName):
      print ('Already computed contact map')
      return 0
    lStructId = self.prefixL+"_l_u.pdb"
    rStructId = self.prefixR+"_r_u.pdb"
    structureL_u =  self.parser.get_structure(lStructId, self.lFname)
    structureR_u =  self.parser.get_structure(rStructId, self.rFname)      
    if self.boundAvailable==False or self.isForPrediction:
      structureL_b= None
      structureR_b= None
    else:   
      try:
        lStructId_b = self.prefix+"_l_b.pdb"
        rStructId_b = self.prefix+"_r_b.pdb"        
        lFname_b= os.path.join( os.path.split(self.lFname)[0], lStructId_b)
        rFname_b= os.path.join( os.path.split(self.rFname)[0], rStructId_b)
        structureL_b =  self.parser.get_structure(lStructId_b, lFname_b)
        structureR_b =  self.parser.get_structure(rStructId_b, rFname_b)
      except IOError as e: # in this case there are just unbound pdbs available
        structureL_b= None
        structureR_b= None

    if self.isForPrediction:
      positiveContacts= None
      chainsNotContactR=set([])
      chainsNotContactL=set([])
    elif structureL_b is None or structureR_b is None: #Compute contacs in bound structures
      positiveContacts, chainsNotContactL, chainsNotContactR= self.getPairsOfResiduesInContact(structureL_u, structureR_u)
    else: #Compute contacs in unbound structures
      positiveContacts, chainsNotContactL, chainsNotContactR= self.getPairsOfResiduesInContact(structureL_b, structureR_b)

    if JUST_INTERACTING_CHAINS==False:
      chainsNotContactR=set([])
      chainsNotContactL=set([])     

    rResDict= self.mapBoundToUnbound(structureR_u, structureR_b, skipBoundChainsIds=chainsNotContactR)          
    lResDict= self.mapBoundToUnbound(structureL_u, structureL_b, skipBoundChainsIds=chainsNotContactL)
    nResiduesL= len(lResDict)
    nResiduesR= len(rResDict)
    if not (self.minNumResiduesPartner< nResiduesL < self.maxNumResiduesPartner):
      raise BadNumberOfResidues(nResiduesL, "1")
    if not (self.minNumResiduesPartner< nResiduesR < self.maxNumResiduesPartner):
      raise BadNumberOfResidues(nResiduesR, "2")
    
    outFile= open(outName,"w")
    outFile.write("chainIdL structResIdL resNameL chainIdR structResIdR resNameR categ\n")
#    print(sorted(lResDict, key= lambda x: x.get_id()))
#    a= raw_input()
    try:
      for resL_bound in sorted(lResDict, key= lambda x: x.get_full_id()):
  #      print(resL_bound.get_full_id())
        resL_unbound= lResDict[resL_bound]
        pdbIdL, modelL, chainIdL, resIdL= resL_unbound.get_full_id()
        resIdL= self.makeStrResId(resIdL)

        try:
          letraL= three_to_one(resL_unbound.resname)
          if letraL!= three_to_one(resL_bound.resname): continue
        except KeyError:
            continue  
        for resR_bound in sorted(rResDict, key= lambda x: x.get_full_id()):
          resR_unbound= rResDict[resR_bound]
          pdbIdR, modelR, chainIdR, resIdR = resR_unbound.get_full_id()
          try:
            letraR= three_to_one(resR_unbound.resname)
            if letraR!= three_to_one(resR_bound.resname): continue
          except KeyError:
            continue
          if self.isForPrediction:
            categ= np.nan
          elif (resL_bound.get_full_id(),resR_bound.get_full_id()) in positiveContacts:
            categ= 1
          else:
            categ= -1
          resIdR= self.makeStrResId(resIdR)
          if chainIdL==" ": chainIdL="*"
          if chainIdR==" ": chainIdR="*"
  #        print("%s %s %s %s %s %s %s\n" %(chainIdL, resIdL, letraL, chainIdR, resIdR, letraR, categ) )
  #        raw_input("enter")
          outFile.write("%s %s %s %s %s %s %s\n" %(chainIdL, resIdL, letraL, chainIdR, resIdR, letraR, categ))
      outFile.close()
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%outName)
      tryToRemove(outName)    
      raise

  def makeStrResId(self, resId):
    valList= [ str(elem) for elem in resId[1:]]
    finalId= "".join(valList).strip()
    return finalId

def testModulePredict():
#  rFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_r_u.pdb"
#  lFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_l_u.pdb"
  
  rFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1A22.pdb"
  lFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1ABR.pdb"
  
  computedFeatsRootDir="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/tmpComputedFeatures"
  boundAvailable= False
  res2res_dist= 6.0
  isForPrediction= True
  comput= ContactMapper(rFname, lFname, computedFeatsRootDir, boundAvailable, res2res_dist, isForPrediction)
  comput.computeFun()
  raise FeatureComputerException("Debug Mode. Comment out testModule() line 277")

def testModuleTrain():
  pdbsIndir="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/"
  prefix= "1ACB"
#  prefix= "1N2C"
#  prefix= "1KKL"  #Este pdb es mucha caca porque es un homotrimero
#  prefix= "1N8O"  #Este pdb es mucha caca porque mezcla varias chainIds en una sola cadena
  
#  rFname= os.path.join(pdbsIndir,"%s_r_u.pdb"%prefix)
#  lFname= os.path.join(pdbsIndir,"%s_l_u.pdb"%prefix)

  rFname=  "/home/rsanchez/tmp/1acbT0_r.pdb"
  lFname=  "/home/rsanchez/tmp/1acbT0_l.pdb"

  computedFeatsRootDir="/home/rsanchez/tmp/tmpRRI"
  boundAvailable= True
  res2res_dist= 6.0
  isForPrediction= False
  comput= ContactMapper(rFname, lFname, computedFeatsRootDir, boundAvailable, res2res_dist, isForPrediction)
  comput.computeFun()

#  evaluateWithPymol( os.path.join(pdbsIndir, "%s_r_u.pdb"%prefix), os.path.join(pdbsIndir, "%s_l_u.pdb"%prefix),
#                     computedFeatsRootDir, pdbsIndir)  
  raise FeatureComputerException("Debug Mode. Comment out testModule() line 277")
      
            
if __name__=="__main__":

#  testModulePredict()
  testModuleTrain()
  
  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3"
  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures"

#  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/b5structures"
#  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/computedFeatures/"
  
  
  FeaturesComputer.computeFeaturesAllComplexes(ContactMapper,pdbsIndir= pdbsIndir, computedFeatsRootDir= computedFeatsRootDir)



