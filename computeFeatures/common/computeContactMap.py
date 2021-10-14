from __future__ import absolute_import
import os
import numpy as np
from Bio.PDB.NeighborSearch import NeighborSearch
from Bio.PDB.Polypeptide import is_aa, PPBuilder, CaPPBuilder

from joblib import load, dump

from pythonTools.myPDBParser import myPDBParser
from ..toolManagerGeneric import ToolManager
from .boundUnboundMapper import BoundUnboundMapper
from .homoOligomerFixer import HomoOligomerFinderWithinPartner, HomoOligomerFinderLR
from .computeVoronoi import getVoroNeigs

from utils import myMakeDir

from bigExceptions import NoValidPDBFile, BadNumberOfResidues

JUST_INTERACTING_CHAINS=False
NEIGS_MAX_DIST_FACTOR= 5.0 #max dist for neigs is res2res_dist*NEIGS_MAX_DIST_FACTOR
class ContactMapper(ToolManager):
  '''
    Extends FeaturesComputer class. Extracts res and chainIds for training and predicting and computes contact maps 
    for a given complex. 
  '''
  VAR_LIST= ['categ']
  
  def __init__(self, prefix, computedFeatsRootDir, areForTrainAndTest=True, boundAvailable= True, 
                     res2res_dist=6.0, statusManager= None):
    '''
      :param prefix: str. An id for the complex
      :param computedFeatsRootDir: str. path where features will be stored

      :param areForTrainAndTest: boolean. True if ligand and receptor are posed in interacting coordinates and thus,
                                we know the label. False if they are for prediction and thus, we cannot assign labels.
      :param: boundAvailable. True if there is a bound and unbound pdb for each complex. False otherwise
      :param res2res_dist: float. max distance between any heavy atoms of 2 amino acids to be considered as interacting
                                  (Angstrom)
      :param statusManager: class that implements .setStatus(msg) to communicate
    '''
    ToolManager.__init__(self, computedFeatsRootDir, statusManager=statusManager)
    self.prefix= prefix
    self.areForTrainAndTest= areForTrainAndTest
    self.res2res_dist= res2res_dist
    self.boundAvailable= boundAvailable
    assert not ( boundAvailable==True and areForTrainAndTest==False), "Error parameters in CMap: boundAvailable==True and areForTrainAndTest==False"
    self.ppb=PPBuilder( radius= 200) # radius set to 200 to not worry about broken chains
    self.outPath= myMakeDir( computedFeatsRootDir, "common")
    self.outPathCM= myMakeDir( self.outPath, "contactMaps")
    self.outPathResDict= myMakeDir( self.outPath, "includedResidues")
    self.outPathNeigs= myMakeDir( self.outPath, "voroNeigs")


  def build_peptides(self, structure):
    pp_list= self.ppb.build_peptides(structure, aa_only= False)      
    if len(pp_list)==0: #case of failure
      pp_list= CaPPBuilder().build_peptides(structure, aa_only= False)
    return pp_list
    
  def mapBoundToUnbound(self, structureUnbound, structureBound, skipBoundChainsIds=set([])):
    '''
      Obtains correspondence between unbound structure and bound structure when available. Returns a dictionary
      that maps bound_residue --> equivalent unbound_residue
      
      :param structureUnbound: Bio.PDB.Structure. Structure in bound state
      :param structureBound:   Bio.PDB.Structure. Structure in unbound state
      :param skipBoundChainsIds:   Set of Chars. Set of chain ids that will be skipped for calculations.
      :return bound2UnboundMapDict: Dict {Bio.PDB.Residue (from bound structure): Bio.PDB.Residue (from unbound structure)}
      :return setOfUnboundResiduesMapped: {Bio.PDB.Residue (from unbound structure) }
      
    '''
    bound2UnboundMapDict={}
    pp_list_unbound= self.build_peptides(structureUnbound)    
    if structureBound is None: # if there is no bound structure, use just unbound.
      boundToUnboundMap= lambda x: x  #For a given residue will return the same residue
      pp_list_bound= pp_list_unbound
    else:
      pp_list_bound= self.build_peptides(structureBound)    
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
    return bound2UnboundMapDict, set(bound2UnboundMapDict.values())
        
  def fixHomoOligomers(self, structureL, structureR, positiveContacts, chainsNotContactL, chainsNotContactR, fixLRHomo=False, fixWithinPartner=False):
    '''
    Option fixLRHomo:
      For each interacting pair of residues (resL_1, resR_2), it will add to positiveContacts (resR_2', res_1', resR_2)

    Option fixWithinPartner:
      For each interacting pair of residues (resL_1, resR_2), it will add to positiveContacts (res_1L', resR_2) and/or
      (resL_1, resR_2') where resL_1' is an equivalent residue in homooligomers of ligand
      
      :param structureL: Bio.PDB.Structure. Structure of ligand
      :param structureR:   Bio.PDB.Structure. Structure of receptor
      :param positiveContacts:  [(ResidueL, ResidueR)]: ligandResId and receptorResIds are full_ids of Bio.PDB.Residue
      :param chainsNoContactL:  Ligand chains that are not in contact with Recetpor
      :param chainsNoContactR:  Recetpor chains that are not in contact with Ligand
      :return positiveContacts, chainsNoContactL, chainsNoContactR. Updated with equivalent residues interactions added
      
    '''
    positiveContactsOri= positiveContacts.copy()
    pp_list_l= self.build_peptides(structureL)
    pp_list_r= self.build_peptides(structureR)

    #Fix for homo-interactions
    if fixLRHomo:
      equivalentLmapper= HomoOligomerFinderLR(pp_list_l, pp_list_r, positiveContacts)
      positiveContacts, (chainsInContactL, chainsInContactR) = equivalentLmapper.update_interactions()
      chainsNotContactL= chainsNotContactL.difference(chainsInContactL)
      chainsNotContactR= chainsNotContactR.difference(chainsInContactR)

    #Fix within partner
    if fixWithinPartner:
      equivalentRmapper= HomoOligomerFinderWithinPartner(pp_list_r, positiveContacts, chainType="r")
      positiveContacts, (chainsInContactL, chainsInContactR)= equivalentRmapper.update_interactions()
      chainsNotContactR= chainsNotContactR.difference(chainsInContactR)

      equivalentLmapper= HomoOligomerFinderWithinPartner(pp_list_l, positiveContacts, chainType="l")
      positiveContacts, (chainsInContactL, chainsInContactR)= equivalentLmapper.update_interactions()
      chainsNotContactL= chainsNotContactL.difference(chainsInContactL)

    homoCorrectionSuccess= positiveContacts != positiveContactsOri
    print("Will contact map be corrected for homo-oligomers? %s"% homoCorrectionSuccess)
    if homoCorrectionSuccess:
      print("%d new pairs added (from %d) "%(len(positiveContacts) - len(positiveContactsOri), len(positiveContactsOri) ))
    return positiveContacts, chainsNotContactL, chainsNotContactR
    
  def getPairsOfResiduesInContact(self, structureL, structureR):
    '''
      Computes which amino acids of ligand are in contact with which amino acids of receptor
      
      :param structureL: Bio.PDB.Structure. Structure of ligand unbound state if available
      :param structureR:   Bio.PDB.Structure. Structure of receptor unbound state if available.
      :return positiveContacts, chainsNotContactL, chainsNotContactR
      
               positiveContacts:  Set {( Bio.PDB.Residue.fullResId (from bound structure structureL), 
                                        Bio.PDB.Residue.fullResId (from bound structure structureR)  )
                                      }
              chainsNotContactL: Set { Bio.PDB.Chain.get_id()}  for ligand chains that are not in contact
              chainsNotContactR: Set { Bio.PDB.Chain.get_id()}  for receptor chains that are not in contact              
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
    positiveContactsResidues= set([])
    chainsInContactL= set([])
    chainsInContactR= set([])
    for res1,res2 in allNeigs:
      pdbId1, modelId1, chainId1, resId1 = res1.get_full_id()
      pdbId2, modelId2, chainId2, resId2 = res2.get_full_id()
      if pdbId1 == lStructId and pdbId2 == rStructId:
        positiveContactsResidues.add( (res1, res2))
        chainsInContactL.add(chainId1)
        chainsInContactR.add(chainId2)
      elif pdbId1 == rStructId and pdbId2 == lStructId:
        positiveContactsResidues.add( (res2, res1))
        chainsInContactL.add(chainId2)
        chainsInContactR.add(chainId1)

    allChainsL= set([ elem.get_id() for elem in  structureL[0].get_list()])
    allChainsR= set([ elem.get_id() for elem in  structureR[0].get_list()])
    chainsNotContactL= allChainsL.difference(chainsInContactL)
    chainsNotContactR= allChainsR.difference(chainsInContactR)        
    return positiveContactsResidues, chainsNotContactL, chainsNotContactR

  def getFNames(self, prefix):
    '''
    Returns a list that contains the fnames that will be used by computeContactMap
    :param prefix: str. prefix for output fnames. Example 1A2K_.
    :return [fname1,fnam2...]
    '''
    outNameCM= os.path.join(self.outPathCM, prefix+".cMap.tab.gz")
    outNameResDict= os.path.join(self.outPathResDict, prefix+".resDict.joblib.pkl.gz")
    outNameNeigs= os.path.join(self.outPathNeigs, prefix+".neigs.joblib.pkl.gz")
    
    return [outNameCM, outNameResDict, outNameNeigs ]
    
  def prepareDataForCmapTrain(self, fnameL, fnameR, structureL, structureR):
    structureL_u, fnameL_u= structureL, fnameL
    structureR_u, fnameR_u= structureR, fnameR      
    if self.boundAvailable==False or not self.areForTrainAndTest:
      structureL_b= None
      structureR_b= None
    else:
      myParserObj= myPDBParser(QUIET=True)
      try:
        lStructId_b = self.prefix+"_l_b.pdb"
        rStructId_b = self.prefix+"_r_b.pdb"        
        fnameL_b= os.path.join( os.path.split(fnameL_u)[0], lStructId_b)
        fnameR_b= os.path.join( os.path.split(fnameR_u)[0], rStructId_b)
        structureL_b =  myParserObj.get_structure( fnameL_b)
        structureR_b =  myParserObj.get_structure( fnameR_b)
      except IOError as e: # in this case there are just unbound pdbs available
        structureL_b= None
        structureR_b= None

    if structureL_b is None or structureR_b is None: #Compute contacts in unbound structures
      positiveContactsRes, chainsNotContactL, chainsNotContactR= self.getPairsOfResiduesInContact(structureL_u, structureR_u)
      lUnb_reslist= structureL_u[0].get_residues()
      rUnb_reslist= structureR_u[0].get_residues()
    else: #Compute contacts in bound structures
      positiveContactsRes, chainsNotContactL, chainsNotContactR= self.getPairsOfResiduesInContact(structureL_b, structureR_b)
      lResBound2Unbo, lUnb_reslist= self.mapBoundToUnbound(structureL_u, structureL_b, skipBoundChainsIds=chainsNotContactL)
      rResBound2Unbo, rUnb_reslist= self.mapBoundToUnbound(structureR_u, structureR_b, skipBoundChainsIds=chainsNotContactR)
  
      positiveContactsRes= set([ (lResBound2Unbo[resL], rResBound2Unbo[resR]) for resL, resR in positiveContactsRes
                                    if resL in lResBound2Unbo and resR in rResBound2Unbo]) # contacting residues mapped to unbound
#    print( [ (key.get_full_id()[2:], key.resname, "-->", rResBound2Unbo[key].get_full_id()[2:], rResBound2Unbo[key].resname) 
#                              for key in sorted(rResBound2Unbo)] ); raw_input("enter")    
    if JUST_INTERACTING_CHAINS==False:
      chainsNotContactR=set([])
      chainsNotContactL=set([])

    residuesL_u= [res for res in structureL_u.get_residues() if is_aa(res, standard= self.filterOutNoStandard) and res in lUnb_reslist ]
    residuesR_u= [res for res in structureR_u.get_residues() if is_aa(res, standard= self.filterOutNoStandard) and res in rUnb_reslist ]
    return  residuesL_u, residuesR_u, positiveContactsRes, chainsNotContactL, chainsNotContactR
    
  def prepareDataForCmapPred(self, fnameL, fnameR, structureL, structureR):

    isSeqL= structureL is None
    isHomo= fnameR is None
    isSeqR= not isHomo and structureR is None

    if isSeqL:
      structureL= self.getStructFromFasta(fnameL, "l") #create fake pdb from fasta file

    if isSeqR:
      structureR = self.getStructFromFasta(fnameR, "r") #create fake pdb from fasta file
    elif isHomo:
      structureR = structureL.copy()
      if isSeqL: #if seq, chains  are named L and R
        structureR[0]["L"].id = "R"
        structureR._reset_full_id()

    # if structureR is None:#then R is a hetero-sequence or is homo
    #   if fnameR: #then is a hetero_sequence, create fake pdb
    #     structureR= self.getStructFromFasta(fnameR, "r")
    #   else: # then is homo, copy the pdb
    #     structureR= structureL.copy()
    #     if fnameR is None and isSeqL:#in the case of homo seq rename of chain receptor is needed
    #       structureR[0]["L"].id= "R"
    #       structureR._reset_full_id()

    residuesL= structureL.get_residues()
    residuesR= structureR.get_residues()
      
    positiveContactsRes= None
    chainsNotContactR=set([])
    chainsNotContactL=set([])
      
    residuesL= [res for res in residuesL if is_aa(res, standard= self.filterOutNoStandard) ]
    residuesR= [res for res in residuesR if is_aa(res, standard= self.filterOutNoStandard) ]
    return  residuesL, residuesR, positiveContactsRes, chainsNotContactL, chainsNotContactR    

        
  def computeComplex(self, fnameL, fnameR, structureL, structureR, isHomoLR=False):
    '''
      Computes the contact map of a complex. Initial input for complex codification. Contact map is a file written at
      self.computedFeatsRootDir/common/contactMaps/ with name perfix_.cMap.tab where prefix is the common name 
      of ligand and receptor pdb files
      1A2K_l_.pdb and 1A2K_r_.pdb   --> 1A2K_.cMap.tab
      
      Computes the neighbours of each residue using voronoi diagrams and stores them in 
      :param fnameL: str. fname to pdbfile or fname to fasta of ligand
      :param fnameR: str. fname to pdbfile or fname to fasta of receptor.
      :param structureL: Bio.PDB.Structure.Structure. Structure of ligand protein (contained in fnameL). None if fasta
      :param structureR: Bio.PDB.Structure.Structure. Structure of receptor protein (contained in fnameR). None if fasta
      :param isHomoLR. True is l and r structures are copies of the same protein.
      :return  resDictL, resDictR  Mappings from Bio.PDB.Residue.Residue to resIds in my format and viceversa
          resDictL =  {resIdL: resObjL, ... resObjL: resIdL}
          resDictR =  {resIdR: resObjR, ... resObjR: resIdR}          
    '''
    prefixToWrite= self.makePrefixExtendedPairwise( self.prefix)
    outNameCM, outNameResDict, outNameNeigs= self.getFNames(prefixToWrite)
    if self.checkAlreayComputed( prefixToWrite):
      dictResDict= load( outNameResDict)
      print ('Already computed contact map %s'%(prefixToWrite))
      return  dictResDict["resDictL"], dictResDict["resDictR"]
    
    print("Computing cmap for %s"%outNameCM)

    if self.areForTrainAndTest:
      residuesL_u, residuesR_u, positiveContactsRes, chainsNotContactL, chainsNotContactR= \
                                            self.prepareDataForCmapTrain(fnameL, fnameR, structureL, structureR)
    else:
      residuesL_u, residuesR_u, positiveContactsRes, chainsNotContactL, chainsNotContactR= \
                                            self.prepareDataForCmapPred(fnameL, fnameR, structureL, structureR)

    if isHomoLR:
      positiveContactsRes, chainsNotContactL, chainsNotContactR= self.fixHomoOligomers(structureL, structureR, positiveContactsRes, chainsNotContactL, chainsNotContactR, fixLRHomo=True)

    nResiduesL= len(residuesL_u)
    nResiduesR= len(residuesR_u)
    if not (self.minNumResiduesPartner< nResiduesL < self.maxNumResiduesPartner):
      raise BadNumberOfResidues(nResiduesL, "1") #1 for ligand
    if not (self.minNumResiduesPartner< nResiduesR < self.maxNumResiduesPartner):
      raise BadNumberOfResidues(nResiduesR, "2") #2 for receptor
      
    fromRes2ChainResIdAndNameF= self.fromRes2ChainResIdAndName
    
    idsTuplesL= [ fromRes2ChainResIdAndNameF(res) for res in residuesL_u]
    residues_idsTuplesL= [ (res, idTuple) for res, idTuple in zip(residuesL_u, idsTuplesL) if 
                            idTuple!=None and idTuple[-1]!="X"]

    idsTuplesR= [ fromRes2ChainResIdAndNameF(res) for res in residuesR_u]
    residues_idsTuplesR= [ (res, idTuple) for res, idTuple in zip(residuesR_u, idsTuplesR) if 
                            idTuple!=None and idTuple[-1]!="X"]
                            
    if positiveContactsRes is None:
      areInContactF= lambda resL, resR: np.nan
    else:
      areInContactF= lambda resL, resR: 1 if (resL, resR) in positiveContactsRes else -1 
      
    try:
      listOfRows=[]
      for resL, elementsOfIdL in residues_idsTuplesL:
        chainIdL, resIdL, resNameL = elementsOfIdL
        for resR, elementsOfIdR in residues_idsTuplesR:
          elementsOfIdR= fromRes2ChainResIdAndNameF(resR)
          chainIdR, resIdR, resNameR = elementsOfIdR
          areInContact=  areInContactF(resL, resR)
          listOfRows.append("%s %s %s %s %s %s %s" %( chainIdL, resIdL, resNameL, chainIdR, 
                                                         resIdR, resNameR, areInContact))
  
      self.writeResultsFromDataDictPairL2R( listOfRows, outNameCM)

      if not structureL is None:
        neigsDictL= getVoroNeigs(residuesL_u, self.res2res_dist*NEIGS_MAX_DIST_FACTOR)
        neigsDictL= { fromRes2ChainResIdAndNameF(resObj)[:-1] : set([fromRes2ChainResIdAndNameF(neig)[:-1] for neig 
                          in neigsDictL[resObj] ]) for resObj in neigsDictL }
      else:
        neigsDictL= None
      if not structureR is None:
        neigsDictR= getVoroNeigs(residuesR_u, self.res2res_dist*NEIGS_MAX_DIST_FACTOR)
        neigsDictR= { fromRes2ChainResIdAndNameF(resObj)[:-1] : set([fromRes2ChainResIdAndNameF(neig)[:-1] for neig 
                          in neigsDictR[resObj] ]) for resObj in neigsDictR }
      else:
        neigsDictR= None
        
      dump( {"neigsDictL":neigsDictL, "neigsDictR":neigsDictR}, outNameNeigs, compress=3)
  
      resDictL=  {}
      for resL, elementsOfIdL in residues_idsTuplesL:
        resDictL[resL]= elementsOfIdL
        resDictL[elementsOfIdL]= resL
      resDictR=  {}
      for resR, elementsOfIdR in residues_idsTuplesR:
        resDictR[resR]= elementsOfIdR
        resDictR[elementsOfIdR]= resR
        
      dump( {"resDictL":resDictL, "resDictR":resDictR}, outNameResDict, compress=3)
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%outNameCM)
      self.tryToRemoveAllFnames( prefixToWrite)    
      raise
    return  resDictL, resDictR


if __name__=="__main__":
  import sys

  fnameL= sys.argv[1]
  fnameR= sys.argv[2]
  computedFeatsRootDir = sys.argv[3]
  prefix= os.path.basename(fnameL).split(".")[0].split("_")[0]
  isHomoLR = False

  parser=myPDBParser()
  structureL= parser.get_structure(fnameL)
  structureR = parser.get_structure(fnameR)

  cmapper= ContactMapper(prefix, computedFeatsRootDir, areForTrainAndTest=True, boundAvailable= False)
  cmapper.computeComplex(fnameL, fnameR, structureL, structureR, isHomoLR=isHomoLR)

  '''
python -m computeFeatures.common.computeContactMap ~/tmp/BIPSPI/cmapTrial/pdbs/3hzr-BA_l_u.pdb  ~/tmp/BIPSPI/cmapTrial/pdbs/3hzr-BA_r_u.pdb  ~/tmp/BIPSPI/cmapTrial/wdir/

python -m computeFeatures.common.computeContactMap ../../data/develData/pdbFiles/1A2K_l_u.pdb ../../data/develData/pdbFiles/1A2K_r_u.pdb  ~/tmp/BIPSPI/cmapTrial/wdir/
  '''