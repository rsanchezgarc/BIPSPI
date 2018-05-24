from __future__ import absolute_import
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np

SCORE_MAT= matlist.blosum62

IDENTITY_THR= 0.95

class HomoOligomerFinder(object):
  '''
  Used to compute resIds mapping from homooligomers proteins
  '''

  def __init__(self, pp_list, interactingPairsResidues, chainType):
    '''
      @param pp_list: [ Bio.PDB.Polypeptide]. list of Bio.PDB.Polypeptide's from either receptor or ligand structure
      @param interactingPairsResidues:  [(ligandResId, receptorResId)]: ligandResId and receptorResIds are full_ids of Bio.PDB.Residue
      @param chainType: str. "l" or "r"
    '''
    self.pp_list= pp_list
    self.interactingPairsResidues= interactingPairsResidues
    self.sequencesList= [ str(onePP.get_sequence()).replace("-","") for onePP in pp_list]
    self.chainIdsList=[onePP[0].get_parent().get_id() for onePP in pp_list]
    self.chainType=  chainType
    
    if chainType=="l":
      self.interactingResidues= {}
      for resL, resR in  self.interactingPairsResidues:
        if not resL in self.interactingResidues:
          self.interactingResidues[resL]=[]
        self.interactingResidues[resL].append(resR)
    elif chainType=="r":
      self.interactingResidues= {}
      for resL, resR in  self.interactingPairsResidues:
        if not resR in self.interactingResidues:
          self.interactingResidues[resR]=[]
        self.interactingResidues[resR].append(resL)
    else:
      raise ValueError("Error, chainType must be 'l' or 'r'")
#    print(self.interactingResidues)
#    raw_input("press enter")
    
  def getIdentityPercent(self, ali0, ali1):
    '''
      @param ali0: str: sequence 0 after alignment
      @param ali1: str: sequence 1 after alignment
    '''
    nElems=len(ali0)
    matches= sum(( 1.0 if elem0==elem1 else 0.0 for elem0, elem1 in zip(ali0, ali1) ))
    return matches/nElems

  def update_interactions(self):
    '''
      fill in self.equivalentInteractingResidues  {res_id_interacting -->res_id_equivalent}
      To do so first computes _all_against_all_ali() and them matches chains that
      have the best structural aligment score. After that, residues of matches chains 
      are mapped thanks to a sequence alignment 
    '''  

    newInteractingPairs=set([])
    for pp0, seq0, chain0 in zip(self.pp_list, self.sequencesList, self.chainIdsList):
      for pp1, seq1, chain1 in zip(self.pp_list, self.sequencesList, self.chainIdsList):
        if chain0== chain1: continue
        alignments = pairwise2.align.localds(seq0, seq1, SCORE_MAT, -11, -0.5)
        ali0= alignments[0][0]
        ali1= alignments[0][1]
        identity= self.getIdentityPercent(ali0, ali1)
        if identity<= IDENTITY_THR: continue
        index0= -1
        index1=-1    
        for letter1, letter0 in zip(ali1, ali0):
          if letter1 !="-":
            index1+=1
          if letter0 !="-":
            index0+=1
          if letter1 !="-" and letter0 !="-":
            res0= pp0[index0].get_full_id()
            res1= pp1[index1].get_full_id()
            if res1 in self.interactingResidues:
#              print(res0, res1)
              if self.chainType=="l":
                for partnerRes in self.interactingResidues[res1]:
                  newInteractingPairs.add((res0, partnerRes ))
              else:
                for partnerRes in self.interactingResidues[res1]:
                  newInteractingPairs.add((partnerRes, res0 ))
#    print( len(newInteractingPairs), len(self.interactingPairsResidues))
    updatedPairs= self.interactingPairsResidues.union(newInteractingPairs)
    if self.chainType=="l":
      involvedResidues= zip(* self.interactingPairsResidues)[0]
    else:
      involvedResidues= zip(* self.interactingPairsResidues)[1]

    chainsInContact= set( [full_res_id[2] for full_res_id in involvedResidues])
#    print(len(updatedPairs))
#    raw_input("enter for next structure")
    return updatedPairs, chainsInContact

