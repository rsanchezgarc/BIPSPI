from __future__ import absolute_import
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from computeFeatures.toolManagerGeneric import threeLetterAA_to_one

SCORE_MAT= matlist.blosum62
IDENTITY_THR= 0.95

class HomoOligomerFinderBase(object):

  def _getIdentityPercent(self, ali0, ali1):
    '''
      :param ali0: str: sequence 0 after alignment
      :param ali1: str: sequence 1 after alignment
    '''
    nElems=len(ali0)
    matches= sum(( 1.0 if elem0==elem1 else 0.0 for elem0, elem1 in zip(ali0, ali1) ))
    return matches/nElems

  def _formatInteractingPairsAsDict(self):
    interactingL={}
    interactingR={}
    for resL, resR in self.interactingPairsResidues:
      if not resL in interactingL:
        interactingL[resL] = []
      interactingL[resL].append(resR)
      if not resR in interactingR:
        interactingR[resR] = []
      interactingR[resR].append(resL)

    return interactingL, interactingR

  def _fromPolypeptideToSequence(self, pp):
    return "".join([ threeLetterAA_to_one(aa.resname) for aa in pp ]).replace("-","")

  def _fromPairsInContactToChains(self, involvedPairs):
    chainsInContactL, chainsInContactR= [], []
    for resL, resR in involvedPairs:
      chainsInContactL.append( resL.get_full_id()[2] )
      chainsInContactR.append(resR.get_full_id()[2])
    return chainsInContactL, chainsInContactR

  def _alignSeqs(self, seq0, seq1):
    return pairwise2.align.localds(seq0, seq1, SCORE_MAT, -11, -0.5)

  def _getChainId(self, pp):
    return pp[0].get_parent().get_id()

  def findEquivalentResidues(self, pp_list0, pp_list1, excludeSameChainId=False):
    eqivalentResidues=[]
    for pp0 in pp_list0:
      seq0= self._fromPolypeptideToSequence(pp0)
      chainId0= self._getChainId(pp0)
      # print(seq0)
      for pp1 in pp_list1:
        seq1= self._fromPolypeptideToSequence(pp1)
        # print(seq1)
        chainId1 = self._getChainId(pp1)
        # print(chainId0, chainId1)
        if excludeSameChainId and chainId0==chainId1:
          continue
        alignments = self._alignSeqs(seq0, seq1)
        if len(alignments)==0:
          continue
        ali0 = alignments[0][0]
        ali1 = alignments[0][1]
        identity = self._getIdentityPercent(ali0, ali1)
        # print("%.2f%%"%(100*identity))
        if identity <= IDENTITY_THR:
          continue
        index0 = -1
        index1 = -1

        for letter1, letter0 in zip(ali1, ali0):
          if letter1 != "-":
            index1 += 1
          if letter0 != "-":
            index0 += 1
          if letter1 != "-" and letter0 != "-":
            res0 = pp0[index0]
            res1 = pp1[index1]
            eqivalentResidues.append( (res0, res1) )

    return eqivalentResidues

class HomoOligomerFinderLR(HomoOligomerFinderBase):
  '''
  Used to add Interaction (resL, resR) if (resR, resL) is present and Ligand and Receptor are the same protein
  '''

  def __init__(self, pp_listL, pp_listR, interactingPairsResidues):
    '''
      :param pp_listL: [ Bio.PDB.Polypeptide]. list of Bio.PDB.Polypeptide's from Ligand
      :param pp_listR: [ Bio.PDB.Polypeptide]. list of Bio.PDB.Polypeptide's from Receptor

      :param interactingPairsResidues:  [(Bio.PDB.ResidueL, Bio.PDB.ResidueR)]: ligandResidueand receptorResidue [ Bio.PDB.Residue ]
    '''

    self.pp_listL= pp_listL
    self.pp_listR= pp_listR
    self.interactingPairsResidues= interactingPairsResidues

  def update_interactions(self):
    equivResiduesLR= self.findEquivalentResidues(self.pp_listL, self.pp_listR)
    equivResiduesLR_dict= dict(equivResiduesLR)
    equivResiduesRL_dict= {v: k for k, v in equivResiduesLR_dict.iteritems()}
    newInteractingSet=set()
    for resL, resR in self.interactingPairsResidues:
      try:
        reversed_resL= equivResiduesRL_dict[resR]
        reversed_resR= equivResiduesLR_dict[resL]
        newInteractingRes = (reversed_resL, reversed_resR)
        newInteractingSet.add( newInteractingRes )
      except KeyError:
        continue

    self.interactingPairsResidues= self.interactingPairsResidues.union(newInteractingSet)
    chainsInContactL, chainsInContactR = self._fromPairsInContactToChains(self.interactingPairsResidues)

    return self.interactingPairsResidues, (chainsInContactL, chainsInContactR )

class HomoOligomerFinderWithinPartner(HomoOligomerFinderBase):
  '''
  Used to compute resIds mapping from homooligomers proteins
  '''

  def __init__(self, pp_list, interactingPairsResidues, chainType):
    '''
      :param pp_list: [ Bio.PDB.Polypeptide]. list of Bio.PDB.Polypeptide's from receptor only or ligand only
      :param interactingPairsResidues:  [(Bio.PDB.ResidueL, Bio.PDB.ResidueR)]: ligandResidue and receptorResidue [ Bio.PDB.Residue ]
      :param chainType: str. "l" or "r"
    '''
    self.pp_list= pp_list
    self.interactingPairsResidues= interactingPairsResidues
    self.sequencesList= [ str(onePP.get_sequence()).replace("-","") for onePP in pp_list]
    self.chainIdsList=[onePP[0].get_parent().get_id() for onePP in pp_list]
    self.chainType=  chainType

    self.interactingResidue2ResidueMap = {}
    if chainType=="l":
        self.interactingResidue2ResidueMap= self._formatInteractingPairsAsDict()[0]
    elif chainType=="r":
      self.interactingResidue2ResidueMap = self._formatInteractingPairsAsDict()[1]
    else:
      raise ValueError("Error, chainType must be 'l' or 'r' or None")

  def update_interactions(self):
    '''
      fill in self.equivalentInteractingResidues  {res_id_interacting -->res_id_equivalent}
      To do so first computes _all_against_all_ali() and them matches chains that
      have the best structural aligment score. After that, residues of matches chains 
      are mapped thanks to a sequence alignment 
    '''  

    equivResidues= self.findEquivalentResidues(self.pp_list, self.pp_list, excludeSameChainId=True )

    if self.chainType == "l":
      def createPair(res0, otherRes):
        return (res0, otherRes)
    else:
      def createPair(res0, otherRes):
        return (otherRes, res0)

    newInteractingSet= set([])
    for res0, res1 in equivResidues:
      if res1 in self.interactingResidue2ResidueMap:
        # print(res0.get_full_id()[2:], res0.resname, res1.get_full_id()[2:], res1.resname)
        # print(sorted([ x.get_full_id()[2:] for x in self.interactingResidue2ResidueMap[res1]], key=lambda x: [1]))
        for partnerRes in self.interactingResidue2ResidueMap[res1]: #Check if this is ok
          newInteractingSet.add( createPair(res0, partnerRes) )

    self.interactingPairsResidues= self.interactingPairsResidues.union(newInteractingSet)
    chainsInContactL, chainsInContactR = self._fromPairsInContactToChains(self.interactingPairsResidues)

    return self.interactingPairsResidues, (chainsInContactL, chainsInContactR )


