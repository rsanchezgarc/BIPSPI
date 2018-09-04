from __future__ import absolute_import
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio.SVDSuperimposer import SVDSuperimposer
import numpy as np


class BoundUnboundMapper(object):
  '''
  Used to compute resIds mapping from bound state to unbound state
  '''
  scoreMat= matlist.blosum62
  def __init__(self, pp_list_unbound, pp_list_bound):
    '''
      @param pp_list_unbound: [ Bio.PDB.Polypeptide]. list of Bio.PDB.Polypeptide's from unbound structure
      @param pp_list_bound: [ Bio.PDB.Polypeptide]. list of Bio.PDB.Polypeptide's from bound structure
    '''
    self.pp_list_unbound= pp_list_unbound
    self.pp_list_bound= pp_list_bound
    self.sequencesUnbound= [ str(elem.get_sequence()).replace("-","") for elem in pp_list_unbound]
    self.sequencesBound= [ str(elem.get_sequence()).replace("-","") for elem in pp_list_bound]

    self.boundToUnboundDict= {}
    self.boundToUnboundDictFromChainResId= None
    self.unboundToBoundToDictFromChainResId= None

  def build_correspondence(self):
    '''
      fill in self.boundToUnboundDict  {res_id_bound -->res_id_unbound}
      To do so first computes _all_against_all_ali() and them matches chains that
      have the best structural aligment score. After that, residues of matches chains 
      are mapped thanks to a sequence alignment 
    '''  
    aligU2BDictsTable, aligU2BScores = self._all_against_all_alig()
    n2Take= max(aligU2BScores.shape)
    maxVal= np.max(aligU2BScores)+ 1e-4
    bestScore= np.min(aligU2BScores)
    bestScoreThr= max(10, bestScore)
#    print(n2Take, bestScoreThr)
#    print(aligU2BScores)    
    while n2Take>0 and bestScore<= bestScoreThr:
      argbest= np.where(aligU2BScores==bestScore)
      unboundInd= argbest[0][0]
      boundInd=   argbest[1][0]
#      print(aligU2BScores)
#      print("Map in this step:", unboundInd, boundInd)
      aligU2BScores[unboundInd, boundInd]+= maxVal
      n2Take-=1
      bestScore= np.min(aligU2BScores)
      self.boundToUnboundDict.update(aligU2BDictsTable[unboundInd][boundInd])

#    print(aligU2BScores)
####    print(self.boundToUnboundDict)
#    raw_input("enter")
    return None
    
  def build2SeqsDictMap(self,  nSeqUnbound, seq_u, nSeqBound, seq_b):
    '''
      reads a seq aligment between 2 sequences and creates 2 dictionaries
      that maps bound to unbound residues and CA atoms boundToUnboundResDict, atomBoundToUnboundMap.
      @param nSeqUnbound: int. The index of the bound sequence that will be aligned
      @param seq_u: str. The alignment result for unbound sequence number nSeqUnbound
      @param nSeqBound: int. The index of the bound sequence that will be aligned
      @param seq_b: str. The alignment result for bound sequence number nSeqBound
      @return boundToUnboundResDict. {Bio.PDB.Residue_bound --> Bio.PDB.Residue_unbound}
      @return atomBoundToUnboundMap. {Bio.PDB.Atom_bound --> Bio.PDB.Atom_unbound} CA atoms    
    '''    
    pp_u= self.pp_list_unbound[nSeqUnbound]
    pp_b= self.pp_list_bound[nSeqBound]
    boundToUnboundResDict= {}
    atomBoundToUnboundMap=[]
    b_index= -1
    u_index=-1    
    for u_letter, b_letter in zip(seq_u, seq_b):
  #        print(u_letter, b_letter, u_index, b_index)
  #        raw_input()
      if u_letter !="-":
        u_index+=1
      if b_letter !="-":
        b_index+=1
      if u_letter !="-" and b_letter !="-":
        boundToUnboundResDict[ pp_b[b_index] ]= pp_u[u_index]
#        print pp_b[b_index], pp_u[u_index]
        atom_b= pp_b[b_index]["CA"] if "CA" in pp_b[b_index] else None
        atom_u= pp_u[u_index]["CA"] if "CA" in pp_u[u_index] else None
        if not atom_b is None and not atom_u is None:    
          atomBoundToUnboundMap+= [ (atom_b, atom_u) ]
    return boundToUnboundResDict, atomBoundToUnboundMap
    
  def _all_against_all_alig(self):
    '''
    Computes rmsd for all bound chains against all unbound chains and builds a table of
    rmsd (aligU2BScores)
    Do seq aligment of all bound seqs agains all unbound seqs and builds a table of
    seq aligments aligU2BTable
    @return aligU2BTable: nxm list of floats, being each element a dict that maps each of the residues of
                           chain_u_i with their equivalent residues chain_b_i
    @return aligU2BScores: nxm list of floats, being each element rmsd of structural aligment of
                           chain_u_i chain_b_i
    '''     
    aligU2BTable= []
    aligU2BScores=[]
    for nSeqUnbound, seqUnbound in enumerate(self.sequencesUnbound):
      aligU2BTable.append([])
      aligU2BScores.append([])
      for nSeqBound,seqBound in enumerate(self.sequencesBound):
        try:
          seqUnboundAli, seqBoundAli, SeqAligScore = self._alig_seq(seqUnbound, seqBound)
        except IndexError:
          continue
        aligScore, boundToUnboundDict= self.getRMSD(nSeqUnbound, seqUnboundAli, nSeqBound, seqBoundAli)
        aligU2BTable[-1].append( boundToUnboundDict )
        aligU2BScores[-1].append(aligScore)

    return aligU2BTable, np.array(aligU2BScores)

  def getRMSD(self, nSeqUnbound, seqUnboundAli, nSeqBound, seqBoundAli):
    '''
    Computes rmsd for nSeqUnbound chain unbound and nSeqBound bound chain
    @param nSeqUnbound: int. The index of the bound sequence that will be aligned
    @param seqUnboundAli: str. The alignment result for unbound sequence number nSeqUnbound
    @param nSeqBound: int. The index of the bound sequence that will be aligned
    @param seqBoundAli: str. The alignment result for bound sequence number nSeqBound
    @return rmsd. float. Root mean square deviation of CA of both imput chains
    @return boundToUnboundResDict. {Bio.PDB.Residue_bound --> Bio.PDB.Residue_unbound}
    '''     
    boundToUnboundResDict, atomBoundToUnboundMap= self.build2SeqsDictMap(nSeqUnbound, seqUnboundAli, nSeqBound, seqBoundAli)
    atoms_x, atoms_y= zip(* atomBoundToUnboundMap)
    coords_x= np.array([elem.get_coord() for elem in atoms_x])
    coords_y= np.array([elem.get_coord() for elem in atoms_y])
    sup= SVDSuperimposer()
    rmsd= sup._rms(coords_x, coords_y)
#    print(boundToUnboundResDict)
    return rmsd, boundToUnboundResDict
    
  def _alig_seq(self, seq1, seq2):
    '''
    aligns seq1 against seq2.
    @param seq1: str. Sequence 1
    @param seq1: str. Sequence 2
    @return seq1Aligment. str.
    @return seq2Aligment. str
    @return seqScore. float    
    '''      
    alignments = pairwise2.align.localds(seq1, seq2,BoundUnboundMapper.scoreMat, -11, -0.5)
#    print(alignments[0][0])
#    print(alignments[0][1])
#    print("+++++++++++++++++++++++++++++++++++++++++++++")
#    raw_input()
    return alignments[0][0], alignments[0][1], alignments[0][2]

  def mapBoundToUnbound(self, boundRes):
    '''
    Gets the equivalent unbound residue for a given bound residue
    
    @param boundRes: Bio.PDB.Residue of bound structure
    @return None if there is no equivalent residue or Bio.PDB.Residue if 
              there is an equivalent unbound residue
    '''
    try:
      return self.boundToUnboundDict[boundRes]
    except KeyError:
      return None

  def convertDictResis2Ids(self, residuesDict):
    boundToUnboundDictFromChainResId= {  res.get_full_id(): residuesDict[res].get_full_id() 
                                                  for res in residuesDict}
    boundToUnboundDictFromChainResId= \
        {(res_id[2], str(res_id[3][1])+ res_id[3][2].strip() ): (boundToUnboundDictFromChainResId[res_id][2],
            str(boundToUnboundDictFromChainResId[res_id][3][1])+ boundToUnboundDictFromChainResId[res_id][3][2].strip()) 
                                    for res_id in boundToUnboundDictFromChainResId}  
    return boundToUnboundDictFromChainResId
    
  def mapBoundToUnboundUsingId(self, chainId, resId):
    '''
    Gets the equivalent unbound residue for a given bound residue
    
    @param chainId: str. Chain id for bound residue
    @param resId: str. res id for bound residue
    @return None if there is no equivalent residue or Bio.PDB.Residue if 
              there is an equivalent unbound residue    
    '''
    
    if self.boundToUnboundDictFromChainResId is None:
      self.boundToUnboundDictFromChainResId= self.convertDictResis2Ids(self.boundToUnboundDict)
    try:
      return self.boundToUnboundDictFromChainResId[(chainId, resId)]
    except KeyError:
      return None
      
  def mapUnboundToBoundUsingId(self, chainId, resId):
    '''
    Gets the equivalent bound residue for a given unbound residue
    
    @param chainId: str. Chain id for unbound residue
    @param resId: str. res id for unbound residue
    @return None if there is no equivalent residue or Bio.PDB.Residue if 
              there is an equivalent unbound residue    
    '''
    if self.boundToUnboundDictFromChainResId is None:
      self.boundToUnboundDictFromChainResId= self.convertDictResis2Ids(self.boundToUnboundDict)
      
    if self.unboundToBoundToDictFromChainResId is None:
      self.unboundToBoundToDictFromChainResId= { self.boundToUnboundDictFromChainResId[res_id]: res_id
                                                  for res_id in self.boundToUnboundDictFromChainResId}
    try:
      return self.unboundToBoundToDictFromChainResId[(chainId, resId)]
    except KeyError:
      return None

