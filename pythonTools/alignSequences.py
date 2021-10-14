#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 20 13:11:17 2020

@author: ruben
"""
import sys
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo

def getMatchingSeqsIndices( targetSeq, referenceSeq, maxMismatchFracAllowed=0.1):
  '''
    targetSeq= "aa0aa1aa2..."  #Sequences whose indices that match referenceSeq we want to obtain
    referenceSeq= "aa0aa1aa2..."

    return: matchingIndices  target  --> ref
  '''
  targetSeq= "".join( targetSeq).replace("-","X")
  referenceSeq= "".join( referenceSeq).replace("-","X")
  matchingIndices= {}
  alignments = pairwise2.align.localds(targetSeq, referenceSeq, MatrixInfo.blosum62, -11, -0.5)
  # print(pairwise2.format_alignment(*alignments[0])); input("enter")
  target_ali, ref_ali = alignments[0][0], alignments[0][1]
  target_ix=-1
  ref_ix=-1
  nMatches=0
  nMismatches=0
  for aa_target, aa_ref in zip(target_ali, ref_ali):
    if aa_target !="-":
      target_ix+=1
    if aa_ref !="-":
      ref_ix+=1
    if aa_target !="-" and aa_ref !="-":
      matchingIndices[target_ix]= ref_ix
      if aa_target==aa_ref:
        nMatches+=1
      else:
        nMismatches+=1

  seqIdentity= nMatches/float(nMatches+nMismatches)
  nXs= sum([1 for elem in targetSeq if elem=="X"])
  if nXs> len(targetSeq)*0.5: #skip seq if contains more than 50% of non aminoacids
    raise ValueError("To many X's in sequence")
  assert abs(len(matchingIndices) - len(referenceSeq)) < maxMismatchFracAllowed*len(referenceSeq
                                    ), "Error, excesive mismatch between \n%s\n+++++++++++++++++\n%s\n"%(targetSeq, referenceSeq)
  return matchingIndices, alignments[0][2], seqIdentity

if __name__=="__main__":
  targetSeq, referenceSeq = sys.argv[1:]
  matchingIdxs= getMatchingSeqsIndices(targetSeq, referenceSeq)
  print(matchingIdxs)