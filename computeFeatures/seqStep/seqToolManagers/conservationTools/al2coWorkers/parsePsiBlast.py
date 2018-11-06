'''
Tested using PSIBLAST 2.2.31+ 

'''
import sys, os
from utils import openForReadingFnameOrGz

def parsePsiBlast( inputSeq, psiBlastOut, seqIdlowCut=30, seqIdUpCut=95, evalueLowCut=1e-4, lengthlocut=11, qfraclocut=66):
  '''
  @param inputSeq: str. An str that represents the input sequence
  @param psiBlastOut: str. Fname to psiblast alignments output
  @param seqIdlowCut: float. Lower cut-off for sequence identity with query (0 to 100%)
  @param seqIdUpCut: float. Upper cut-off for sequence identity with query (0 to 100%)
  @param evalueLowCut: float. Lower cut-off E-value
  @param lengthlocut: float. Lower cut-off for filtering alignments by length (in amino acids)
  @param qfraclocut: float. Lower cut-off for filtering alignments by percentage of query length (0 to 100%)
  
  @return List of psiblast hits [ hit ]
  
    hit: Dictionary that contains information about one pairwise alignment between the query and the target
      e.g.  { 'targetId': 'UniRef90_K3ZEZ8', 'targetStartNum': 8, 'targetEndNum': 73, 'queryStartNum': 1, 'queryEndNum': 70, 
              'targetSeq': 'TEETSGWTSWPEVVGMSVEEAKKVILKDKSDADIVVLPVGSPVALDLRLDRVRIFVD----TVAQTPHVG', 
              'querySeq':  'TEFGSELKSFPEVVGKTVDQAREYFTLHYPQYDVYFLPEGSPVTLDLRYNRVRVFYNPGTNVVNHVPHVG', 
              'evalue': 4e-19, 'target_full_id': 'UniRef90_K3ZEZ8_8_73', 'psiBlastRound': 3}
    
  '''
  inputSeqLen= len(inputSeq)
  with openForReadingFnameOrGz(psiBlastOut) as f:
    psiBlastRound=0
    hitsById={}
    allHits=[]

    querySeq=""
    targetSeq=""
    targetStartNum=99999
    targetEndNum= -99999
    queryStartNum=99999
    queryEndNum= -99999
    targetId=None
    identity, evalue, positives= [None]*3
    for line in f:
      if line.startswith("Results from round "):
        psiBlastRound+=1
        
      if psiBlastRound>0:
        if line.startswith(">"):
          qlen= len(querySeq.replace("-",""))
          if (identity >= seqIdlowCut and  identity <= seqIdUpCut and evalue <= evalueLowCut and 
                  qlen >= lengthlocut and   100.*qlen/inputSeqLen >= qfraclocut): #valid pairwise alig params
            target_full_id= targetId+"_%s_%s"%(targetStartNum, targetEndNum)
            targetInfo= {"querySeq":querySeq, "targetSeq":targetSeq, "targetStartNum":targetStartNum,
                         "targetEndNum":targetEndNum, "targetId":targetId, "psiBlastRound":psiBlastRound,
                         "queryStartNum":queryStartNum, "queryEndNum":queryEndNum, "evalue":evalue, 
                         "target_full_id":target_full_id}
            addValidTargetsWithId(targetInfo, hitsById, allHits)
          
          #Reset for new target
          querySeq=""
          targetSeq=""
          targetStartNum=99999
          targetEndNum= -99999
          queryStartNum=99999
          queryEndNum= -99999
          targetId= line.split()[1]

        if line.startswith(" Score = "):
          lineArray= line.split()
          score=  float(lineArray[2])
          evalue= float(lineArray[7][:-1])
          
        if line.startswith(" Identities ="):
          lineArray= line.split()
          identity=  int(lineArray[3][1:-3])
          positives= int(lineArray[7][1:-3])

        if line.startswith("Query "):
          lineArray= line.split()
          if len(lineArray)<3:
            querySeq+=lineArray[1]
          else:
            querySeq+=lineArray[2]          
            queryStartNum= min(int(lineArray[1]), queryStartNum)
            queryEndNum= max(int(lineArray[-1]), queryEndNum)
          
        if line.startswith("Sbjct "):
          lineArray= line.split()
          if len(lineArray)<3:
            querySeq+=lineArray[1]
          else:
            targetSeq+=lineArray[2]        
            targetStartNum= min(int(lineArray[1]), targetStartNum)
            targetEndNum= max(int(lineArray[-1]), targetEndNum)

  if (identity >= seqIdlowCut and  identity <= seqIdUpCut and evalue <= evalueLowCut and 
          qlen >= lengthlocut and   100.*qlen/float(inputSeqLen) >= qfraclocut): #valid pairwise alig params
    target_full_id= targetId+"_%s_%s"%(targetStartNum, targetEndNum)
    targetInfo= {"querySeq":querySeq, "targetSeq":targetSeq, "targetStartNum":targetStartNum,
                 "targetEndNum":targetEndNum, "targetId":targetId, "psiBlastRound":psiBlastRound,
                 "queryStartNum":queryStartNum, "queryEndNum":queryEndNum, "evalue":evalue, 
                 "target_full_id":target_full_id}
    addValidTargetsWithId(targetInfo, hitsById, allHits )
#  print( len(allHits)); raw_input("enter")
  return allHits
    


def overlap(previousEntry, currentEntry):
  '''
    Computes if two 
    @param: entry. Represents a pairwise alignment
  '''
  oldStart, oldEnd= previousEntry["targetStartNum"], previousEntry["targetEndNum"]
  newStart, newEnd= currentEntry["targetStartNum"], currentEntry["targetEndNum"]
  result =False
  maxAllowedOverlap = 0.5
  if not ( newEnd < oldStart) or ( newStart > oldEnd):
    startOverlap = oldStart
    if(newStart > oldStart):
      startOverlap = newStart
    endOverlap = oldEnd
    if(newEnd < oldEnd):
      endOverlap = newEnd
    lengthOld = oldEnd - oldStart
    lengthOverlap = endOverlap - startOverlap

    if lengthOverlap/float(lengthOld) > maxAllowedOverlap:
      result=True
  return result

def addValidTargetsWithId(targetInfo, hitsById, allHits):
  validTargets= []
  
  targetId= targetInfo["targetId"]
  targetStartNum, targetEndNum= targetInfo["targetStartNum"], targetInfo["targetEndNum"]
  if targetId in hitsById:
    otherTargetsSameId= hitsById[targetId]
    for otherTarget in otherTargetsSameId:
      if overlap(otherTarget, targetInfo):
        #redundant entry
        allHits.remove(otherTarget)
      else:
        validTargets+= [otherTarget]

  validTargets+= [ targetInfo ]

  hitsById[targetId]= validTargets
  allHits.append( targetInfo)
  
if __name__=="__main__":
  querySeq="TEFGSELKSFPEVVGKTVDQAREYFTLHYPQYDVYFLPEGSPVTLDLRYNRVRVFYNPGTNVVNHVPHVG"
  psiblastOutName= "/home/rsanchez/Tesis/rriPredMethod/data/develData/computedFeatures/seqStep/conservation/psiblast/1ACB_l_+_u_.psiblast"

#  querySeq="CGVPAIQPVLSGLSRIVNGEEAVPGSWPWQVSLQDKTGFHFCGGSLINENWVVTAAHCGVTTSDVVVAGEFDQGSSSEKIQKLKIAKVFKNSKYNSLTINNDITLLKLSTAASFSQTVSAVCLPSASDDFAAGTTCVTTGWGLTRYTNANTPDRLQQASLPLLSNTNCKKYWGTKIKDAMICAGASGVSSCMGDSGGPLVCKKNGAWTLVGIVSWGSSTCSTSTPGVYARVTALVNWVQQTLAAN"
#  psiblastOutName= "/home/rsanchez/Tesis/rriPredMethod/data/develData/computedFeatures/seqStep/conservation/psiblast/1ACB_r_B_u_.psiblast"

  parsePsiBlast( querySeq, psiblastOutName, seqIdlowCut=30, seqIdUpCut=95, evalueLowCut=1e-4, lengthlocut=30, qfraclocut=66)
      

