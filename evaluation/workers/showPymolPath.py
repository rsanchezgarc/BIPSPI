import pymol
import os
from Bio.PDB.Polypeptide import one_to_three, PPBuilder
from computeFeatures.structStep.myPDBParser import myPDBParser as PDBParser
from computeFeatures.common.boundUnboundMapper import BoundUnboundMapper


def getBoundResList(fname_bound, fname_unbound, listOfDictsChainToResId):
  parser= PDBParser(QUIET=True)
  structureUnbound= parser.get_structure(fname_unbound, fname_unbound)
  structureBound= parser.get_structure(fname_bound, fname_bound)  
  ppb= PPBuilder()
  pp_list_unbound= ppb.build_peptides(structureUnbound, aa_only= False)    
  pp_list_bound=   ppb.build_peptides(structureBound, aa_only= False)    
  mapper= BoundUnboundMapper( pp_list_unbound,pp_list_bound)
  mapper.build_correspondence()
  newDictsList=[]

  for dictOfChainsToRes in listOfDictsChainToResId:
    tempDict={}
    for chainId_u in dictOfChainsToRes:
      for resId_u in sorted(dictOfChainsToRes[chainId_u]):
        chainId_b_resId_b= mapper.mapUnboundToBoundUsingId(" " if chainId_u=="*" else chainId_u, resId_u)
#        print(chainId_u, resId_u, chainId_b_resId_b)
        if chainId_b_resId_b is None: continue
        chainId_b, resId_b= chainId_b_resId_b
        if not chainId_b in tempDict:
          tempDict[chainId_b]=[]
        tempDict[chainId_b].append( resId_b)
    newDictsList.append( tempDict)
  return newDictsList
    
def showPDB_Worker_Interface(dictOfChains, structIdPymol, isLigand=True, isBound=False, color="green", name="unknown"):

  colors=["red","blue", "magenta","darksalmon","lightblue", "orange", "wheat", "purpleblue","warmpink"]
  n_colors= len(colors)
  patchWord= "Lpatch"+name if isLigand else "Rpatch"+name
  if isBound:
    patchWord+="_b_"
  commandOneChain=""
  for i,chain in enumerate(dictOfChains):
#        print( "chain:", repr(chain))
    commandOneChain= ", model "+ structIdPymol +" and (chain "+chain+" and resi "
    residuesByChain=""
    for res in dictOfChains[chain]:
      residuesByChain+= res+"+"
    residuesByChain= residuesByChain[:-1] +" )"
    selectStringTmp = commandOneChain+ residuesByChain
    print("SelectStr:", patchWord+chain, selectStringTmp)
    pymol.cmd.select(patchWord+chain, selectStringTmp)
#    pymol.cmd.show("sticks",patchWord+chain)
##    pymol.cmd.color( colors[i], patchWord+chain )
    pymol.cmd.color( color, patchWord+chain )
    pymol.cmd.deselect()

def fromListToDict(interfaceList):
  dictRes={}
  for residue in interfaceList:
    chain, resIndex = residue
    try:
      dictRes[chain].append(str(resIndex))
    except KeyError:
      dictRes[chain]=[str(resIndex)]
  return dictRes
  
def showPDB_patches_all(pdbPath, complexName, res_pred_l, res_pred_r, res_true_l, res_true_r):

  rnameFile= os.path.join(pdbPath,complexName+"_r_u.pdb")
  rname= complexName+"_r_u"

  lnameFile= os.path.join(pdbPath,complexName+"_l_u.pdb")
  lname= complexName+"_l_u"

###   Load Structures
  pymol.finish_launching()

  pymol.cmd.load(lnameFile,lname)
  pymol.cmd.show_as("cartoon",lname)

  pymol.cmd.load(rnameFile,rname)
  pymol.cmd.show_as("cartoon",rname)

  pymol.cmd.select("ml"," model "+ lname)
  pymol.cmd.color("grey","ml")
  pymol.cmd.select("mr"," model "+ rname)
  pymol.cmd.color("palegreen","mr")

#  pymol.cmd.util.cbc("*")
  res_pred_l, res_true_l= set(res_pred_l), set(res_true_l)
  res_pred_r, res_true_r= set(res_pred_r), set(res_true_r)
  
  if len(res_true_l)>0:
    l_truePos= fromListToDict(res_pred_l.intersection(res_true_l))
    l_falsPos= fromListToDict(res_pred_l.difference(res_true_l))
    l_falsNeg= fromListToDict(res_true_l.difference(res_pred_l))
  else:
    l_truePos= fromListToDict(res_pred_l)

  if len(res_true_r)>0:      
    r_truePos= fromListToDict(res_pred_r.intersection(res_true_r))
    r_falsPos= fromListToDict(res_pred_r.difference(res_true_r))
    r_falsNeg= fromListToDict(res_true_r.difference(res_pred_r)) 
  else:
    r_truePos= fromListToDict(res_pred_r)

  showPDB_Worker_Interface(l_truePos, lname, isLigand=True, color= "blue", name="truePos")
  if len(res_true_l)>0:  
    showPDB_Worker_Interface(l_falsNeg, lname, isLigand=True, color= "cyan", name="falseNeg")
    showPDB_Worker_Interface(l_falsPos, lname, isLigand=True, color= "lime", name="falsePos")
  
  showPDB_Worker_Interface(r_truePos, rname, isLigand=False, color="red", name="truePos")
  if len(res_true_r)>0:  
    showPDB_Worker_Interface(r_falsNeg, rname, isLigand=False, color="violet", name="falseNeg")  
    showPDB_Worker_Interface(r_falsPos, rname, isLigand=False, color="yelloworange", name="falsePos")

  
  rnameFile_bound= os.path.join(pdbPath,complexName+"_r_b.pdb")
  rname_bound= complexName+"_r_b"

  lnameFile_bound= os.path.join(pdbPath,complexName+"_l_b.pdb")
  lname_bound= complexName+"_l_b"
  if os.path.isfile(rnameFile_bound) and os.path.isfile(lnameFile_bound):
    pymol.cmd.load(lnameFile_bound,lname_bound)
    pymol.cmd.show_as("cartoon",lname_bound)

    pymol.cmd.load(rnameFile_bound,rname_bound)
    pymol.cmd.show_as("cartoon",rname_bound)

    pymol.cmd.select("ml_bound"," model "+ lname_bound)
    pymol.cmd.color("grey80","ml_bound")
    pymol.cmd.select("mr_bound"," model "+ rname_bound)
    pymol.cmd.color("smudge","mr_bound")


    if len(res_true_l)>0:      
      l_truePos, l_falsPos, l_falsNeg= getBoundResList(lnameFile_bound, lnameFile, [l_truePos, l_falsPos, l_falsNeg])  
    else:
      l_truePos = getBoundResList(rname_bound, rname, [l_truePos])[0]
    

    if len(res_true_r)>0:      
      r_truePos, r_falsPos, r_falsNeg= getBoundResList(rnameFile_bound, rnameFile, [r_truePos, r_falsPos, r_falsNeg])  
    else:
      r_truePos = getBoundResList(rname_bound, rname, [r_truePos])[0]
      
    showPDB_Worker_Interface(l_truePos, lname_bound, isLigand=True, isBound=True, color= "blue", name="truePos")
    if len(res_true_l)>0:  
      showPDB_Worker_Interface(l_falsNeg, lname_bound, isLigand=True, isBound=True, color= "cyan", name="falseNeg")
      showPDB_Worker_Interface(l_falsPos, lname_bound, isLigand=True, isBound=True, color= "lime", name="falsePos")
    
    showPDB_Worker_Interface(r_truePos, rname_bound, isLigand=False, isBound=True, color="red", name="truePos")
    if len(res_true_r)>0:  
      showPDB_Worker_Interface(r_falsNeg, rname_bound, isLigand=False, isBound=True, color="violet", name="falseNeg")  
      showPDB_Worker_Interface(r_falsPos, rname_bound, isLigand=False, isBound=True, color="yelloworange", name="falsePos")



def showPDB_interactions(pdbPath, complexName, l_r_pairsres , showPairs=True):
  rnameFile= os.path.join(pdbPath,complexName+"_r_u.pdb")
  rname= complexName+"_r_u"

  lnameFile= os.path.join(pdbPath,complexName+"_l_u.pdb")
  lname= complexName+"_l_u"

###   Load Structures
  pymol.finish_launching()

  pymol.cmd.load(lnameFile,lname)
  pymol.cmd.show_as("cartoon",lname)

  pymol.cmd.load(rnameFile,rname)
  pymol.cmd.show_as("cartoon",rname)

  pymol.cmd.select("ml"," model "+ lname)
  pymol.cmd.color("grey","ml")
  pymol.cmd.select("mr"," model "+ rname)
  pymol.cmd.color("palegreen","mr")
  if showPairs:
    for i, ((chainL, resL, resNameL), (chainR, resR, resNameR)) in enumerate(l_r_pairsres):
      resNameL=  one_to_three(resNameL)
      resNameR=  one_to_three(resNameR)
      pymol.cmd.distance( "dist%d"%i, " model %s and chain %s and resi %s and resn %s and name CA"%(lname, chainL, resL, resNameL),
                                " model %s and chain %s and resi %s and resn %s and name CA"%(rname, chainR, resR, resNameR),
                    )
      print( (chainL, resL, resNameL), (chainR, resR, resNameR))
  else:
    res_pred_l, res_pred_r = zip(* l_r_pairsres)
    res_pred_l= [ elem[:2] for elem in res_pred_l]
    res_pred_r= [ elem[:2] for elem in res_pred_r]
    showPDB_patches_all(pdbPath, complexName, res_pred_l, res_pred_r, res_true_l={}, res_true_r={})
    
