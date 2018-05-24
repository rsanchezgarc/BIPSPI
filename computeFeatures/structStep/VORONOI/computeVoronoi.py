import os,sys
from computeFeatures.structStep.myPDBParser import myPDBParser as PDBParser
from Bio.PDB.Polypeptide import is_aa
from subprocess import Popen, PIPE
from scipy.spatial import distance

from ..StructFeatComputer import StructFeatComputer
from utils import myMakeDir, tryToRemove #utils is at the root of the package


class VORONOIComputer(StructFeatComputer):
  '''
  Extends StructFeatComputer class. Computes VORONOI neighbours for receptor and ligand independently
  '''
  def __init__(self, rFname, lFname, computedFeatsRootDir= None, maxDist= 30, statusManager=None):
    '''
      @param rFname: str. path to receptor pdb file
      @param lFname: str. path to ligand pdb file      
      @param computedFeatsRootDir: str. path where features will be stored
      @param statusManager: class that implements .setStatus(msg) to communicate
    '''
    StructFeatComputer.__init__(self, rFname, lFname, computedFeatsRootDir, statusManager=statusManager)
    self.maxDist= maxDist
    self.outPath= myMakeDir(self.computedFeatsRootDir, "VORONOI/VORONOI_"+str(self.maxDist))
    
    self.qhullExeDir = os.path.join(self.qhullBinDir,"qdelaunay")
   
  def computeOneFile(self, fileName):
    '''
      Computes VORONOI neighbours for a given pdb file
      @param fileName: str. fname to pdb file
    '''
    
    voro_file= os.path.split(fileName)[-1]
    voro_file= voro_file.split(".")[0]+ ".voro"
    outName= os.path.join(self.outPath, voro_file )
    if os.path.isfile(outName):
      print ('Already computed VORONOI')
      return 0
    
    struct= PDBParser(QUIET= True).get_structure("oneStruct", fileName)
    
    ids= []
    coords= []
    for res in struct[0].get_residues():
      if not is_aa(res, standard= True): continue
      structName, modelId, chainId, resIdTuple = res.get_full_id()
      if resIdTuple[2]!= " ": continue
      if chainId==" ": chainId="*"
      resIdRepr= str(resIdTuple[1]) +"_"+ chainId
      try:
        coords.append( res["CA"].get_coord())
        ids.append(resIdRepr)
      except KeyError:
        try:
          coords.append( res["CB"].get_coord())
          ids.append(resIdRepr)
        except KeyError:
          coords.append( res.get_list()[0].get_coord())
          ids.append(resIdRepr)

    inStream= "\n".join([ " ".join([str(elem) for elem in coord_line]) for coord_line in  coords ])
    inStream= "3  #sample 3-d input\n"+str(len(coords))+"\n"+inStream
 
    proc= Popen( "%s Fv Qt | cut -d\" \" -f2-5"%( self.qhullExeDir), shell=True, stdin=PIPE, stdout=PIPE,
                 stderr=PIPE)
    out= proc.communicate( inStream)

    if not out[0][0].isdigit() or len(out[1])>0:
      raise FeatureComputerException("Error in qhull execution for %s: %s"%(fileName, out))
    out= out[0]
    lines= out.split("\n")
    try:
      with open(outName, "w") as f:
        for line in lines[1:]:

          splitLine= [int(elem) for elem in line.split()]
          if len(splitLine)<1: 
            continue
            print ("badLine")
          alreadyPrinted= set([])
          for resAndChainInd1 in splitLine:
            for resAndChainInd2 in splitLine:
              if resAndChainInd1 != resAndChainInd2:
                dist= distance.euclidean(coords[resAndChainInd1], coords[resAndChainInd2])
                if dist< self.maxDist:
                  if (ids[resAndChainInd1], ids[resAndChainInd2]) not in alreadyPrinted:
                    f.write("%s\t%s\n"%(ids[resAndChainInd1], ids[resAndChainInd2]))
                    alreadyPrinted.add((ids[resAndChainInd1], ids[resAndChainInd2]))

                  if (ids[resAndChainInd2], ids[resAndChainInd1]) not in alreadyPrinted:                
                    f.write("%s\t%s\n"%(ids[resAndChainInd2], ids[resAndChainInd1]))
                    alreadyPrinted.add((ids[resAndChainInd2], ids[resAndChainInd1]))
    except (KeyboardInterrupt, Exception):
      print("Exception happend computing %s"%outName)
      tryToRemove(outName)    
      raise                        
    return 0

                
def testModule():
#  rFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_r_u.pdb"
#  lFname="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3/1A2K_l_u.pdb"
  
  rFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1A22.pdb"
  lFname="/home/rsanchez/Tesis/managePDB/data/dimersAndTrimersRawPDBs/1ABR.pdb"
  
  computedFeatsRootDir="/home/rsanchez/Tesis/rriPredMethod/data/ppdockingBenchData/tmpComputedFeatures"
  voroComp= VORONOIComputer(rFname, lFname, computedFeatsRootDir)
  voroComp.computeFun()
  raise FeatureComputerException("Debug Mode. Comment out testModule() line 128")                
                    
if __name__=="__main__":
#  testModule()
  
  maxDist= 10
  
  pdbsIndir= "~/Tesis/rriPredMethod/data/ppdockingBenchData/BenchmarkData/rawBenchmark3"
  computedFeatsRootDir= "~/Tesis/rriPredMethod/data/ppdockingBenchData/computedFeatures"
  
#  pdbsIndir= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/b5structures/"
#  computedFeatsRootDir= "/home/rsanchez/Tesis/rriPredMethod/data/bench5Data/computedFeatures"
  

  StructFeatComputer.computeAllComplexes(VORONOIComputer,pdbsIndir= pdbsIndir ,computedFeatsRootDir= computedFeatsRootDir,
            classArgs= (maxDist,) )


