from scipy.spatial import distance
from scipy.spatial import Delaunay
from ..featuresComputer import FeatureComputerException
     
def getVoroNeigs(listOfResidues, maxDist):
  '''
  given a list of residues computes voronoi neighbours closer than maxDist angstroms
  :param listOfResidues: [ Bio.PDB.Residue ]
  :param maxDist: float. In angstrongs
  :return:  { Bio.PDB.Residue: set( [ Bio.PDB.Residue ] } Dict from residue to neigs residues
  
  '''
  coords= []
  residues= []
  for res in listOfResidues:
    try:
      coords.append( res["CA"].get_coord())
      residues.append(res)
    except KeyError:
      try:
        coords.append( res["CB"].get_coord())
        residues.append(res)
      except KeyError:
        coords.append( res.get_list()[0].get_coord())
        residues.append(res)
        
  delaTriang= Delaunay(coords, furthest_site=False, incremental=False, qhull_options=None)
  indices, indptr= delaTriang.vertex_neighbor_vertices
  neigsDict={}
  for k in range(len(coords)):
    residue= residues[k]
    neigsIdx= indptr[indices[k]:indices[k+1]]
    neigsDict[residue] = set( [residues[neigIdx] for neigIdx in neigsIdx 
                              if distance.euclidean(coords[k], coords[neigIdx])< maxDist])
  return neigsDict


