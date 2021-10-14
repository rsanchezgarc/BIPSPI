import argparse

from Bio.PDB import PDBIO
from Bio.PDB.vectors import Vector, rotaxis
import numpy as np

from .myPDBParser import myPDBParser

def rotateTranslatePdb(fname, rotX, rotY, rotZ, transX, transY, transZ, fnameOut=None):
  print(fname, rotX, rotY, rotZ, transX, transY, transZ, fnameOut)
  struct= myPDBParser(QUIET=True).get_structure(fname)

  rotationX = rotaxis( -rotX, Vector(1, 0, 0))
  rotationY = rotaxis( rotY, Vector(0, 1, 0))
  rotationZ = rotaxis( -rotZ, Vector(0, 0, 1))

  translation = np.array((transX, transY, transZ), 'f')

  rotation = rotationX.dot(rotationY).dot(rotationZ)
  struct.transform(rotation, translation)

  if fnameOut is not None:
    fnameOut= fnameOut
    pdbWriter= PDBIO()
    pdbWriter.set_structure(struct)
    pdbWriter.save(fnameOut)

  return struct

if __name__=="__main__":

  description="Translates and rotates a pdb file using the patchdock axes convention"
  parser = argparse.ArgumentParser(description)
  parser.add_argument('-i', '--pdbIn', type=str, nargs=1, required=True,
                      help='pdb to rotate\n')
  parser.add_argument('-r', '--rotXYZ', type=str, required=False, default="0,0,0",
                      help='rotation angles for each axis. E.g.: "-r 2,0,0". '
                           'Rotation of 2 rads over axis X')

  parser.add_argument('-t', '--transXYZ', type=str, required=False, default="0,0,0",
                      help='translation for each axis E.g.: "-t 20,0,0". '
                           'Translation of 20 A along axis X')

  parser.add_argument('-d', '--anglesInDegrees', action='store_true', required=False,
                      default=False, help='flag to consider angles as degrees')

  parser.add_argument('-o', '--pdbOut', type=str, nargs=1, required=False,default=None,
                      help='file to store rotated pdb\n')
  args= parser.parse_args()

  rots= [ float(elem) for elem in args.rotXYZ.split(",") ]
  transX, transY, transZ = [float(elem) for elem in args.transXYZ.split(",")]
  print(args)
  if args.anglesInDegrees:
    rots=[ np.pi*elem/180. for elem in rots]
  rotX, rotY, rotZ= rots
  fnameOut = None if args.pdbOut is None else args.pdbOut[0]

  rotateTranslatePdb(args.pdbIn[0], rotX, rotY, rotZ, transX, transY, transZ, fnameOut)