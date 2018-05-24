'''
from https://stackoverflow.com/questions/11685716/how-to-extract-chains-from-a-pdb-file
'''
import os
from Bio import PDB
try:
  from computeFeatures.structStep.myPDBParser import myPDBParser as PDBParser
except ImportError:
  from Bio.PDB.PDBParser import PDBParser

class ChainSplitter:
    def __init__(self, out_dir=None):
        """ Create parsing and writing objects, specify output directory. """
        self.parser = PDBParser(QUIET=True)
        self.writer = PDB.PDBIO()
        if out_dir is None:
            out_dir = os.path.join(os.getcwd(), "chain_PDBs")
        self.out_dir = out_dir

    def make_pdb(self, pdb_path, chain_letters, overwrite=False, struct=None, rejectInsteadAccept=False):
        """ Create a new PDB file containing only the specified chains.

        Returns the path to the created file.

        :param pdb_path: full path to the crystal structure
        :param chain_letters: iterable of chain characters (case insensitive)
        :param overwrite: write over the output file if it exists
        :param rejectInsteadAccept: use SelectChainsReversed instead SelectChains
        """
        chain_letters = [chain.upper() for chain in chain_letters]

        # Input/output files
        (pdb_dir, pdb_fn) = os.path.split(pdb_path)
        pdb_id = os.path.split(pdb_fn)[-1].split(".")[0]
        sep= "!" if rejectInsteadAccept else ":" 
        out_name = "%s%s%s.pdb" % (pdb_id, sep, "".join(chain_letters))
        out_path = os.path.join(self.out_dir, out_name)
        print ("OUT PATH: %s"%out_path)
        plural = "s" if (len(chain_letters) > 1) else ""  # for printing

        # Skip PDB generation if the file already exists
        if (not overwrite) and (os.path.isfile(out_path)):
            print("Chain%s %s of '%s' already extracted to '%s'." %
                    (plural, ", ".join(chain_letters), pdb_id, out_name))
            return out_path

        print("Extracting chain%s %s from %s..." % (plural,
                ", ".join(chain_letters), pdb_fn))

        # Get structure, write new file with only given chains
        if struct is None:
            struct = self.parser.get_structure(pdb_id, pdb_path)
        abortExecution=True
        for chain in struct.get_chains():
          if chain.get_id() in chain_letters:
            abortExecution=False
        if abortExecution: raise ValueError("Error: chains '%s' are not contained in pdb"%(", ".join(chain_letters)))
        self.writer.set_structure(struct)
        if rejectInsteadAccept:
            self.writer.save(out_path, select=SelectChainsReversed(chain_letters))
        else:
            self.writer.save(out_path, select=SelectChains(chain_letters))

        return out_path


class SelectChains(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (chain.get_id() in self.chain_letters)


class SelectChainsReversed(PDB.Select):
    """ Only accept the specified chains when saving. """
    def __init__(self, chain_letters):
        self.chain_letters = chain_letters

    def accept_chain(self, chain):
        return (not chain.get_id() in self.chain_letters)
        
if __name__ == "__main__":
    """ Parses PDB id's desired chains, and creates new PDB structures. """
    import sys
    if len(sys.argv) != 5:
        print "Usage: $ python %s 'pdb_in' 'chain' 'pdb_outDir' '0/1 for exclude/include'" % __file__
        sys.exit()
    else:
        pdb_path= sys.argv[1]
        chain= sys.argv[2]
        pdbOutDir= sys.argv[3]
        rejectInsteadAccept= True if sys.argv[4]=="1" else False
    splitter = ChainSplitter(pdbOutDir)
    pdbOutName= splitter.make_pdb(pdb_path, chain, rejectInsteadAccept=rejectInsteadAccept)


