# Modified from BioPython
# Copyright (C) 2002, Thomas Hamelryck (thamelry@binf.ku.dk)
# This code is part of the Biopython distribution and governed by its
# license.  Please see the LICENSE file that should have been included
# as part of this package.

# This is a fix for HSExposureCB that crashes where central residue is Gly and there is no N or CA.


"""Half-sphere exposure and coordination number calculation."""

from __future__ import print_function

import warnings
from math import pi

from Bio.PDB.AbstractPropertyMap import AbstractPropertyMap
from Bio.PDB.Polypeptide import CaPPBuilder, is_aa

from Bio.PDB.HSExposure import _AbstractHSExposure, HSExposureCA

class HSExposureCB(_AbstractHSExposure):
    """Class to calculate HSE based on the real CA-CB vectors."""

    def __init__(self, model, radius=12, offset=0):
        """Initialize class.

        @param model: the model that contains the residues
        @type model: L{Model}

        @param radius: radius of the sphere (centred at the CA atom)
        @type radius: float

        @param offset: number of flanking residues that are ignored in the calculation of the number of neighbors
        @type offset: int
        """
        _AbstractHSExposure.__init__(self, model, radius, offset,
                'EXP_HSE_B_U', 'EXP_HSE_B_D')

    def _get_cb(self, r1, r2, r3):
        """Method to calculate CB-CA vector.

        @param r1, r2, r3: three consecutive residues (only r2 is used)
        @type r1, r2, r3: L{Residue}
        """
        if r2.get_resname() == 'GLY':
            gly_cb_vector= self._get_gly_cb_vector(r2)
            if gly_cb_vector is None:
              return None
            else:
              return gly_cb_vector, 0.0
        else:
            if r2.has_id('CB') and r2.has_id('CA'):
                vcb = r2['CB'].get_vector()
                vca = r2['CA'].get_vector()
                return (vcb - vca), 0.0
        return None


class ExposureCN(AbstractPropertyMap):
    """Residue exposure as number of CA atoms around its CA atom."""

    def __init__(self, model, radius=12.0, offset=0):
        """Initialize.

        A residue's exposure is defined as the number of CA atoms around
        that residues CA atom. A dictionary is returned that uses a L{Residue}
        object as key, and the residue exposure as corresponding value.

        @param model: the model that contains the residues
        @type model: L{Model}

        @param radius: radius of the sphere (centred at the CA atom)
        @type radius: float

        @param offset: number of flanking residues that are ignored in the calculation of the number of neighbors
        @type offset: int

        """
        assert(offset >= 0)
        ppb = CaPPBuilder()
        ppl = ppb.build_peptides(model)
        fs_map = {}
        fs_list = []
        fs_keys = []
        for pp1 in ppl:
            for i in range(0, len(pp1)):
                fs = 0
                r1 = pp1[i]
                if not is_aa(r1) or not r1.has_id('CA'):
                    continue
                ca1 = r1['CA']
                for pp2 in ppl:
                    for j in range(0, len(pp2)):
                        if pp1 is pp2 and abs(i - j) <= offset:
                            continue
                        r2 = pp2[j]
                        if not is_aa(r2) or not r2.has_id('CA'):
                            continue
                        ca2 = r2['CA']
                        d = (ca2 - ca1)
                        if d < radius:
                            fs += 1
                res_id = r1.get_id()
                chain_id = r1.get_parent().get_id()
                # Fill the 3 data structures
                fs_map[(chain_id, res_id)] = fs
                fs_list.append((r1, fs))
                fs_keys.append((chain_id, res_id))
                # Add to xtra
                r1.xtra['EXP_CN'] = fs
        AbstractPropertyMap.__init__(self, fs_map, fs_keys, fs_list)
