patchDockTemplate='''########################################################################
#    PatchDock Parameter File for %(rPdbFname)s-%(lPdbFname)s
########################################################################

#    File Names:
receptorPdb %(rPdbFname)s
ligandPdb %(lPdbFname)s
receptorActiveSite %(patchDockWdir)s/rBindingSite.tab
ligandActiveSite %(patchDockWdir)s/lBindingSite.tab
#receptorBlockingSite site1.txt
#ligandBlockingSite site2.txt
protLib %(patchDockRootDir)s/chem.lib
log-file %(patchDockWdir)s/patch_dock.log
log-level 2

#   Distance constraint parameters:
#       distanceConstraints <rec_atom_index> <lig_atom_index> <dist_thr>
# <rec_atom_index> - receptor atom used for constraint
# <lig_atom_index> - ligand atom used for constraint
# <dist_thr> - maximum allowed distance between the specified atom centers
#distanceConstraints rec_atom_index lig_atom_index dist_thr
#pointDistanceConstraints rec_coord lig_coord min_dist max_dist
#distanceConstraintsFile file_name
#membraneLine rec_atom_index lig_atom_index

#    Surface Segmentation Parameters:
#       receptorSeg <low_patch_thr> <high_patch_thr> <prune_thr>
#                   <knob> <flat> <hole>
#                   <hot spot filter type>
#    <low_patch_thr>,<high_patch_thr> - min and max patch diameter
#    <prune_thr> - minimal distance between points inside the patch
#    <knob> <flat> <hole> - types of patches to dock (1-use, 0-do not use) (may need tuning)
#    <hot spot filter type> :None - 0, Antibody - 1, Antigen - 2
#                             Protease - 3, Inhibitor - 4, Drug - 5
receptorSeg 10.0 20.0 1.5 1 0 1 0
ligandSeg 10.0 20.0 1.5 1 0 1 0

#    Scoring Parameters:
#        scoreParams <small_interfaces_ratio> <max_penetration> <ns_thr>
#                    <rec_as_thr> <lig_as_thr> <patch_res_num> <w1 w2 w3 w4 w5>
#    <small_interfaces_ratio> - the ratio of the low scoring transforms to be removed
#    <max_penetration> - maximal allowed penetration between molecules surfaces
#    <ns_thr> - normal score threshold
#    <rec_as_thr> <lig_as_thr> - the minimal ratio of the active site area in the solutions
#    <patch_res_num> - number of results to consider in each patch
#    <w1 w2 w3 w4 w5> - scoring weights for ranges:
#                [-5.0,-3.6],[-3.6,-2.2],[-2.2,-1.0],[-1.0,1.0],[1.0-up]
scoreParams 0.3 -5.0 0.5 0.0 0.0 1500 -8 -4 0 1 0

#    Desolvation Scoring Parameters:
#        desolvationParams <energy_thr> <cut_off_ratio>
#    <energy_thr> - remove all results with desolvation energy higher than threshold
#    <cut_off_ratio> - the ratio of low energy results to be kept
#    First filtering with energy_thr is applied and the remaining results
#    can be further filtered with cut_off_ratio.
desolvationParams 500.0 1.0

########################################################################
#   Advanced Parameters
########################################################################

#    Clustering Parameters:
#    clusterParams < rotationVoxelSize > < discardClustersSmaller > < rmsd > < final clustering rmsd >
clusterParams 0.1 4 2.0 4.0

#    Base Parameters:
#    baseParams <min_base_dist> <max_base_dist> <# of patches for base: 1 or 2>
baseParams 4.0 13.0 2

#    Matching Parameters:
#  matchingParams <geo_dist_thr> <dist_thr> <angle_thr> <torsion_thr> 
#     <angle_sum_thr>
matchingParams 1.5 1.5 0.4 0.5 0.9
# 1 - PoseClustering (default), 2 - Geometring Hashing
matchAlgorithm 1

#    Grid Parameters:
#      receptorGrid <gridStep> <maxDistInDistFunction> <vol_func_radius>
#      NOTE: the vol_func_radius of small molecules and peptides should be 3A!
receptorGrid 0.5 6.0 6.0
ligandGrid 0.5 6.0 6.0

#    Connolly Surface Parameters:
#      receptorMs <density> <probe_radius>
receptorMs 10.0 1.8
ligandMs 10.0 1.8
'''
