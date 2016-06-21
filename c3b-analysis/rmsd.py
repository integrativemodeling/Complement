from __future__ import print_function
import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import IMP.atom
import glob

rmfs=glob.glob("kmeans_weight_0_500_1//cluster.0/*.rmf3")


selection_dictionary={"beta_ANT_MG6beta":[(1,645,"beta"),(727,745,"alpha"),(746,805,"alpha")],
                  "ANA":[(650,726,"alpha")],
                  "MG7":[(806,911,"alpha")],
                  "CUBf_CUBg":[(912,962,"alpha"),(1269,1330,"alpha")],
                  "TED":[(963,1268,"alpha")],
                  "MG8":[(1331,1474,"alpha")],
                  "Anchor_C345C":[(1475,1495,"alpha"),(1496,1641,"alpha")]}

model=IMP.Model()


frames=[0]*len(rmfs)

model=IMP.Model()
pr=IMP.pmi.analysis.Precision(model,1,
                              selection_dictionary=selection_dictionary)
#pr.set_precision_style('pairwise_drmsd_k')

pr.add_structures(zip(rmfs,frames), 'all')



refrmf='../c3b-native/native.rmf3'
pr.set_reference_structure(refrmf,0)

print(pr.get_rmsd_wrt_reference_structure_with_alignment('all',"beta_ANT_MG6beta"))