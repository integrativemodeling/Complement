import IMP
import IMP.pmi1
import IMP.pmi1.analysis
import IMP.pmi1.output
import IMP.atom
import glob
import itertools

is_mpi=False
test_mode=False  # runs on the first 10 structures to test if it runs smoothly

# specify the cluster directory to be analysed

root_cluster_directory='kmeans_weight_0_500_2'

# choose whatever selection for the precision calculation
selection_dictionary={"beta_ANT_MG6beta":[(1,645,"beta"),(727,745,"alpha"),(746,805,"alpha")],
                  "ANA":[(650,726,"alpha")],
                  "MG7":[(806,911,"alpha")],
                  "CUBf_CUBg":[(912,962,"alpha"),(1269,1330,"alpha")],
                  "TED":[(963,1268,"alpha")],
                  "MG8":[(1331,1474,"alpha")],
                  "Anchor_C345C":[(1475,1495,"alpha"),(1496,1641,"alpha")],
                  "beta":["beta"],
                  "alpha":["alpha"]}



#####################################################################
# don't change anything below
rmfs=[]
frames=[]
clusterdirectories=glob.glob(root_cluster_directory+'/cluster.*/')

if test_mode:
  # runs on the first 10 structures to test if it runs smoothly
  for clusterdirectory in clusterdirectories:
      rmfs.append(glob.glob(clusterdirectory+'/*.rmf3')[0::10])
      frames.append([0]*len(rmfs[-1]))
else:
  for clusterdirectory in clusterdirectories:
      rmfs.append(glob.glob(clusterdirectory+'/*.rmf3')[0::1])
      frames.append([0]*len(rmfs[-1]))
 
model=IMP.Model()
pr=IMP.pmi1.analysis.Precision(model,1,selection_dictionary=selection_dictionary)
pr.set_precision_style('pairwise_rmsd')

for n in range(len(rmfs)):
    pr.add_structures(zip(rmfs[n],frames[n]),clusterdirectories[n])


for pair in itertools.product(range(len(rmfs)), repeat=2):
    clus1=pair[0]
    clus2=pair[1]
    pr.get_precision(
                     clusterdirectories[clus1],
                     clusterdirectories[clus2],
                     outfile=root_cluster_directory+"/precision."+str(clus1)+"."+str(clus2)+".out",
                     skip=1,
                     selection_keywords=selection_dictionary.keys())

for n in range(len(rmfs)):
    pr.get_rmsf(clusterdirectories[n],clusterdirectories[n]+"/",skip=1,set_plot_yaxis_range=(0,20.0))
