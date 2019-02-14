import IMP
import IMP.pmi1
import IMP.pmi1.macros
import sys
import glob
import ihm.location

is_mpi=True

model=IMP.Model()

# initialize the macro








# fields that have to be extracted for the stat file

feature_list=["CrossLinkingMassSpectrometryRestraint_Distance_",
              "CrossLinkingMassSpectrometryRestraint_Score_",
              "ISDCrossLinkMS_Data_Score",
              "GaussianEMRestraint_None",
              "SimplifiedModel_Linker_Score_None",
              "ISDCrossLinkMS_Psi",
              "ISDCrossLinkMS_Sigma"]

# Dictionary of densities to be calculated
# the key is the name of the file and the value if the selection
# example:
#              {"med17-CTD":[(200,300,"med17")],"med17-CTD.med14":[(200,300,"med17"),"med14"]   }


reduced_density_dict_0={"beta_ANT_MG6beta":[(1,645,"beta"),(727,745,"alpha"),(746,805,"alpha")],
                  "ANA":[(650,726,"alpha")],
                  "MG7":[(806,911,"alpha")],
                  "CUBf_CUBg":[(912,962,"alpha"),(1269,1330,"alpha")],
                  "TED":[(963,1268,"alpha")],
                  "MG8":[(1331,1474,"alpha")],
                  "Anchor_C345C":[(1475,1495,"alpha"),(1496,1641,"alpha")]}



# list of component names needed to calculate the RMSD for the clustering





alignment_components_names={}
alignment_components_names["beta"]="beta"

rmsd_components_names={}
rmsd_components_names["alpha"]="alpha"
rmsd_components_names["beta"]="beta"

nclusters=2


mc=IMP.pmi1.macros.AnalysisReplicaExchange0(model,
                                        stat_file_name_suffix="stat",     # don't change
                                        merge_directories=["../ic3.1","../ic3.2"],
                                        global_output_directory="./output")  # don't change

if '--mmcif' in sys.argv:
    mc.test_mode = simo1.dry_run
    for po in simo1.protocol_output:
        mc.add_protocol_output(po)

                                      # number of clusters needed by kmeans
mc.clustering("SimplifiedModel_Total_Score_None",  # don't change, field where to find the score
              "rmf_file",                          # don't change, field where to find the path for the rmf_file
              "rmf_frame_index",                   # don't change, field for the frame index
              prefiltervalue=50.0,               # prefilter the models by score
              number_of_best_scoring_models=200,   # number of models to be clustered
              alignment_components=alignment_components_names,           # don't change, (list of proteins you want to use for structural alignment
              rmsd_calculation_components=rmsd_components_names, # list of proteins used to calculated the rmsd
              distance_matrix_file="distance.rawmatrix.pkl", # save the distance matrix
              outputdir="kmeans_weight_0_500_"+str(nclusters)+"/",  # directory name for the clustering
              feature_keys=feature_list,                     # extract these fields from the stat file
              state_number=0,
              first_and_last_frames=[0.0,1.0],
              load_distance_matrix_file=True,                # skip the matrix calcuklation and read the precalculated matrix
              skip_clustering=False,                         # skip clustering
              display_plot=False,                            # display the heat map plot of the distance matrix
              exit_after_display=False,                      # exit after having displayed the distance matrix plot
              get_every=1,                                   # skip structures for faster computation
              number_of_clusters=nclusters,                  # number of clusters to be used by kmeans algorithm
              voxel_size=2.0,                                # voxel size of the mrc files
              density_custom_ranges=reduced_density_dict_0)    # setup the list of densities to be calculated

if '--mmcif' in sys.argv:
    # Point to deposited ensembles in DCD format
    dcds = []
    for i in range(nclusters):
        r = ihm.location.Repository(doi=doi,
                                    url="%s/ic3-cluster%d.dcd" % (url_top, i))
        dcds.append(ihm.location.OutputFileLocation(path='.', repo=r,
                            details="All iC3 models in cluster %d" % (i+1)))
    for po in simo1.protocol_output:
        for i, dcd in enumerate(dcds):
            po.set_ensemble_file(i, dcd)
