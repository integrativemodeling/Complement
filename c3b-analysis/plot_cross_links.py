import IMP
import IMP.pmi
import IMP.pmi.output
import IMP.pmi.io
import IMP.pmi.io.crosslink
import operator
import sys

stat_file=sys.argv[1] #="kmeans_weight_0_500_10/cluster.3/stat.out"
output_file=sys.argv[2] #=Distances


po=IMP.pmi.output.ProcessOutput(stat_file)
#po=IMP.pmi.output.ProcessOutput("../../../AKAP79_Calmodulin_modeling/output/stat.0.out")
l0=[p for p in po.get_keys() if "CrossLinkingMassSpectrometryRestraint_Distance" in p]

d0=po.get_fields(l0)

from IMP.pmi.io.crosslink import FilterOperator as FO
cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Peptide Chain 1")
cldbkc.set_protein2_key("Peptide Chain 2")
cldbkc.set_residue1_key("Residue 1")
cldbkc.set_residue2_key("Residue 2")
cldbkc.set_psi_key("PSI")

cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file("../data/QCLMS_iC3-Domain-Architecture_Rappsilber_TableS2-1.csv")
cldb.create_new_keyword(cldb.state_key)
cldb.set_value(cldb.state_key,0)

fo=FO('Quantitation summary',operator.eq,'C3b-iC3 shared')| \
   FO('Quantitation summary',operator.eq,'C3b-C3 shared')| \
   FO('Quantitation summary',operator.eq,'Common')| \
   FO('Quantitation summary',operator.eq,'C3b Unique')

filtered_cldb=cldb.filter(fo)

tuplelist=[]
for kl in l0:
    tuplelist.append([[],"None","None"])
    tuplelist[-1][0]=map(float,d0[kl])
    for xl in filtered_cldb:
        if filtered_cldb.get_short_cross_link_string(xl) in kl:
           uid=xl[filtered_cldb.unique_sub_id_key]
           r1=str(xl[filtered_cldb.residue1_key])
           r2=str(xl[filtered_cldb.residue2_key])
           tuplelist[-1][1]=uid+"|"+r1+"|"+r2
           tuplelist[-1][2]=float(uid)

tuplelist.sort(key=lambda x: x[2])

IMP.pmi.output.plot_fields_box_plots(output_file,[p[0] for p in tuplelist]
                                                  ,range(len([p[0] for p in tuplelist])),frequencies=True,xlabels=[p[1] for p in tuplelist], scale_plot_length=0.3)

#IMP.pmi.output.plot_fields_box_plots("Distances",dists1
#                                                  ,range(len(dists1)),frequencies=True,xlabels=xlabels1, scale_plot_length=0.3)
