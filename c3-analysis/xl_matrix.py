import IMP.pmi1
import IMP.pmi1.io
import IMP.pmi1.io.crosslink
import IMP.pmi1.io.xltable

import operator

from IMP.pmi1.io.crosslink import FilterOperator as FO
cldbkc=IMP.pmi1.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Peptide Chain 1")
cldbkc.set_protein2_key("Peptide Chain 2")
cldbkc.set_residue1_key("Residue 1")
cldbkc.set_residue2_key("Residue 2")
cldbkc.set_psi_key("PSI")

cldb=IMP.pmi1.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file("../data/QCLMS_iC3-Domain-Architecture_Rappsilber_TableS2-1.csv")
cldb.create_new_keyword(cldb.state_key)
cldb.set_value(cldb.state_key,0)

fo=FO('Quantitation summary',operator.eq,'C3 Unique')| \
   FO('Quantitation summary',operator.eq,'C3-iC3 shared')|  \
   FO('Quantitation summary',operator.eq,'Common')| \
   FO('Quantitation summary',operator.eq,'C3b-C3 shared')

filtered_cldb=cldb.filter(fo)

xlt=IMP.pmi1.io.xltable.XLTable(35)
prots = ["alpha","beta"]



xlt.load_sequence_from_fasta_file("../data/C3-iC3-C3b_sequence.fasta","C3|iC3_full","alpha")
xlt.load_sequence_from_fasta_file("../data/C3-iC3-C3b_sequence.fasta","C3|iC3_full","beta")
xlt.load_rmf_coordinates("kmeans_weight_0_500_3/cluster.1/28.rmf3",0,prots)
xlt.load_crosslinks(filtered_cldb)
xlt.setup_contact_map()

xlt.plot_table(prot_listx=prots,
   prot_listy=prots,
   alphablend=0.4,
   scale_symbol_size=1.1,
   gap_between_components=50,
   filename="XL_table.pdf",
   contactmap=True,
   crosslink_threshold=35.0,
   display_residue_pairs=False)
