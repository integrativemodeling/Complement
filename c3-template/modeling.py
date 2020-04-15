import IMP
import IMP.core
import IMP.algebra
import IMP.atom
import IMP.container

import IMP.pmi1.mmcif
import ihm.location
import ihm.model
import ihm.cross_linkers
try:
    import ihm.reference
except ImportError:
    pass
import IMP.pmi1.restraints.crosslinking_new
import IMP.pmi1.restraints.stereochemistry
import IMP.pmi1.restraints.em
import IMP.pmi1.restraints.basic
import IMP.pmi1.representation
import IMP.pmi1.tools
import IMP.pmi1.samplers
import IMP.pmi1.output
import IMP.pmi1.macros
import IMP.pmi1.io

import os
import sys
import operator


'''

Suggested rigid body
Beta-chain
ANA
MG7
CUB
TED
MG8
C345C

Domain	start	end	chain	Domain				
MG1	    1	    105	    A	MG1				
MG2	    106	    205	    B	MG2				
MG3	    206	    330	    C	MG3				
MG4	    331	    440	    D	MG4				
MG5	    441	    534	    E	MG5				
MG6alpha	535	577	    F	MG6alpha				
LNK	    578	    645	    G	LNK				
ANA	    650	    726	    H	ANA	Cleaved in C3b, however retained in iC3			
ANT	    727	    745	    I	ANT				
MG6beta	746	    805	    J	MG6beta				
MG7	    806	    911	    K	MG7				
CUBf	912	    962	    L	CUBf				
TED	    963	    1268	M	TED				
CUBg	1269	1330	N	CUBg				
MG8	    1331	1474	O	MG8				
Anchor	1475	1495	P	Anchor				
C345C	1496	1641	Q	C345C				


'''




# setting up parameters

rbmaxtrans = 4.00
fbmaxtrans = 5.00
rbmaxrot=0.03
outputobjects = []
sampleobjects = []

# setting up topology

m = IMP.Model()
simo1 = IMP.pmi1.representation.Representation(m,upperharmonic=True,disorderedlength=False)
simo1.state.short_name = 'C3'
simo1.state.long_name = 'Human complement component C3'

fastadirectory="../data/"
pdbdirectory="../data/"
xlmsdirectory="../data/"

compactrepresentation=False

if '--mmcif' in sys.argv:
    # Record the modeling protocol to an mmCIF file
    po = IMP.pmi1.mmcif.ProtocolOutput(open('complement.cif', 'w'))
    simo1.add_protocol_output(po)
    po.system.title = ('Structure of complement C3(H2O) revealed by '
                       'quantitative cross-linking/mass spectrometry '
                       'and modeling')
    # Add publication
    po.system.citations.append(ihm.Citation.from_pubmed_id(27250206))

    # Point to repositories where files are deposited
    zenodo_id = '1285940'
    doi = '10.5281/zenodo.' + zenodo_id
    url_top = 'https://zenodo.org/record/%s/files' % zenodo_id
    simo1.add_metadata(ihm.location.Repository(
                       doi=doi, root="..",
                       url="%s/integrativemodeling/Complement-v1.0.zip"
                                            % url_top,
                       top_directory="integrativemodeling-Complement-a6a1494"))


simo1.dry_run = '--dry-run' in sys.argv

if compactrepresentation:
       # compname  hier_name    color         fastafile              fastaid          pdbname      chain    resrange      read    "BEADS"ize rigid_body super_rigid_body emnum_components emtxtfilename  emmrcfilename chain of super rigid bodies
    domains=   [
    ("beta",  "beta",    0.0,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3|C3b_beta_1-645", pdbdirectory+"/2A73.pdb" ,   "A",   (1,645,0),  None,        10,      0,         [0],     0,   None,  None,   [0]),
    ("alpha",  "ANA",    0.10,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (650,726,0),  None,        10,      1,         [1],     0,   None,  None,   [1]),
    ("alpha",  "ANT",    0.20,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (727,745,0),  None,        10,      0,         [0],     0,   None,  None,   [1]),
    ("alpha",  "MG6beta",    0.30,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (746,805,0),  None,        10,  0,         [0],     0,   None,  None,   [1]),
    ("alpha",  "MG7",    0.40,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (806,911,0),  None,        10,      2,         [1],     0,   None,  None,   [1]),
    ("alpha",  "CUBf",    0.50,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (912,962,0),  None,        10,     3,         [1],     0,   None,  None,   [1]),
    ("alpha",  "TED",    0.60,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (963,1268,0),  None,        10,     4,         [1],     0,   None,  None,   [1]),
    ("alpha",  "CUBg",    0.70,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (1269,1330,0),  None,        10,   3,         [1],     0,   None,  None,   [1]),
    ("alpha",  "MG8",    0.80,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (1331,1474,0),  None,        10,    5,         [1],     0,   None,  None,   [1]),
    ("alpha",  "Anchor",    0.90,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (1475,1495,0),  None,        10, 6,         [1],     0,   None,  None,   [1]),
    ("alpha",  "C345C",    1.0,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (1496,1641,0),  None,        10,  6,         [1],     0,   None,  None,   [1]),
    ]

else:
     domains=   [
    ("beta",  "beta",    0.0,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3|C3b_beta_1-645", pdbdirectory+"/2A73.pdb" ,   "A",       (1,645,0),  None,       1,     0,         [0],     0,   None,  None,   [1]),
    ("alpha",  "ANA",    0.10,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",      (650,726,0),  None,     1,     1,         [1],     0,   None,  None,   [1]),
    ("alpha",  "ANT",    0.20,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",      (727,745,0),  None,     1,     0,         [0],     0,   None,  None,   [1]),
    ("alpha",  "MG6beta",    0.30,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",  (746,804,0),  None,     1,     0,         [0],     0,   None,  None,   [1]),
    ("alpha",  "MG6beta_MG7_Link",    0.30,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", "BEADS" ,   None,         (805,806,0),  None,     1,     None,         [0],     0,   None,  None,   [1]),
    ("alpha",  "MG7",    0.40,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",      (807,910,0),  None,     1,     2,         [1],     0,   None,  None,   [1]),
    ("alpha",  "MG7_CUBf_Link",    0.40,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", "BEADS" ,   None,            (911,912,0),  None,     1,     None,         [1],     0,   None,  None,   [1]),
    ("alpha",  "CUBf",    0.50,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",     (913,961,0),  None,     1,     3,         [1],     0,   None,  None,   [1]),
    ("alpha",  "CUBf_TED_Link",    0.50,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", "BEADS" ,   None,            (962,963,0),  None,     1,     None,         [1],     0,   None,  None,   [1]),
    ("alpha",  "TED",    0.60,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",      (964,1267,0),  None,    1,     4,         [1],     0,   None,  None,   [1]),
    ("alpha",  "TED_CUBg_Link",    0.60,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", "BEADS" ,   None,            (1268,1269,0),  None,   1,     None,         [1],     0,   None,  None,   [1]),
    ("alpha",  "CUBg",    0.70,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",     (1270,1329,0),  None,   1,     3,         [1],     0,   None,  None,   [1]),
    ("alpha",  "CUBg_MG8_Link",    0.70,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", "BEADS" ,   None,            (1330,1331,0),  None,   1,     None,         [1],     0,   None,  None,   [1]),
    ("alpha",  "MG8",    0.80,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",      (1332,1473,0),  None,   1,     5,         [1],     0,   None,  None,   [1]),
    ("alpha",  "MG8_Anchor_Link",    0.80,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", "BEADS" ,   None,          (1474,1475,0),  None,   1,     None,         [1],     0,   None,  None,   [1]),
    ("alpha",  "Anchor_C345C",    1.0,  fastadirectory+"/C3-iC3-C3b_sequence.fasta",  "C3|iC3_full_alpha", pdbdirectory+"/2A73.pdb" ,   "B",   (1476,1641,0),  None,   1,     6,         [1],     0,   None,  None,   [1]),
    ]   

bm1=IMP.pmi1.macros.BuildModel1(simo1)
bm1.build_model(domains,sequence_connectivity_scale=1.0,sequence_connectivity_resolution=1)
bm1.scale_bead_radii(40,0.8)
#bm1.save_rmf("simo1.rmf3")
bm1.set_coordinates("MG6beta_MG7_Link",IMP.pmi1.tools.get_terminal_residue_position(simo1,bm1.domain_dict["MG6beta"], terminus="C", resolution=1))
bm1.set_coordinates("MG7_CUBf_Link",IMP.pmi1.tools.get_terminal_residue_position(simo1,bm1.domain_dict["MG7"], terminus="C", resolution=1))
bm1.set_coordinates("CUBf_TED_Link",IMP.pmi1.tools.get_terminal_residue_position(simo1,bm1.domain_dict["CUBf"], terminus="C", resolution=1))
bm1.set_coordinates("TED_CUBg_Link",IMP.pmi1.tools.get_terminal_residue_position(simo1,bm1.domain_dict["TED"], terminus="C", resolution=1))
bm1.set_coordinates("CUBg_MG8_Link",IMP.pmi1.tools.get_terminal_residue_position(simo1,bm1.domain_dict["CUBg"], terminus="C", resolution=1))
bm1.set_coordinates("MG8_Anchor_Link",IMP.pmi1.tools.get_terminal_residue_position(simo1,bm1.domain_dict["MG8"], terminus="C", resolution=1))


# randomize the initial configuration

# simo1.shuffle_configuration(100)


simo1.set_rigid_bodies_max_rot(rbmaxrot)
simo1.set_floppy_bodies_max_trans(fbmaxtrans)
simo1.set_rigid_bodies_max_trans(rbmaxtrans)

outputobjects.append(simo1)
sampleobjects.append(simo1)

simo1.show_component_table("alpha")
simo1.show_component_table("beta")


# scoring function

ev = IMP.pmi1.restraints.stereochemistry.ExcludedVolumeSphere(simo1,resolution=10)
ev.set_weight(1.0)
ev.add_to_model()
outputobjects.append(ev)


from IMP.pmi1.io.crosslink import FilterOperator as FO
cldbkc=IMP.pmi1.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Peptide Chain 1")
cldbkc.set_protein2_key("Peptide Chain 2")
cldbkc.set_residue1_key("Residue 1")
cldbkc.set_residue2_key("Residue 2")
cldbkc.set_psi_key("PSI")
cldb=IMP.pmi1.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file("../data/QCLMS_iC3-Domain-Architecture_Rappsilber_TableS2-1.csv")

fo=FO('Quantitation summary',operator.eq,'C3 Unique')| \
   FO('Quantitation summary',operator.eq,'C3-iC3 shared')|  \
   FO('Quantitation summary',operator.eq,'Common')| \
   FO('Quantitation summary',operator.eq,'C3b-C3 shared')

filtered_cldb=cldb.filter(fo)

xl = IMP.pmi1.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(representation=simo1,
                            CrossLinkDataBase=filtered_cldb,
                            length=21.0,
                            slope=0.01,
                            resolution=1.0,
                            label="XL")
# Point to the raw mass spec data and peaklists used to derive the crosslinks.
l = ihm.location.PRIDELocation('PXD003486',
                         details='All raw mass spectrometry files and '
                                 'peaklists used in the study')
xl.dataset.add_primary(ihm.dataset.MassSpecDataset(location=l))

xl.add_to_model()
sigma=xl.sigma_dictionary["SIGMA"][0]
psic=xl.psi_dictionary["Confident"][0]
psiq=xl.psi_dictionary["Questionable"][0]
sigma.set_scale(11.0)
psic.set_scale(0.01)
psiq.set_scale(0.01)
xl.set_psi_is_sampled(True)
xl.set_sigma_is_sampled(True)
sampleobjects.append(xl)
outputobjects.append(xl)

'''
xld = IMP.pmi1.restraints.crosslinking_new.DisulfideCrossLinkRestraint(
    simo1,
    (851,851,"alpha"),
    (1491,1491,"alpha"),
    label="DisulfideBond1",
    resolution=1,
    slope=0.05)
xld.add_to_model()
outputobjects.append(xld)
'''

#setting up the disulfide bond
xld=IMP.pmi1.restraints.basic.DistanceRestraint(simo1,
                                               (851,851,"alpha"),
                                               (1491,1491,"alpha"),
                                               4.0,
                                               7.0)
xld.add_to_model()
outputobjects.append(xld)


mc1=IMP.pmi1.macros.ReplicaExchange0(m,
                                    simo1,
                                    monte_carlo_sample_objects=sampleobjects,
                                    output_objects=outputobjects,
                                    crosslink_restraints=[xl,xld],
                                    monte_carlo_temperature=1.0,
                 		            simulated_annealing=False,
                                    simulated_annealing_minimum_temperature=1.0,
                                    simulated_annealing_maximum_temperature=2.5,
                                    simulated_annealing_minimum_temperature_nframes=100,
                                    simulated_annealing_maximum_temperature_nframes=10,
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=2.5,
                                    number_of_best_scoring_models=50,
                                    monte_carlo_steps=10,
                                    number_of_frames=100000,
                                    write_initial_rmf=True,
                                    initial_rmf_name_suffix="initial",
                                    stat_file_name_suffix="stat",
                                    best_pdb_name_suffix="model",
                                    do_clean_first=True,
                                    do_create_directories=True,
                                    global_output_directory="output",
                                    rmf_dir="rmfs/",
                                    best_pdb_dir="pdbs/",
                                    replica_stat_file_suffix="stat_replica",
                                    test_mode=simo1.dry_run)
mc1.execute_macro()

if '--mmcif' in sys.argv:
    # Link entities to UniProt
    if hasattr(ihm, 'reference'):
        lpep = ihm.LPeptideAlphabet()
        sd_beta = [ihm.reference.SeqDif(292, lpep['P'], lpep['L'])]
        for subunit, accession, db_align_begin, seq_dif in (
                ('beta', 'P01024', 23, sd_beta), ('alpha', 'P01024', 672, [])):
            ref = ihm.reference.UniProtSequence.from_accession(accession)
            ref.alignments.append(ihm.reference.Alignment(
                db_begin=db_align_begin, seq_dif=seq_dif))
            e = po.asym_units[subunit].entity.references.append(ref)

    # Add clustering info to the mmCIF file
    os.chdir('../c3-analysis')
    loc = ihm.location.WorkflowFileLocation('clustering.py',
                      details='Clustering and analysis script for C3 state')
    simo1.add_metadata(loc)
    with open('clustering.py') as fh:
        exec(fh.read())

    # Add info on the C3b and iC3 state modeling and clustering
    for state in ('c3b', 'ic3'):
        os.chdir('../%s-template' % state)
        loc = ihm.location.WorkflowFileLocation('modeling.py',
                      details='Main script for %s state modeling' % state)
        simo1.add_metadata(loc)
        with open('modeling.py') as fh:
            exec(fh.read())

        os.chdir('../%s-analysis' % state)
        loc = ihm.location.WorkflowFileLocation('clustering.py',
                 details='Clustering and analysis script for %s state' % state)
        simo1.add_metadata(loc)
        with open('clustering.py') as fh:
            exec(fh.read())

    # Correct # of models (we ran two independent runs)
    for state_group in po.system.state_groups:
        for state in state_group:
            for model_group in state:
                for model in model_group:
                    p = model.protocol
                    p.steps[-1].num_models_end = 200000
                    p.analyses[0].steps[0].num_models_begin = 200000

    # Correct crosslinker type (defaults to restraint label, which is "XL")
    for r in po.system.restraints:
        if hasattr(r, 'linker'):
            r.linker = ihm.cross_linkers.bs3

    # Add reaction cycle between states (Fig 1B)
    c3, c3b, ic3 = po.system.state_groups[0]
    proc = ihm.model.OrderedProcess('steps in a reaction pathway')
    po.system.ordered_processes.append(proc)
    s = ihm.model.ProcessStep(
                    description='proteolytic cleavage of C3 yields C3b')
    for begin_model_group in c3:
        for end_model_group in c3b:
            s.append(ihm.model.ProcessEdge(
                                begin_model_group, end_model_group))
    proc.steps.append(s)

    s = ihm.model.ProcessStep(
                    description='thioester hydrolysis of C3 yields iC3')
    for begin_model_group in c3:
        for end_model_group in ic3:
            s.append(ihm.model.ProcessEdge(begin_model_group, end_model_group))
    proc.steps.append(s)

    # End up in initial working directory
    os.chdir('../c3-template')
    po.flush()
