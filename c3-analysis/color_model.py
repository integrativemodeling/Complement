
# I have to deprecate this function
import IMP
import IMP.atom
import IMP.core
import IMP.display
import IMP.rmf
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.output
import RMF



model=IMP.Model()
rh = RMF.open_rmf_file_read_only('kmeans_weight_0_500_1/cluster.0/52.rmf3')
prots = IMP.rmf.create_hierarchies(rh, model)
rs = IMP.rmf.create_restraints(rh, model)
prot=prots[0]
IMP.rmf.load_frame(rh, 0)
model.update()

reduced_density_dict={"beta_ANT_MG6beta":[(1,645,"beta"),(727,745,"alpha"),(746,805,"alpha")],
                  "ANA":[(650,726,"alpha")],
                  "MG7":[(806,911,"alpha")],
                  "CUBf_CUBg":[(912,962,"alpha"),(1269,1330,"alpha")],
                  "TED":[(963,1268,"alpha")],
                  "MG8":[(1331,1474,"alpha")],
                  "Anchor_C345C":[(1475,1495,"alpha"),(1496,1641,"alpha")]}


color_dict={"beta_ANT_MG6beta":(0.3,0.3,0.3),
                  "ANA":(0.5,0.5,1.0),
                  "MG7":(1.0,0.87,0.37),
                  "CUBf_CUBg":(0.1,0.6,0.6),
                  "TED":(0.75,1.0,0.25),
                  "MG8":(0.555,0.222,0.111),
                  "Anchor_C345C":(1.0,0.5,1.0)}

for p in reduced_density_dict:
    for l in reduced_density_dict[p]:
      s=IMP.atom.Selection(prots,molecule=l[2],residue_indexes=range(l[0],l[1]+1))
      psel=s.get_selected_particles()
      c=color_dict[p]
      color=IMP.display.Color(c[0],c[1],c[2])
      for part in psel:
        IMP.display.Colored(part).set_color(color)



o=IMP.pmi.output.Output()
o.init_rmf("colored.rmf3",[prots[0]],rs)
o.write_rmf("colored.rmf3")
o.close_rmf("colored.rmf3")
