import IMP
import IMP.pmi
import IMP.pmi.analysis
import IMP.pmi.io
import os

model = IMP.Model()

native="../c3-native/native.rmf3"
reference="kmeans_weight_0_500_1/cluster.0/0.rmf3"


native_hiers,rs=IMP.pmi.analysis.get_hiers_and_restraints_from_rmf(model,0,native)
native_hier=native_hiers[0]
reference_hiers=IMP.pmi.analysis.get_hiers_from_rmf(model,0,reference)
reference_hier=reference_hiers[0]

s=IMP.atom.Selection(native_hier,molecule="beta")
ps=s.get_selected_particles()
xyz_native=[IMP.core.XYZ(p).get_coordinates() for p in ps]


s=IMP.atom.Selection(reference_hier,molecule="beta")
ps=s.get_selected_particles()
xyz_reference=[IMP.core.XYZ(p).get_coordinates() for p in ps]

transformation=IMP.algebra.IMP.algebra.get_transformation_aligning_first_to_second(xyz_native,xyz_reference)

rbs = set()
for p in IMP.atom.get_leaves(native_hier):
    if not IMP.core.XYZR.get_is_setup(p):
        IMP.core.XYZR.setup_particle(p)
        IMP.core.XYZR(p).set_radius(0.0001)
        IMP.core.XYZR(p).set_coordinates((0, 0, 0))

    if IMP.core.RigidBodyMember.get_is_setup(p):
        rb = IMP.core.RigidBodyMember(p).get_rigid_body()
        rbs.add(rb)
    else:
        IMP.core.transform(IMP.core.XYZ(p),
                           transformation)
for rb in rbs:
    IMP.core.transform(rb,transformation)

o=IMP.pmi.output.Output()
out_pdb_fn=os.path.join("native.pdb")
out_rmf_fn=os.path.join("native.rmf3")
o.init_pdb(out_pdb_fn,native_hier)
o.write_pdb(out_pdb_fn)

o.init_rmf(out_rmf_fn,[native_hier],rs)
o.write_rmf(out_rmf_fn)
o.close_rmf(out_rmf_fn)
