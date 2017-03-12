#!/usr/bin/python
from pdb_utils_noclass import *
from pyrosetta.toolbox import pose_from_rcsb
import random
from math import *
from xyzMath import *

def ran_range(frange):
    return random.uniform(-frange,+frange)


def get_helix_params(res1,res2):
#    print res1
    stub1 = stub(res1[0],res1[1],res1[2])
    stub2 = stub(res2[0],res2[1],res2[2])
    xform = stub2 * ~stub1
    axis, ang, cen = xform.rotation_axis_center()
    translation_along_axis_of_rotation = axis.dot(xform.t)
    radius_from_CA_1 = projperp(axis, cen - res1[1] ).length()
    print translation_along_axis_of_rotation,radius_from_CA_1,ang
    return

init_pyrosetta()
p=Pose()
rosetta.core.pose.make_pose_from_sequence(p, "AAAAAAAA","fa_standard")

Nsamples=10
omega=180.
frange=30.
for i in range(Nsamples):
    phi1=-120. + ran_range(frange)
    psi1=125. + ran_range(frange)
    phi2=-120. + ran_range(frange)
    psi2=125. + ran_range(frange)
    for j in [1,3,5,7]:
        p.set_phi(j,phi1)
        p.set_psi(j,psi1)
        p.set_omega(j,180.)
        p.set_phi(j+1,phi2)
        p.set_psi(j+1,psi2)
        p.set_omega(j+1,180.)
    p.dump_pdb('test_%s.pdb'%i)
    nres=p.size()
  #  for i in range(nres):
        #print(i+1,p.phi(i+1),p.psi(i+1),p.omega(i+1))    
    res1=[Vec(p.residue(1).xyz("N")),Vec(p.residue(1).xyz("CA")),Vec(p.residue(1).xyz("C"))]
    res2=[Vec(p.residue(3).xyz("N")),Vec(p.residue(3).xyz("CA")),Vec(p.residue(3).xyz("C"))]
#    res2=p.residue(3).xyz("CA")
    print('helix params 1-2, 2-3, and 1-3')
    get_helix_params(res1,res2)
#    res3=[Vec(p.residue(5).xyz("N")),Vec(p.residue(5).xyz("CA")),Vec(p.residue(5).xyz("C"))]
#    get_helix_params(res2,res3)
#    get_helix_params(res1,res3)
#arm = pose_from_rcsb("5AEI")
p=rosetta.core.import_pose.pose_from_file("DHR7.pdb")
res1=[Vec(p.residue(50).xyz("N")),Vec(p.residue(50).xyz("CA")),Vec(p.residue(50).xyz("C"))]
res2=[Vec(p.residue(92).xyz("N")),Vec(p.residue(92).xyz("CA")),Vec(p.residue(92).xyz("C"))]
res3=[Vec(p.residue(134).xyz("N")),Vec(p.residue(134).xyz("CA")),Vec(p.residue(134).xyz("C"))]
print('DHR7 params ')
get_helix_params(res1,res2)
get_helix_params(res2,res3)

p=rosetta.core.import_pose.pose_from_file("ARM_pept.pdb")
#nres=ARM_pept.size()
res1=[Vec(p.residue(3).xyz("N")),Vec(p.residue(3).xyz("CA")),Vec(p.residue(3).xyz("C"))]
res2=[Vec(p.residue(5).xyz("N")),Vec(p.residue(5).xyz("CA")),Vec(p.residue(5).xyz("C"))]
print('ARM peptide params ')
get_helix_params(res1,res2)

p=rosetta.core.import_pose.pose_from_file("5AEI_A.pdb")
print(p.residue(106).name(),p.residue(107).name(),p.residue(149).name(),p.residue(191).name())
res1=[Vec(p.residue(107).xyz("N")),Vec(p.residue(107).xyz("CA")),Vec(p.residue(107).xyz("C"))]
res2=[Vec(p.residue(149).xyz("N")),Vec(p.residue(149).xyz("CA")),Vec(p.residue(149).xyz("C"))]
res3=[Vec(p.residue(191).xyz("N")),Vec(p.residue(191).xyz("CA")),Vec(p.residue(191).xyz("C"))]
print('ARM protein params ')
get_helix_params(res1,res2)
get_helix_params(res2,res3)
#for i in range(nres):
#    print(i+1,ARM_pept.phi(i+1),ARM_pept.psi(i+1),ARM_pept.omega(i+1))


