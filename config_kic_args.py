#!/work/sheffler/anaconda2/bin/ipython 
from __future__ import print_function

from pyrosetta import *
from rosetta import *
from itertools import chain


# ideal values for the peptide bond (taken from the alanine params file)
_peptide_bond_params = {'atom1': 'C', 'atom2': 'N', 'bond_length': 1.328685,
                        'angle1': 116.199993, 'angle2': 121.699997,
                        'torsion': 180.}


_DEBUG = False  # global debug setting


def configure_peptide_stub_mover(residues, resNames='ALA'):
    if type(resNames) != list:
        resNames = [resNames] * len(residues)
    assert(len(residues) == len(resNames))

    psm = protocols.cyclic_peptide.PeptideStubMover()
    for resNo, resName in zip(residues, resNames):
        anchorRes, stub_mode = (resNo - 1, 'Insert') if resNo != residues[-1] \
            else (resNo, 'Prepend')
        if _DEBUG:
            print('PeptideStubMover: {}  {}'.format(stub_mode, anchorRes))
        psm.add_residue(stub_mode, resName, 0, False, '', 1, anchorRes, '')
    return psm


def configure_set_torsion_mover(residues, torsionNames='omega'):
    # Next we need to set the newly-added omegas to 180 degrees.
    if type(torsionNames) != list:
        torsionNames = [torsionNames] * len(residues)
    assert(len(residues) == len(torsionNames))

    st = protocols.simple_moves.SetTorsion()
    # The configuration of this mover interacts with private data, so it can't
    # be set up through PyRosetta. Fortunately, we can set the torsions
    # directly in a for loop, or create a custom, simple Mover to do this step
    '''
    for resNo, torsionName in zip(residues, torsionNames):
        torsionRes = resNo - 1 if resNo != residues[-1] else resNo
        st.residues_.push_back(torsionRes)
        st.torsion_name_.push_back(torsionName)
        st.angle_.push_back('180.')  # weird, but yes this should be a string
    '''
    return st


def set_omegas_to_180(p, residues):
    for resNo in residues:
        torsionRes = resNo - 1 if resNo != residues[-1] else resNo
        if _DEBUG:
            print('SetTorsion: {}  {}'.format('omega', torsionRes))
        p.set_omega(torsionRes, 180.)


def declare_bond(new_bond):
    # Now we need to tell Rosetta that there's a peptide bond between the last
    # of the three residues that we appended and the residuethat we prepended.
    # Note that this bond starts out ridiculously stretched -- it's just a line
    # in space between two atoms that are in entirely the wrong place.
    # GenKIC will later give us good bond geometry:
    db = protocols.cyclic_peptide.DeclareBond()
    if _DEBUG:
        print('DeclareBond: {}  {}'.format(residues[-2], residues[-1]))
    db.set(new_bond[0], _peptide_bond_params['atom1'],
           new_bond[1], _peptide_bond_params['atom2'], True)
    return db


def configure_gen_kic_mover(residues, new_bond, excluded, scorefxn):
    if type(excluded) != list:
        excluded = [excluded]
    pivot_res = [resNo for resNo in residues if resNo not in excluded]

    gk = protocols.generalized_kinematic_closure.GeneralizedKIC()
    gk.set_selector_type('lowest_energy_selector')
    gk.set_selector_scorefunction(scorefxn)

    gk.set_closure_attempts(500)
    gk.set_min_solution_count(5)
    gk.set_ntries_before_giving_up(500)
    # gk.set_preselection_mover(NAME_OF_MOVER)

    gk.add_perturber('randomize_backbone_by_rama_prepro')
    for resNo in residues:
        gk.add_loop_residue(resNo)
        gk.add_residue_to_perturber_residue_list(resNo)

        if resNo not in pivot_res:
            continue

        # The pivot residues are not necessarily in good regions of
        # Ramachandran space, so we should filter by the rama_prepro energy of
        # pivot positions.
        gk.add_filter('rama_prepro_check')
        gk.set_filter_resnum(resNo)
        gk.set_filter_rama_cutoff_energy(2.0)

    gk.add_filter('loop_bump_check')

    args = list(chain.from_iterable(zip(pivot_res, ['CA'] * len(pivot_res))))
    if _DEBUG:
        print('GeneralizedKIC: {}'.format(args))

    gk.set_pivot_atoms(*args)

    gk.close_bond(new_bond[0], _peptide_bond_params['atom1'],
                  new_bond[1], _peptide_bond_params['atom2'],
                  0, '', 0, '',  # optional params -- use default values
                  _peptide_bond_params['bond_length'],
                  _peptide_bond_params['angle1'],
                  _peptide_bond_params['angle2'],
                  _peptide_bond_params['torsion'], True, False)
    return gk

def configure_and_apply_movers(orig_pose, scorefxn, insert_res, KIC_res, new_bond,excluded_from_pivot):
    attempts = 0
    while True:
        p = orig_pose.clone()
 
        # configure movers in the loop to erase accumulated state
        psm = configure_peptide_stub_mover(insert_res)
        # st = configure_set_torsion_mover(residues_to_add)
        db = declare_bond(new_bond)



        gk = configure_gen_kic_mover(KIC_res, new_bond, excluded_from_pivot, scorefxn)

        # apply PeptideStubMover, SetTorsion, DeclareBond, GenKIC (in order)
        psm.apply(p)
        # st.apply(p)
 


       #rebuild fold tree here:

        ft=FoldTree()
        ft.add_edge(1,new_bond[0],-1)
        ft.add_edge(p.size(),new_bond[1],-1)
        ft.add_edge(1,p.size(),1)
        ft.reorder(1)
        p.fold_tree(ft)

        set_omegas_to_180(p, insert_res)
        db.apply(p)
        gk.apply(p)

        if _DEBUG:
            print(gk.get_last_move_status())

        attempts += 1
        if gk.get_last_move_status() == protocols.moves.MoverStatus.MS_SUCCESS:
            print('{} attempts before success!'.format(attempts))
            success=1
        else:
            success=0
        return p, success
        #gk.reset_status()


def main(argv):
  import string
  print('##### arguments ####')
  print('input_file with PDB code and start point for loop on each line')
  print('number of structures to generate for each loop length')
  print('largest size loop to try (starts with 1 res loop)')
  print('number of residues to chew in from termini')
  print('-start_res sets starting residue from command line instead of input fiel')
  args=argv[1:]
  input_data=args[0]
  nstruct=int(args[1])
  max_loop_length=int(args[2])
#  stub1=int(argv[4])
  delta=int(args[3])   # how many residues to chew in from termin
  
  stub1=999
  if args.count('start_res'):
    pos = args.index('start_res')
    stub1=int(args[pos+1])
    del args[pos]
    del args[pos]


  data=map(string.split,open(input_data,'r').readlines())
  for line in data:
    input_pdb=line[0]
    if stub1==999: #not read in command line, so read in from file
      stub1=int(line[1])
  #  input_pdb=argv[1]
  #  stub1=int(argv[2])

    print(input_pdb,stub1)
    init(extra_options='-beta -ignore_waters false')
   
    
    base = string.split(input_pdb,'.')
    fname = string.split(base[0],'/')[-1]
##
    # if are growing from ends, need to insert at least 2 residues to have four residue loop
    for loop_length in range( max(1,2-2*delta),max_loop_length):
        loop_residues=[stub1]
        for i in range(delta): loop_residues.append(stub1+1+i)

        insert_residues=[]
        new_bond=[]
        out_fname = '%s_%s_%s_{:04}.pdb'%(fname,loop_length,delta)        
        print(out_fname)
        for i in range(loop_length):
            loop_residues.append(stub1+1+delta+i)
            insert_residues.append(stub1+1+delta+i)
        loop_residues.append(stub1+1+delta+loop_length)  # 1st residue of 2nd stub
        for i in range(delta): loop_residues.append(stub1+1+delta+loop_length+1+i) # 2nd residue of 2nd stub (intercting sidechain)
        new_bond=[stub1+max(loop_length,1)+delta-1,stub1+max(loop_length,1)+1+delta-1]
        print (' loop length: %s \n'%loop_length)
       # print 'loop_length ',loop_length
        print ('loop: %s \n'%loop_residues)
        print ('insert: %s \n'%insert_residues)
        print ('new_bond: %s \n'%new_bond)
#        tname=('test_%s.pdb'%loop_length)
#        print (tname)


# first residue of insert is (stub1 + delta +1).  last residue is (stub1+delta+loop_length).  midpoint is stub1+delta + 1+(loop_length-1)/2
# also, the number of loop residues is len(loop_residues), so the middle residue is at i=len(loop_residues)/2
        non_pivot = list()
        for i, resNo in enumerate(loop_residues):
          if i == 0 or i >= len(loop_residues) - 1 or i == len(loop_residues)/2:   # change later to middle of loop
            print ('pivot: %s \n'%resNo)
            continue
          non_pivot.append(resNo)

        if len(non_pivot) == 0:
        # is this a three-residue loop? is that a problem?
            sys.exit('No residues excluded from pivot')

   
   
        orig_pose = pose_from_file(input_pdb)
        sf = get_score_function()
        for i in range(nstruct):
            pose, success = configure_and_apply_movers(orig_pose, sf, insert_residues,loop_residues,new_bond,
                                          non_pivot)
            if success: pose.dump_pdb(out_fname.format(i + 1))


if __name__ == '__main__':
    main(sys.argv)

