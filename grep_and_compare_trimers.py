#!/usr/bin/python
from os import system,popen,path
from sys import argv,stderr
import string

print ' this version extracts central three coordinates for each peptide of each type bound to a particular repeat protein \n identifies those with the same motif in multiple places \n and prints out those pairs where the N and C termini of central tripeptide are close enough to be joined by a loop \n'
print 'input pdb list'

input=argv[1]
    
pdbs=open(input,'r').readlines()
coords={}
names={}
for pdb in pdbs:

    t=string.split(pdb[:-1],'/')[-3]
    index = t.find('_')
    design=t[0:index] 

    t=string.split(pdb[:-1],'/')[-1]
    motif=string.split(t,'_')[1]

    lines2=[]
    for n in range(99,102):
#        print 'grep  "Y %3d" %s'%(n,pdb)
        lines2.append(popen('grep  "Y %3d" %s'%(n,pdb)).readlines())

    if (design,motif) in coords.keys():
        coords[(design,motif)].append(lines2)
        names[(design,motif)].append(pdb)
    else:
        coords[(design,motif)]=[]
        names[(design,motif)]=[]
        coords[(design,motif)].append(lines2)
        names[(design,motif)].append(pdb)

#print coords

cutoff=64.0  # distance**2 cutoff between termini 
for key in coords.keys():
    ngood=0
    print key, len(coords[key])
    if len(coords[key])>1:
        ## have more than one hit with this motif on this repeat protein
        coord=coords[(key)]
        nhits=len(coord)
#        print key

        for hit1 in range(nhits):
            min_dist=cutoff
            best_hit=999
            N = string.split(coord[hit1][0][0])[6:9]
            mid1=string.split(coord[hit1][1][1])[6:9]
            #print N
            for hit2 in range(nhits):
                 if hit2 == hit1: continue
 #                print (coord[hit2][-1][-1])[:-1]
                 C=string.split(coord[hit2][-1][-1])[6:9]
                 mid2=string.split(coord[hit2][1][1])[6:9]
                 dist=0.
 #                print N, C
                 dist_mid=0.
                 for i in range(3):
                      dist=dist+( float(N[i])- float(C[i]))**2
                      dist_mid=dist_mid+ ( float(mid1[i])-float(mid2[i]))**2
                 if dist < cutoff and dist_mid > 20. : 
                     ngood=ngood+1
                     print 'good',key, ngood, hit1,hit2,dist
                     name='%s_%s_%s.pdb'%(key[0],key[1],ngood)
                     outfile=open('temp.pdb','w')
                     for residue in coords[key][hit1]:
                         for atom in residue:
                       	    outfile.write(atom)
                     for residue in coords[key][hit2]:
                         for atom in residue:
			   outfile.write(atom) 
                     outfile.close()
                     command = 'cat temp.pdb /work/baker/repeat_peptide/designs/%s.pdb > temp2.pdb'%(key[0])
                     print command
                     system(command)
                     system('cat temp2.pdb | ~krypton/scripts/renum.pl  > %s'%name)
                 if dist < min_dist and dist_mid > 20.:

                        min_dist=dist
                        best_hit=hit2
                        best_key=key
            if min_dist< cutoff:
                 print key,hit1,best_hit,min_dist,'best'
#             for i in range




#     print nres,pdb_tag,native
#     outfile=open('t.pdb','w')
# #    for line in lines1:
# #        print line
# #        outfile.write(line)
#     for res in lines2:
#         for line in res:
#     #        print line
#             outfile.write(line)
#     outfile.close()
#     command = 'cat %s t.pdb > %s_%s'%(native,pdb_tag,native) 
#     print command
#     system(command)





