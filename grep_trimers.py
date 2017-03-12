#!/usr/bin/python
from os import system,popen,path
from sys import argv,stderr
import string

print 'input list of match pdbs and native pdb '

input=argv[1]
nat_pdb=argv[2]
    
pdbs=open(input,'r').readlines()
for pdb in pdbs:
    if nat_pdb=='LIST':
        t=string.split(pdb[:-1],'/')[-3]
        index = t.find('_')
        native=t[0:index] + '.pdb'
    else:
        native=nat_pdb
    tag=string.split(pdb[:-1],'/')[-1]
 
    index= tag.find('.pdb')
    pdb_tag= tag[0:index]
    lines1=popen('grep " X " %s'%pdb).readlines()
   # print lines1
    nres=string.split(popen('grep " X " %s | grep CA | wc'%pdb[0:-1]).readline())[0]
    lines2=[]
    for n in range(99,102):
#        print 'grep  "Y %3d" %s'%(n,pdb)
        lines2.append(popen('grep  "Y %3d" %s'%(n,pdb)).readlines())
    print nres,pdb_tag,native
    outfile=open('t.pdb','w')
#    for line in lines1:
#        print line
#        outfile.write(line)
    for res in lines2:
        for line in res:
    #        print line
            outfile.write(line)
    outfile.close()
    command = 'cat %s t.pdb > %s_%s'%(native,pdb_tag,native) 
    print command
    system(command)
