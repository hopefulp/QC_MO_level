import sys
import re
import os
from atom_valence import *       # fatom_valence
#print "\n".join(sys.path)

######################## MAIN MODULE: draw_MOs_link #####################################
Nmax_lumo=10            # normally 3
Nmax_homo=30
File_list=[]

"""
scan sys.argv here, then provide variables to main module
if atom input exists, it becomes selective in MO drawing
    if natom != 1, mo_type should be input
if not atom input, default is S='all'
    if S != 'all', atom list, mo_type should be coded
"""
#### MO type for selected atoms
#### ONE    : display only one atom
#### 1sub   : if there is main atom, display all the atoms
#### ALL    : display all the atoms in the list
#### NONE   : draw all atoms w.o. list
MO_type_list=['ONE', '1sub', '1sub-a', 'ALL', 'NONE']

# Atom id such as C-1, C-2 are not coded yet
#FL_Atom_types=[             
#    [['Ni','d']],           # 1st type
#    [['C','p'],['O','p']]   # 2nd type
#    ]

#### read qchem output file upto 3
narg=len(sys.argv)
if narg <= 1 :
    print "Error with no arguments: input qchem outfile"
    print "Usage:: sys.argv[0] qchem.out (with job=sp) [L=0|1] [S=a|s]"
    print " in case nfile==3, 1st 3rd are fragments and the middle is complex"
    exit(1)

for arg in sys.argv:
    if re.search("out", arg):
    	File_list.append(arg)
    if re.search("=", arg):
        option=arg.split('=')
        if re.search("L", option[0].upper()):
            Link_tag=int(option[1])
        elif re.search("S", option[0].upper()):
            tag_draw_select=option[1]
        elif re.search("A", option[0].upper()):
            Sel_atoms=option[1]
        elif re.search("T", option[0].upper()): # Atom list 'T'ype
            mo_type=option[1]

Fl_MO_atoms=[]
Fl_MO_type=[]
for fname in File_list:
    #### just to gain atoms in filename for type analysis below else-block
    Adic_fname=fatom_valence(fname)
    #print "in", os.path.basename(__file__), "atom_dic from filename of", fname,":", Adic_fname

    #### if atoms are input, draw method is selective
    if "Sel_atoms" in locals():
        MO_atoms=re.split("\s",Sel_atoms)
        tag_draw_select='s'
        if len(MO_atoms) == 1:
            mo_type='ONE'
        else:
            #### if there are atom_list, and more than one atom, mo_type should be input
            if 'mo_type' not in locals():
                print "Error:: with more than 2 atoms in A='atoms', input atom list type"
                print "Usage:: A='atom1 atom2' T[ype]= all[1sub]"
                exit(11)
        if 'mo_type' in locals() and re.search("s", mo_type):
            mo_type='1sub'
    #### if not atom selection            
    #### for The Present Work Park
    else:
        #### if draw all MO, skip this part
        if 'tag_draw_select' not in locals():
            tag_draw_select='a'
        if re.search('a', tag_draw_select.lower()): 
            pass

        #### FOR INPUT of SEVERAL FILES w. different TAGS, Modify here
        #### if you find 'P' or 'C', 'O' in filename, change the tag
        #### P for Ni caltalyst -CO2 system and A-fragment
        elif 'P' in Adic_fname.keys() and not 'C' in Adic_fname.keys():
            if re.search("A", fname): ####
        #### P for Ni caltalyst -CO2 system and A-fragment
                if 'N' in Adic_fname.keys():
                    #### for 3-PNP
                    MO_atoms=['Ni','P','N']
                    mo_type='1sub'
                else:
                    #### for 1-PP & 2-PPP
                    MO_atoms=['Ni','P']      # draw only for Ni
                    mo_type='1sub'
            #### for Ni-ligand and CO2 complex            
            else:
                #### 1sub-a for additional atom id of C-1
                MO_atoms=['Ni', 'C', 'O']
                mo_type="1sub-a"     #### Ni is necessary and 1 of C or O
        #### for CO2
        elif 'C' in Adic_fname.keys() and 'O' in Adic_fname.keys():
            MO_atoms=['C', 'O-1', 'O-2']
            mo_type="NONE"        # print all the atoms
        else:
            "Error for Molecular type for selective MO drawing"
    ######## The END of Ni-PPP system works in the else block
    if tag_draw_select=="s":
        Fl_MO_atoms.append(MO_atoms)
        Fl_MO_type.append(mo_type)

if tag_draw_select=="s":
    print  "mo_type list: ", MO_type_list
    for type in Fl_MO_type:
        if not type in MO_type_list:
            print "type should be one of ", MO_type_list
            exit(10)
    print "atoms for MO in ", os.path.basename(__file__), Fl_MO_atoms
    print "MO_type in ", os.path.basename(__file__) ,": ", Fl_MO_type

print "input files in module ", os.path.basename(__file__),  ": ", File_list
Nfiles=len(File_list)

'''
1                                       88 homo    
2   92(w) 93(w) 94(s,bo,O1) 99(s,bb,C-O1)           100(s,bb-ab,C-O1
3   CO2                                             lumo
Link between MOs are given manually
indices: homo-1=-1 homo=0, lumo=1, lumo-1=2
MO_link_hl_id=[[id_1st_mol, id_2nd_mol],[id_2nd_mol, id_3rd_mol]]
    1st link between 1-2 mols 2nd link between 2-3 mols
'''

#### Fl_nlumo: file_list for maximum num of lumo for each file
#### MO_link_hl_id: pair between 1-2 and 2-3 fragments
if Nfiles == 1:
    XMAX=1
    Fl_nlumo=[Nmax_lumo]
#### make link if Nfiles >= 2    
elif Nfiles == 2:
    XMAX=2
    Fl_nlumo=[Nmax_lumo,Nmax_lumo]
    MO_link_hl_id=[[[0,0],[1,1]]]    # does it need for 2 files: CO2 vs bent CO2
elif Nfiles == 3:
    XMAX=3
    if re.search("1-P", File_list[0]): 
        #MO_link_hl_id=[[[0,1],[0,0],[-1,-1],[-2,-2],[-3,-3],[-4,-4],[-2,-6]], [[1,1],[0,1],[-1,0],[-2,-1],[-3,-1],[-4,0],[-6,-1],[-32,0]]]    # 1-PP w. CO2-bended
        MO_link_hl_id=[[[0,2],[0,0],[-1,-1],[-2,-2],[-3,-3],[-4,-4]], [[2,1],[0,1],[-1,0],[-2,-1],[-3,-1],[-4,0]]]    # 1-PP G631s
        #MO_link_hl_id=[[[0,2],[0,0]], [[2,1],[0,1]]]    # 1-PP G631s
        Fl_nlumo=[1,2,2]        # cut number of lumo but include degeneracy of lumo of CO2
    elif re.search("2-P", File_list[0]): 
        MO_link_hl_id=[[[0,0],[0,5]],[[0,1],[5,1]]]
        Fl_nlumo=[2,5,2]        # cut number of lumo but include degeneracy of lumo of CO2
    elif re.search("3-P", File_list[0]):
        MO_link_hl_id=[[[0,0],[0,6]],[[0,1],[6,1]]]
        Fl_nlumo=[2,6,2]        # cut number of lumo but include degeneracy of lumo of CO2
    else:
        Fl_nlumo=[Nmax_lumo,Nmax_lumo,Nmax_lumo]
elif Nfiles == 4:
    XMAX=4
    print "error"
    exit(10)
elif Nfiles == 5:
    XMAX=5
    if re.search("1-P", File_list[0]): 
        MO_link_hl_id=[[[0,0],[1,1]],[[0,2],[0,0],[-1,-1],[-2,-2],[-3,-3],[-4,-4]], [[2,1],[0,1],[-1,0],[-2,-1],[-3,-1],[-4,0]],[[-1,0],[0,0],[1,1],[2,1]]]    # 1-PP G631s
        Fl_nlumo=[1,1,2,2,2]        # cut number of lumo but include degeneracy of lumo of CO2
    elif re.search("2-P", File_list[0]): 
        MO_link_hl_id=[[[0,0],[1,1]],[[0,0],[0,5]],[[0,1],[5,1]],[[-1,0],[0,0],[1,1],[2,1]]]
        Fl_nlumo=[2,2,5,2,2]        # cut number of lumo but include degeneracy of lumo of CO2
    elif re.search("3-P", File_list[0]):
        MO_link_hl_id=[[[0,0],[1,1]],[[0,0],[0,6]],[[0,1],[6,1]],[[-1,0],[0,0],[1,1],[2,1]]]
        Fl_nlumo=[2,2,6,2,2]        # cut number of lumo but include degeneracy of lumo of CO2
    else:
        Fl_nlumo=[Nmax_lumo,Nmax_lumo,Nmax_lumo,Nmax_lumo,Nmax_lumo]
else:
    print "Error: too many files of ", Nfiles, " : ", os.path.basename(__file__)
    exit(11)
#print FL_cut_nlumo
#exit(99)






########################## MODULE: mplot_qcdraw ######################################
#### x divides 3 parts for 3 files
XMIN=0
if 'XMAX' not in locals():
    print "Error in ", os.path.basename(__file__), ": no XMAX defined"
    exit(15)

# for Ni-CO2 0.05, 1.0 for CO2 -charges
'''
if re.search("1-P", File_list[0]):
    YMAX    =  0.05
    YMIN    = -0.4
elif re.search("2-P", File_list[0]):
    YMAX    =  0.05
    YMIN    = -0.3
elif re.search("3-P", File_list[0]):   
    YMAX    =  0.2
    YMIN    = -0.3
else:
    YMAX    = 0.5       # 0.2
    YMIN    = -0.6      # -0.7
'''

YMAX    =  0.05
YMIN    = -0.4
