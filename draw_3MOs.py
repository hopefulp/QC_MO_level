#!/home/jackjack5/epd/bin/python
### !/Users/JoonhoPark/anaconda/bin/python
### made by J. Park
### 2016. 10. 25 developed for radical
### 2016. 10. 28 developed for line
### 2016. 11. 14 developed for homo-reference
### 2016. 11. 15 developed for nbo BD
### 2016. 11. 21 developed for MO Coefficients analysis
### 2016. 11. 24 developing for nbo charge
### 2016. 12. 13 developed for MO linkage using hash
### 2017. 01. 24 developing for 6-31G(d)
### 2017. 02. 06 developing for ini.py
### required: draw only a few levels    
### saved in Tmppy/draw_MOs_link_dev1.py
### to control the display of MO levels, CTRL1, CTRL2

import sys
import re
import numpy as np
from mplt_qcdraw import *
from draw_MOs_mod import *
from mo_draw_ini import *

#### print option
ptag_link=1
ptag_molevel=1
print_level=1
print_mo=1

#### analysis option
NBO_anal=1    
if re.search('a', tag_draw_select.lower()):
    lumo_cut=Nmax_lumo
    MOCoeff_anal=0
else:   #### find FILTER for selection
    MOCoeff_anal=1

#### FL stands for variable hase file list [ [f1:...],[f2:...],[f3:...] ]
#### obtain variables from module draw_MO_ini: qchem_out files, nfile, nlumos, link_id, link option
FL_qchem_outf=File_list[:]
#### constants from module draw_MO_ini
nfile=Nfiles
FL_nlumo=Fl_nlumo[:]
if nfile >= 2:
    MOlink_hl_id=MO_link_hl_id[:]
    ######## display link ? : default is yes, if not use 'L=0'
    if 'Link_tag' not in globals():
        Link_tag=1

#### Filter by MO coefficients
MOcoeff_crit=0.15           # 0.15, 0 for all the levels

###
FL_sMO_dic=[]
FL_link_dic=[]
######## Key Word for Detection
####     key word for MO
#KW_MOene_occ_alpha="Alpha MOs"
KW_MOene_occ_alpha="Occupied"
KW_MOene_vir="Virtual"
KW_MOene_beta="Beta MOs"
#key_beta1="unrestricted"
KW_MOene_end="-----"
######## KW for MO Coefficient
KW_MOcoeff="MOLECULAR ORBITAL COEFFICIENTS" # common for alpha & beta
KW_nbasis="basis functions"

FL_homo_ene=[]      # it can have beta, save homo energy [ [f1 a, f1 b],[f2 a, f2 b]... ]
FL_homo_id=[]       # ene and MO id go parallel
FL_MOene_list=[]      # [ [e1, e2... ], [e1, e2...] [...
FL_MOene_select=[]
FL_abMO_dic_list=[]
FL_MOid_select=[]
FL_beta_tag=[]      # 0 for restricted, 1 for beta spin of  unrestricted calculation 

#### obtain atoms for MO and mol types in QChem outfile
FL_MO_type=Fl_MO_type[:]
FL_MO_atoms=Fl_MO_atoms[:]

#### depending on nbasis, num of nbo BD line is different
nbo_atoms=['Ni', 'P', 'C','O']    # atom for 1st fragment, 2nd fragment
nbasis=[22, 8, 9, 9]
nbo_BD_nline_atoms=[]
for nb in nbasis:
    nb_line5=int(nb/5)+1
    nbo_BD_nline_atoms.append(nb_line5)

KW_nao_occ="NATURAL POPULATIONS"
KW_nec="Natrual Electron"           # This has information of num of e in NAO
KW_nbo="Bond orbital"
KW_nbo_BD="BD"
KW_nbo_end="NBO state"

AngSym=('1s', '2s', '2px', '2py', '2pz')

lcmo_mol=1      
lcmo_nlnh=[1,1]
lcmoA=[]
lcmoB=[]
lcmoC=[]

"""
Constants in module
readline file and extract (1) energies, (2) coefficients, (3) nbos
"""
fid=0
# read three files
for file in FL_qchem_outf:
    tag_MO="YES"    # MO ene
    flag_MO="NO"
    flag_vline="NO"
    tag_beta="NO"
    tag_nbo="YES"
    flag_nbo="NO"
    flag_bd="NO"
    tag_moc="YES"
    flag_moc="NO"
    flag_moc_beta="OFF"
    inf=open(file, 'r')
    f_homo_i_line=[]
    f_list_ene_select=[]
    f_list_ene_sel_ab=[]
    f_list_imo_select=[]
    f_list_imo_sel_ab=[]
    f_ene_dic={}
    f_ene_dic_beta={}
    f_dic_ab_list=[]
    disp_MOids=[]
    line_ene_a=[]     # for job=optimization, whenever find new keyword, empty
    line_ene_b=[]
    i=0
    i_ene_line=0
    ibd=0
    ibd_atom=0
    nb_line=[]
    nbasis=0
    temp_id_atom='Z'    # not atom id such as Ni, Ni-1, Ni-2...
    #### for MOCoefficient of BLOCK 2
    imoc=0
    imoc_block=0
    id_tmp=0    # mo-id starts from 1 in MO Coefficient Block
    mo_ene_tmp=-100
    dump_mo_id_block=[]

    s='-'       # atom-id linker
    # read one file
    while 1:
        line=inf.readline()
        i += 1
        if not line: break         # end of while-loop, one file
        #### find nbasis for MO Coefficient here
        if MOCoeff_anal==1 and nbasis==0:
            if re.search(KW_nbasis, line) and re.search("shells", line):
                field=re.split("\s+", line)
                field=[x for x in field if x]
                nbasis=float(field[5])
                #print "found num of basis: ", nbasis
        ###### 1st: MO energies                         
        ################################ BLOCK 1 ###############################################
        #### find KW_MO energy block
        # inside one file
        if flag_MO == "NO" and tag_MO == "YES":
            if re.search(KW_MOene_occ_alpha, line):
                print "################### ", FL_qchem_outf[fid], "#########################################"
                print "###################  BLOCK 1 :: QChem MO energy"
                flag_MO="YES"
                FL_beta_tag.append(0)
                eline_homo_tmp=[]
        #### have found MO
        elif flag_MO == "YES" and tag_MO == "YES":
            '''
            after found occupied orbital energy
            options except energies
            1. detect endline: change tag_MO into NO then this block is skipped
            2. virtual line
            '''
            #### if arrived the END line of MO energy, break
            if re.search(KW_MOene_end, line):
                tag_MO="Done"
            #### symmetry lines between energy lines
            elif re.search("A'", line):
                print "skip symmetry line"
                continue
            #### normally if found VIRTUAL, go to next line,
            elif re.search(KW_MOene_vir, line):
                flag_vline="YES"    # goes to the next line
            #### in case of Beta, there is blank line, how to treat BLANK
            #### have found beta
            elif not line.strip():
                flag_vline="NO"     # virtual energy ended 
            #### have found BETA, the same as above line
            elif re.search(KW_MOene_beta, line):
                tag_beta="YES"
                FL_beta_tag[fid]=1  # BETA is defined here 
                i_ene_line=0
                f_homo_i_line.append(iline_homo)
                print "spin beta in ", FL_qchem_outf[fid]
            #### save energy line
            else:
                # in case of beta, pass KW=="Occupied" again
                if re.search(KW_MOene_occ_alpha, line):
                    pass
                # save energy line before end line
                else:
                    if tag_beta == "YES":
                        line_ene_b.append(line)
                        if flag_vline != "YES":
                            iline_homo=i_ene_line
                    else:
                        line_ene_a.append(line)
                        if flag_vline != "YES":
                            iline_homo=i_ene_line
                            #print "alpha homo", iline_homo
                    #print i_ene_line, "here"
                    i_ene_line += 1
        ###### when MO energy reading is done, come to here
        ################################ BLOCK 2 :: MO Coeff ########################################
        ###### ene & index are obtained here
        if MOCoeff_anal==1 and tag_moc=="YES": # "cut" means alpha has done
            #### Look for KW_MO coeff
            if flag_moc == "NO" :
                if re.search(KW_MOcoeff, line) :    #common for alpha & beta
                    flag_moc="YES"
                    moc_block=[]
                    f_homo_id=FL_homo_id[fid][0]
                    lumo_id_cut=f_homo_id+FL_nlumo[fid]
                    if max_lumo < lumo_id_cut:
                        print "Error:: too small lumos"
                        exit(33)
                    #print "MO cut id: ", f_homo_id, FL_nlumo[fid], lumo_id_cut
                    if flag_moc_beta == "ON":
                        f_homo_id=FL_homo_id[fid][1]
                        print "###################  BLOCK 2 :: MO Coeff beta starts "
                        imoc=0
                        imoc_block=0
                        id_tmp=0            # mo-id starts from 1 in MO Coefficient Block
                        mo_ene_tmp=-100
                        dump_mo_id_block=[]
                        f_list_imo_sel_ab=[]
                        f_list_ene_sel_ab=[]
                    else:
                        print "###################  BLOCK 2 :: MO Coeff starts "
                    #print FL_qchem_outf[fid], ":: MO Levels"
                    continue
            elif flag_moc == "YES":
                #### if end of MO Coefficients
                if line =='\n':
                    if FL_beta_tag[fid]==0:
                        tag_moc="Done"
                    continue
                #### 1. new MOC block and 1st line: 1.1 save MO id
                if imoc==0:
                    #imo=re.split("[\s\(\)]+",line)     # for parenthesis
                    imo=re.split("[\s]+",line)
                    imo=[x for x in imo if x]
                    n_ene_col=len(imo)
                    #print "imoc==0; n_ele_col=", n_ene_col, line 
                #### save energy
                elif imoc==1:
                    imo_ene=re.split("\s+",line)
                    imo_ene=[x for x in imo_ene if x]
                    del imo_ene[0]
                    #print imo_ene
                ###### arrived MO Coefficient line in MOC block: read moc line 
                elif imoc <= 2+nbasis:
                    #moc_line=re.split("[\s\(\)-]+", line) # coeff's are positive now
                    moc_line=re.split("[\s]+", line)    # coeff's might be negative 
                    moc_line=[x for x in moc_line if x]
                    #print moc_line
                    #not storing but display
                    istart = len(moc_line) - n_ene_col 
                    id_atom=moc_line[1]
                    if re.match("\d",moc_line[2]):
                        id_atom+=moc_line[2]        # id_atom == Ni1, Ni2, H30, H31, etc
                    #### FILTER 1  skip core orbital based on atom species in moc_line[1]: follow QChem basis
                    if id_atom != temp_id_atom:
                        temp_id_atom=id_atom
                        ic_atom=0
                        #print moc_line
                        ncore=Atom_Core_631g[Atom_Table.index(moc_line[1])]
                        #### if there is beta electron, makes error
                    else:
                        ic_atom+=1
                    #print "atom ID & num of core atomic orbital: ", id_atom, ncore
                    if ic_atom <= ncore:
                        pass
                    else:
                        #### for each MOC line, SAVE n_ene_col==6 (if last, n_ene_col<=6) in moc_block
                        for j in range(n_ene_col):
                            mo_symbol=[]
                            moc_tmp=[]
                            #### FILTER 2 of MOC line by MOcoeff_crit
                            if abs(float(moc_line[istart+j])) > MOcoeff_crit:
                                mo_symbol.extend((moc_line[1], moc_line[2]))
                                if istart==3:
                                    #if print_level >= 3:
                                    #    print imo[j], imo_ene[j], moc_line[istart+j]
                                    moc_tmp.extend((mo_symbol,imo[j],imo_ene[j],moc_line[istart+j]))
                                elif istart==4:
                                    #if print_level >= 3:
                                    #    print moc_line[3], imo[j], moc_line[istart+j]
                                    mo_symbol.append(moc_line[3])
                                    moc_tmp.extend((mo_symbol,imo[j],imo_ene[j],moc_line[istart+j]))
                                elif istart==5:
                                    #if print_level >= 3:
                                    #    print moc_line[3], imo[j], moc_line[istart+j]
                                    mo_symbol.extend([moc_line[3],moc_line[4]])
                                    moc_tmp.extend((mo_symbol,imo[j],imo_ene[j],moc_line[istart+j]))
                                else:
                                    print "MO Coefficients format error"
                                    exit(11)
                                #print mo_symbol,
                                #### Add moc line in moc_block of 6 mo_id
                                moc_block.append(moc_tmp)
                                #print "1st block in MO Coeff"
                                #print(np.matrix(moc_block))
                    #### if read nbasis + 2 line
                    #### HAVE READ one MOC block with 6 mo_id and TREAT it and INITIALIZE
                    if imoc == 1+nbasis:
                        imoc=0
                        imoc_block+=1
                        #### treat block of MOC here
                        #### FORMAT of moc_block: [ [   P,     1,           s], mo_id, mo_ene, mo_coeff
                        ####                      [ [atom, index, angular mom], mo_id, mo_ene, mo_coeff
                        ####                                         SORTED by  mo_id
                        #print moc_block
                        sorted_moc_block=sorted(moc_block, key=lambda x:float(x[1]))
                        for item in sorted_moc_block:
                            if print_level >=2:
                                print item # [['Ni', 'd', '4'], '100', '-0.016', '0.19560']
                        
                        #### eff_coeff for effective coefficients == filtered by valuable coefficient
                        #### SORTED_MOC_BLOCK with valence shell, effective coefficients
                        for eff_coeff in sorted_moc_block:
                            #print "format:: [[atom_kind id angular_mom] mo_id energy MO_Coeff]"
                            #print eff_coeff
                            ######## FILTER 3: cut lines using FL_nlumo 
                            #print int(eff_coeff[1]) ,"<", lumo_id_cut
                            if int(eff_coeff[1]) > lumo_id_cut:
                                #### dump MO Coefficient before break
                                #print "last: ", f_list_imo_sel_ab
                                imo_select=f_imo_dump(dump_mo_id_block, FL_MO_atoms[fid], FL_MO_type[fid],id_tmp)
                                if imo_select:
                                    f_list_imo_sel_ab.append(id_tmp)
                                    f_list_ene_sel_ab.append(mo_ene_tmp)
                                    #print "last-1: ", f_list_imo_sel_ab
                                #### save imo & ene for alpha or beta
                                f_list_imo_select.append(f_list_imo_sel_ab)
                                f_list_ene_select.append(f_list_ene_sel_ab)
                                #print f_list_imo_select[0]
                                if FL_beta_tag[fid]==0:
                                    tag_moc="Done"
                                elif FL_beta_tag[fid]==1 and flag_moc_beta=="ON":
                                    tag_moc="Done"
                                #### run beta again
                                elif FL_beta_tag[fid]==1 and flag_moc_beta=="OFF":
                                    tag_moc="YES"
                                    flag_moc_beta="ON"
                                    flag_moc="NO"       # find KW for beta
                                #### FINISH MO Coeff selection
                                break
                            ######## DUMP for selected atom line in selected 6 mo-id's ############################
                            ######## FILTER 4       :: cut by atom species here
                            ####     f_imo_dump()   :: decision for dumping the imo-block
                            #### cf. eff_coeff== [['Ni', 'd', '4'], '100', '-0.016', '0.19560']
                            #### 'pass' writes the line in dump string      of str_dump
                            #### 'continue' skips the line before writing   in str_dump
                            if FL_MO_type[fid] == "ONE":
                                if eff_coeff[0][0] == FL_MO_atoms[fid][0]:
                                    pass
                                else:
                                    continue
                            #### if eff_coeff[0][0] is one of atom_list, save
                            #### 'ALL' and '1sub' works the same but different in f_imo_dump()
                            #### more detail screening difference between 'ALL' and '1sub', refer to f_imo_dump
                            elif FL_MO_type[fid] == "ALL" or FL_MO_type[fid] == "1sub":    # 'all' for the selected atom 
                                f_tag=0
                                for atom in FL_MO_atoms[fid]:
                                    if eff_coeff[0][0] == atom:
                                        f_tag=1
                                        break
                                if f_tag==1:
                                    pass
                                else:
                                    continue
                            #### for more options for atom id w. index
                            elif FL_MO_type[fid] == "1sub-a":
                                if (eff_coeff[0][0] == FL_MO_atoms[fid][0]) or (eff_coeff[0][0]==FL_MO_atoms[fid][1] and eff_coeff[0][1]=="1") or (eff_coeff[0][0]==FL_MO_atoms[fid][2]):
                                    pass
                                else:
                                    #if eff_coeff[0][0]=='Ni' and eff_coeff[0][1]=='s':
                                    continue    # inside inner for-eff_coeff
                            elif FL_MO_type[fid] == 'NONE': # for CO2
                                pass
                            else:
                                "Error:: there is no atom type"
                                exit(11)
                            #print eff_coeff
                            imo=int(eff_coeff[1])
                            mo_ene=eff_coeff[2]

                            #### if index MO is the same, save for write
                            if imo == id_tmp:
                                if float(eff_coeff[3]) < 0:
                                    dump_str="      "+s.join(eff_coeff[0])+"\t    "+eff_coeff[3]
                                else: 
                                    dump_str="      "+s.join(eff_coeff[0])+"\t     "+eff_coeff[3]
                                dump_mo_id_block.append(dump_str)
                            #### if index MO is changed, dump
                            else:
                                #### dump the previous saved MO Coefficients
                                #### FILTER 4.5 :: the imo-block will be dumped?
                                #### imo_select returns 1 if dump or 0 for skip
                                imo_select=f_imo_dump(dump_mo_id_block, FL_MO_atoms[fid], FL_MO_type[fid], id_tmp)
                                #### RECORD and INITIALIZE
                                #print  id_tmp, mo_ene_tmp, imo_sel
                                #### imo_select was decided to be printed so save
                                if imo_select:
                                    #print id_tmp
                                    f_list_imo_sel_ab.append(id_tmp)
                                    f_list_ene_sel_ab.append(mo_ene_tmp)
                                    #print f_list_imo_sel_ab
                                #### initialize mo id-1-block
                                dump_mo_id_block=[]
                                id_tmp=imo
                                mo_ene_tmp=float(mo_ene)
                                #### selection position 1
                                #f_list_imo_sel_ab.append(imo)
                                #f_list_ene_sel_ab.append(mo_ene)
                                #### this line is added for stdout
                                #if ptag_molevel >= 2:
                                dump_str="MO level  "+str(imo)+"   "+mo_ene
                                dump_mo_id_block.append(dump_str)
                                if float(eff_coeff[3]) < 0:
                                    dump_str="      "+s.join(eff_coeff[0])+"\t    "+eff_coeff[3]
                                else:
                                    dump_str="      "+s.join(eff_coeff[0])+"\t     "+eff_coeff[3]
                                dump_mo_id_block.append(dump_str)

                        #### end of MO coeff one block
                        moc_block=[]
                        continue
                    #### have read one block of nbasis lines
                #### inside MOC block
                imoc+=1
                continue
                
        ###### Block 3 :: NBO BonDing extraction
        ######################### BLOCK 3 #####################################
        if NBO_anal==1:
            if nfile==3 and fid == 1:   #### As for complex molecule
                if flag_nbo == "NO" and tag_nbo == "YES":
                    if re.search(KW_nbo, line):
                        flag_nbo = "YES"
                        print "#############  NBO Analysis  ##########"
                        continue
                #### found NBO and the next line is "---------------------------------"                
                elif flag_nbo == "YES" and tag_nbo == "YES":
                    if flag_bd=="NO":
                        #### Only Ni is checked
                        #if re.search(KW_nbo_BD, line) and re.search('Ni', line):
                        if re.search(KW_nbo_BD, line) and (re.search('Ni', line) or re.search('O', line)):
                            flag_bd="YES"
                            bd_line=re.split("[\s\(\)-]+", line)
                            bd_line=[x for x in bd_line if x]
                            print bd_line[0], bd_line[1], bd_line[2], bd_line[3], bd_line[4], bd_line[5], bd_line[6], bd_line[7] 
                            nb_line.append(nbo_BD_nline_atoms[nbo_atoms.index(bd_line[4])])
                            nb_line.append(nbo_BD_nline_atoms[nbo_atoms.index(bd_line[6])])
                            #print nb_line
                            ibd += 1
                        elif re.search(KW_nbo_end, line):
                            tag_nbo = "Done"
                        else:
                            pass
                        continue
                    ### in the NBO-BD block, write
                    else:
                        if ibd == 1:
                            #print line
                            a=line[15:24]
                            b=line[27:33]
                            c=line[34:36]
                            d=line[40:41]
                            e=line[41:50]
                            f=line[50:51]
                            g=line[56:65]
                            if len(line) > 70:
                                h=line[65:66]
                                hh=line[71:80]
                                print "\t", a, b, c, d, e, f, g, h, hh
                            else:
                                print "\t", a, b, c, d, e, f, g
                            #print line
                            ibd+=1
                            continue
                        if ibd < nb_line[ibd_atom]+1:
                            ibd+=1
                            continue
                        else:
                            if ibd_atom==0:
                                ibd=1       #### to read second atom part
                                ibd_atom+=1
                                continue
                            else:
                                ibd=0
                                ibd_atom=0
                                nb_line=[]
                                flag_bd="NO"
                                continue
        ############################ BLOCK 1.1 #################################################
        ###### Block 4 :: extract energy from MO energy lines : this is done after Block 1
        ######              Homo Energy and ID
        if tag_MO=="Done":
            print "################### BLOCK 1.1 : make list of MO_ene"
            #print FL_qchem_outf[fid], "MO obtain HOMO Energy and id"
            f_homo_i_line.append(iline_homo)
            #print f_homo_i_line
            if print_level >= 3:
                print FL_qchem_outf[fid], "homo energy line number: ", f_homo_i_line[fid]
            #### f is for 1 file, 2d means alpha, beta
            f_ene_line2d=[]
            f_ene_line2d.append(line_ene_a)
            if tag_beta == "YES":
                f_ene_line2d.append(line_ene_b)
            if print_level >= 3:
                print f_ene_line2d
            
            f_ene_tmp=[]

            list_ene=[]
            list_ene2d=[]   # 2d for [[alpha], [beta]]
            e_homo_tmp=[]
            e_homo_id_tmp=[]
            iab=0
            # split alpha beta
            for ab_line_list in f_ene_line2d:     # alpha beta [ A[line...],B[line...]]
                list_ene_tmp=[]
                j=0     # line index
                ihomo=0
                for e_line in ab_line_list:     # line index, it might be list[0] / [ list[0], list[1]]
                    #print j, e_line
                    ene=re.split("\s+", e_line)
                    one_line_ene=[]
                    for x in ene:
                        if x :                   # remove empty entry
                            list_ene.append(x)
                            one_line_ene.append(x)
                            #list_ene_tmp.append(x)
                    #### finding homo ene level
                    if j == f_homo_i_line[iab]:
                        ehomo=list_ene[len(list_ene)-1]
                        ihomo=j*8+len(one_line_ene)
                        e_homo_tmp.append(ehomo)
                        e_homo_id_tmp.append(ihomo)
                        #print j, ehomo, ihomo
                    j+=1
                #list_ene2d.append(list_ene_tmp)
                iab += 1
            max_lumo=len(list_ene)
            FL_homo_ene.append(e_homo_tmp)
            FL_homo_id.append(e_homo_id_tmp)
            FL_MOene_list.append(list_ene)
            print "HOMO energy, ID, beta tag:: ", FL_homo_ene[fid], FL_homo_id[fid], FL_beta_tag[fid]
            #if 1 in  FL_beta_tag:
            #    print "Orbital energies in ", FL_qchem_outf[fid], ":: alpha, beta are combined"
            #### first block and last 4-th block are skipped
            tag_MO="Finish" 

        # end of last if-block
    # end of read 1 file
    inf.close() 
######## FILTER 5: retailing plot list in each file
#### cut lower energy limit for drawing in figure
    #### make dictionary for each file
    #### change into all levels for tag_draw_select=='all'
    if tag_draw_select == 'a':
        j=1
        tmp_list=[]
        tmp_elist=[]
        not_digit=[]
        #print FL_homo_id[fid]
        nlumo_cut=min(max_lumo, FL_homo_id[fid][0]+Nmax_lumo)
        nhomo_cut=max(1,FL_homo_id[fid][0]-Nmax_homo)
        print "draw all level :", nhomo_cut, nlumo_cut
        for ene in list_ene:
            if re.search("\*", ene):
                pass
            elif nhomo_cut <= j and j <= nlumo_cut:
                tmp_elist.append(float(ene))
                tmp_list.append(j)
            j+=1
        
        f_list_ene_select.append(tmp_elist[:])
        f_list_imo_select.append(tmp_list[:])
        #print f_list_ene_select[0]
    if print_level >= 3:
        print "After imo selection, make dictionary: ", f_list_imo_select
        print f_list_ene_select
    print "\"Draw MO levels of ", FL_qchem_outf[fid], "\""
    print "MO_level_ID   Energy"
    #print tag_draw_select
    #print f_list_imo_select[0]
    #print f_list_ene_select[0]
    for imo, moe in zip(f_list_imo_select[0],f_list_ene_select[0]) :
        i_imo=imo
        f_mo=moe
        print "%11d %8.3f" % (i_imo, moe)
        f_ene_dic[i_imo]=f_mo
    f_dic_ab_list.append(f_ene_dic)
    if FL_beta_tag[fid] == 1:
        #### change key for alpha, beta mo-id
        for imo, moe in zip(f_list_imo_select[1],f_list_ene_select[1]) :
            i_imo=imo
            f_mo=moe
            print i_imo, moe
            f_ene_dic_beta[i_imo]=f_mo
        f_dic_ab_list.append(f_ene_dic_beta)
    #### code for beta is developed up to here
    #print f_ene_dic
    FL_MOid_select.append(f_list_imo_select[0])
    FL_MOene_select.append(f_list_ene_select[0])
    FL_sMO_dic.append(f_ene_dic)
    FL_abMO_dic_list.append(f_dic_ab_list)
    fid += 1
######################################################
###### end of reading all  the files
print "############################################"
print "############### Result #####################"
print "        Filename        HOMO_energy & ID  beta_tag "
for i in range(nfile): 
    print "%20s" % FL_qchem_outf[i], 
    print FL_homo_ene[i], FL_homo_id[i], FL_beta_tag[i]

################# Block 5 ::PLOT using MPL ###########
str="--------"
"""
MO_link_id=(94,99)
                    for left-side link      for right-side link
                 1st-mol 2nd-mol            2nd-mol 3rd-mol
MOlink_hl_id=[[[ 0       ,0  ],[0 ,   1]], [[0,       1],[    1,     1]]]
We need          y-val  y-val y-val y-val 
"""
print "##### Block5: Final MO energies for mpl plot:: "


"""
Link is only working for non-degenerate MOs
"""
MO_link_hl_ene=[]
#### Make MO-link pair
if nfile >= 2 and 'MOlink_hl_id' in locals() and Link_tag:
    iside_link=0     # 0 for left side 1 for right side
    ihomo_f1=FL_homo_id[0][0]       # due to alpha, beta
    ihomo_f2=FL_homo_id[1][0]
    if nfile ==3: ihomo_f3=FL_homo_id[2][0]
    for lr_link in MOlink_hl_id:
        link_y=[]
        for link in lr_link:
            print "Link ID:: ", 
            print link
            #### find MO id using homo-id
            #### between MO 1 and MO 2
            if iside_link==0:
                key_f1=ihomo_f1+link[0]
                key_f2=ihomo_f2+link[1]
                if key_f1 not in FL_sMO_dic[0]:
                    print "Link dictionary error:: exit(44)"
                    print "use option L=0"
                    exit(44)
                if key_f2 not in FL_sMO_dic[1]:
                    print key_f2, FL_sMO_dic[1]
                    print "Link dictionary error:: exit(45)"
                    print "use option L=0"
                    exit(45)
                li_ene1=FL_sMO_dic[0][key_f1]
                li_ene2=FL_sMO_dic[1][key_f2]
                link_ene=[li_ene1, li_ene2]
                link_y.append(link_ene)
            #### between MO 2 and MO 3
            elif iside_link==1:
                key_f2=ihomo_f2+link[0]
                key_f3=ihomo_f3+link[1]
                if key_f2 not in FL_sMO_dic[1]:
                    print "Link dictionary error:: exit(46)"
                    print "use option L=0"
                    exit(46)
                if key_f3 not in FL_sMO_dic[2]:
                    print "Link dictionary error:: exit(47)"
                    print "use option L=0"
                    exit(47)
                li_ene2=FL_sMO_dic[1][key_f2]
                li_ene3=FL_sMO_dic[2][key_f3]
                link_ene=[li_ene2, li_ene3]
                link_y.append(link_ene)
            else: 
                print "Error"
                exit(48)
        MO_link_hl_ene.append(link_y)
        iside_link+=1
        #print link_y
#### 3D list, draw with x-coord's
if ptag_link >= 2:
    print "Link energy:: "
    print MO_link_hl_ene
"""
    type ab_draw(int, int, int, int, hash(int, float), float)
        sub function: 
"""

fid=0
#### There might be alpha, beta two dictionaris for one file
for f_abMO_dic in FL_abMO_dic_list:
    #### draw one or two dictionalry
    beta_tag=FL_beta_tag[fid]
    if beta_tag==0:
        # number of files, file id, beta?, alpha or beta, dic[a or b], homo(a or b)
        mplt_ab_draw(nfile, fid, beta_tag, 0, f_abMO_dic[0], float(FL_homo_ene[fid][0]))
    else:           # 0 for alpha, 1 for beta
        mplt_ab_draw(nfile, fid, beta_tag, 0, f_abMO_dic[0], float(FL_homo_ene[fid][0]))
        mplt_ab_draw(nfile, fid, beta_tag, 1, f_abMO_dic[1], float(FL_homo_ene[fid][1]))
    fid+=1

#### draw MPL link line
if "Link_tag" in globals():
    if nfile == 2:
        i=0
        x=Xrange_nf(nfile, 0, 0, 0, 1)
        x2=x[1]
        x=Xrange_nf(nfile, 1, 0, 0, 1)
        x3=x[0]
        xl_link=[x2, x3]
        for pair in MO_link_hl_ene:
            mplt_line_link(xl_link, pair)
            i+=1
            

    elif nfile == 3:
        i=0
        x=Xrange_nf(nfile, 0, 0,0, 1)
        x2=x[1]
        x=Xrange_nf(nfile, 1, 0,0, 1)
        x3=x[0]
        x4=x[1]
        x=Xrange_nf(nfile, 2, 0,0, 1)
        x5=x[0]
        xl_link=[x2, x3]
        xr_link=[x4, x5]
        for rl_link in MO_link_hl_ene:    
            for pair in rl_link:
                if i==0: x=xl_link
                else:    x=xr_link
                mplt_line_link(x, pair)
            i+=1
plt.show()
print "normal stop"	

