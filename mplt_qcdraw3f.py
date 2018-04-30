import matplotlib.pyplot as plt
from mo_draw_ini import *

#### print option
tag_print=1

### line draw function
TEMP_ENE=-1000
print_draw=0


#### main module cannot give a value to sub-module
#### another module is required to get Y values

if 'YMIN' not in globals():
    print "YMIN is not defined in globals"
    XMIN    = 0
    XMAX    = 3
    YMAX    =  0.5
    YMIN    = -1.0

X_SHIFT = [0,1,2]
AB_GAP  = 0.4

Linewidth=2.0
Linewidth1=1.0

fig=plt.figure()
ax=plt.subplot(111)
ax.set_xlim([XMIN,XMAX])
ax.set_ylim([YMIN,YMAX])

def mplt_line(degeneracy, x, y, ehomo):
    """
    mplt_line arges::  degeneracy, x value list, y value list, homo energy
    draw one level with a few lines depending on degeneracy
    degeneracy depends on the size of x_list
    maximum degeneracy is 3
    """
    if degeneracy == 1:
	    if float(y[0]) <= ehomo:
	        plt.plot(x, y, 'r' , lw=Linewidth)
	    else:
	        plt.plot(x, y, 'b', lw=Linewidth )
    elif degeneracy == 2:
        #print x, y
        x1=x[0:2]
        x2=x[2:4]
        if y[0] <= ehomo:
	        plt.plot(x1, y, 'r', x2, y, 'r', lw=Linewidth)
        else:
	        plt.plot(x1, y, 'b', x2, y, 'b', lw=Linewidth)
    elif degeneracy == 3:
	    x1=x[0:2]
	    x2=x[2:4]
	    x3=x[4:6]
	    if y[0] <= ehomo:
	        plt.plot(x1, y, 'r', x2, y, 'r', x3, y, 'r', lw=Linewidth)
   	    else:
	        plt.plot(x1, y, 'b', x2, y, 'b', x3, y, 'b', lw=Linewidth)
    else:
        print "Error in n_degeneracy in function mplt_line"
        exit(55)
    return 0

def mplt_line_link(x, y):
    plt.plot(x, y, 'g', ls='dashed', lw=Linewidth1)

    return 0
### 
#key_word="Alpha MOs"

def Xrange_nf(nfile,ifile,beta,ab, degen):
    """
    Xrange_nf args:: n_files, index_file, only_alpha(0)_or_beta_exists(1), alpha(0)_or_beta(1), degeneracy
    xrange: 0 ~ to XMAX
    divide xrange with n_files
        dx: xrange in each file is divided by degeneracy
    ifile indicates each section by xmin, xmax
    if beta, xmin~xmax is divided by half
    alpha locates 1st half and beta locates 2nd half
        dx: half of each file is divided by degeneracy
    """
    #### xrange for files
    xrange=float(XMAX-XMIN)/nfile
    #### devide xrange into partsf for degeneracy
    npart=2*degen+1
    dx=float(xrange)/float(npart)   # devide Xrange into parts

    #### xmin, xmax for a file
    if nfile==1:
        xmin=XMIN
        xmax=XMAX
    elif nfile==2:
        if ifile==0:
            xmin=XMIN
            xmax=xmin+float(xrange)
        elif ifile==1:
            xmin=XMIN+float(xrange)
            xmax=XMAX
    elif nfile==3:
        if ifile==0:
            xmin=XMIN
            xmax=xmin+float(xrange)
        elif ifile==1:
            xmin=XMIN+float(xrange)
            xmax=xmin+float(xrange)
        elif ifile==2:
            xmin=XMIN+float(xrange)*2
            xmax=xmin+float(xrange)
    #### xmin, xmax for alpha or beta in the file            
    if beta:
        x_2=float(xmax-xmin)/2
        if ab==0:
            #xmin=xmin
            xmax=xmin+x_2
        elif ab==1:
            xmin+=x_2
            #xmax=xmax
        else:
            print "error in function draw"
            exit(97)
        dx=float(x_2)/float(npart)

    x1=xmin+dx
    x2=x1+dx
    x3=x2+dx
    x4=x3+dx
    x5=x4+dx
    x6=x5+dx
    x7=x6+dx
    x8=x7+dx
    if degen == 1:
        return x1, x2
    elif degen == 2:
        return x1, x2, x3, x4
    elif degen==3:
        return x1, x2, x3, x4, x5, x6
    elif degen==4:
        return x1, x2, x3, x4, x5, x6, x7, x8
        
    return 0

def mplt_ab_draw(nfile, fid, beta_tag, ab_tag, diction, e_homo):
    """
    type ab_draw(int, int, int, int, hash(int, float), float)
    args:: n_files, index_file, beta_exists?, alpha_or_beta, imo_ene_dictionary, homo_energy
        sub function: Xrange_nf, mplt_line
        degeneracy is calculated using dictionary
    """
    #### change into sorted dictionary: but list of (key, value) tuples
    if tag_print >= 2:
        print "func: ab_draw ", e_homo
    l_dic=sorted(diction.items())     
    nkeys=len(l_dic)

    #### initialize variables for each file drawing
    tag_draw_beta="NO"
    #### MO level ragne is classified by nfiles
    i=0
    degen_eline=1
    #e_homo=float(FL_homo_ene[fid][0])
    #### for H radical
    if nkeys == 1:
        y=[l_dic[0][1], l_dic[0][1]]
        mplt_line(1, x_3kinds[0], y, e_homo)
        fid+=1
        #continue
    temp_ene=TEMP_ENE
    tag_last_draw="NO"

    #if re.search("\*", v_ene): # in qchem outfile, 1st low number might be ******
    #### When use dictionary, be needed to sort keys 
    for tup in l_dic:
        i+=1
        v_ene=tup[1]
        if print_draw:
            print v_ene
        if temp_ene==TEMP_ENE:         # for initial energy value, just save 
            temp_ene=v_ene
            continue
        #### draw the previous energy y
        y=[temp_ene, temp_ene]
        #### if energies are same, increase multi degeneracy 
        if temp_ene == v_ene:
            degen_eline += 1
            ### if not last energy
            if i != nkeys:
                continue
            ### draw 3 plots if last ene and finish
            else:
                x_r=Xrange_nf(nfile, fid, beta_tag, ab_tag, degen_eline)
                mplt_line(degen_eline, x_r, y, e_homo)
            #print "finish drawing"
        ### Draw if energies are different 
        else:
            ### draw 3 kind of plots with degeneracy for previous energies
            if degen_eline > 4:
                print "ERROR:: so many degeneracy in drawing"
                exit(10)
            #print "degen_eline = ", degen_eline
            #print nfile, fid, beta_tag, ab_tag
            x_r=Xrange_nf(nfile, fid, beta_tag, ab_tag, degen_eline)
            #print "ab_tag: ", ab_tag,  x_r
            mplt_line(degen_eline, x_r, y, e_homo)
            # if last ene and no degenracy, plot additional energy
            if i == nkeys:
                y=[v_ene, v_ene]
                #print nfile, fid, beta_tag, ab_tag
                #print Xrange_nf(nfile, fid, beta_tag, ab_tag, 1)
                x_r=Xrange_nf(nfile, fid, beta_tag, ab_tag, 1)
                #print x_r
                mplt_line(1, x_r, y, e_homo )
        temp_ene=v_ene
        degen_eline = 1
        #### end of one energy value
    #### end of energy list
    return 0



