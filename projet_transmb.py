import argparse 
import Bio 
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from collections import defaultdict
import mpl_toolkits.mplot3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from numpy import pi, cos, sin, arccos, arange
import os
import pandas as pd

#import mayavi.mlab as m
#from ipywidgets import interact, interactive, fixed, interact_manual
#import ipywidgets as widgets

def extracts_calpha(model, res_dico, hphobic_dico) : 
    """Extracts carbone alpha coordinates of chain A protein from the PDB structure, using Biopython module. 
    Parameters
    ----------
    model : 
        From DSSP 
    hphobic_dico : dict 
        Contains the hydrophobic residu information (chain id, residu letter code, RSA (relative solvent accesible surface)). 
    Returns 
    -------
    ca_coord : array
        Contains x, y and z coordinates array of each 'exposed' hydrophobic residus
    """
    ca_coord = list()
    for chain in model: 
        if chain.id == "A" : # select only chain A 
            for res in chain : 
                res_tupl = chain.id,res.id[1]
                for atom in res :
                    if atom.get_name() == "CA" : # select only carbone alpha
                        ca_coord.append(atom.get_coord())
                        res_dico[res_tupl].extend(atom.get_coord())
      
    res_df = pd.DataFrame(res_dico).T
    res_df.columns = ["Xcoord","Ycoord","Zcoord"] # replace the columns name of the 'res_df' dataframe.
    
    # merging the two dataframes (hphobic_df and res_df) by coordinates of the 3D coordinates of CA exposed hydrophobic residu atoms.  
    # 'hphobic_df' only contains hydrophobic exposed residu information,  now it also contains their 3D coordinates.
    hphobic_df = pd.DataFrame.from_dict(hphobic_dico).T
    hphobic_df = hphobic_df.merge(res_df,left_index = True, right_index = True) 
    hphobic_df.columns = ['RESCODE', 'RSA', "Xcoord","Ycoord","Zcoord"] # replace the columns name of the 'hphobic_df' dataframe.
    ca_coord = hphobic_df[["Xcoord","Ycoord","Zcoord"]] # extracts 3D coordinates into an array. 

    # compute the center of mass 
    center_mass = ca_coord.mean(axis=0)
    print("CENTER OF MASS", center_mass)
     
    # substract center of mass values to coordinates, to translate to (0,0,0) for hphobic_df 
    hphobic_df["Xcoord"] = hphobic_df["Xcoord"]- center_mass[0] 
    hphobic_df["Ycoord"] = hphobic_df["Ycoord"]- center_mass[1] 
    hphobic_df["Zcoord"] = hphobic_df["Zcoord"]- center_mass[2]
    
    # substract center of mass values to coordinates, to translate to (0,0,0) for res_df 
    res_df["Xcoord"] = res_df["Xcoord"]- center_mass[0] 
    res_df["Ycoord"] = res_df["Ycoord"]- center_mass[1] 
    res_df["Zcoord"] = res_df["Zcoord"]- center_mass[2]
                    
    return np.array(ca_coord), res_df, hphobic_df


def determine_best_direction(npts) : 
    """Search the "best" direction which passes through the center of mass, among npts points (from sampling) that maximizes the hydrophocity score
    Parameter
    ---------
    npts : int 
        Defines the number of vector (direction)
    Returns 
    -------
    best_hphobic_score : float 
        The hightest hydrophobicity score that passes through the center of mass, among npts points 
    best_hphobic_inside_df : dataframe 
        Contains rows of CA hydrophobic atom that are inside the membrane, used for compute the "best_hphobic_score". 
    pa, pb, pc, : float, float, float 
        Defines the parameter of the equation plane : ax+by+cz+d = 0, that gives the "best_hphobic_score". 
    best_coord_inside : array 
        Contains the 3D coordinates of CA hydrophobic atom that are inside the membrane, and gives the "best_hphobic_score".
    best_pas : float
        Defines the parameter d in equation cartesian of plane, that gives the "best_hphobic_score".
    x,y,z : array, array, array
        Defines the coordinates x,y and z of sampling points on the sphere
    """  
    
    best_hphobic_score = 0
    best_coord_inside = list() 
    pas = 0
    
    # generating the point on the sphere (this following code is found in the Internet) 
    indices = arange(0, npts, dtype=float) + 0.5
    phi = arccos(1 - 2*indices/npts)
    theta = pi * (1 + 5**0.5) * indices
    x, y, z = cos(theta) * sin(phi), sin(theta) * sin(phi), cos(phi);
    
    
    for i in range(npts): 
        xs = np.linspace(-4,4,2)
        #ys = np.linspace(0,0.1,2)
        ys = np.linspace(-3,3,2)
        X,Y = np.meshgrid(xs,ys)
        a,b,c,d = x[i], y[i], z[i], 0
        coord_inside, pas  = determine_hphobic_inside_mb(hphobic_df, a,b,c,pas) 
        hphobic_inside_df, hphobic_score = compute_hydrophobic_score(hphobic_df, coord_inside)
        if best_hphobic_score < hphobic_score : 
            best_hphobic_score = hphobic_score 
            best_pas = pas 
            pa, pb, pc = a,b,c
            best_hphobic_inside_df =  hphobic_inside_df
            best_coord_inside = coord_inside
        print("BEST HSCORE",best_hphobic_score)
        
    return best_hphobic_score, best_hphobic_inside_df, pa, pb, pc, best_pas, np.array(best_coord_inside), x,y,z

def determine_hphobic_inside_mb(hphobic_df, a,b,c,pas) : 
    """Determine the hydrophic CA atom that are inside the membrane. 
    Parameters 
    ----------
    hphobic_df : dataframe 
    a,b,c : float, float, float 
        Defines the parameter of the equation plane : ax+by+cz+d = 0
    pas : int 
        Define the translation step (in angstrom) of the membrane (defined by two planes).  
    Returns
    -------
    coord_inside : array 
        Contains x,y and z coordinates of CA hydrophobic atom that are inside the membrane
    pas : int 
        Define the translation step (in angstrom) of the membrane (defined by two planes).  
    """
    coord_inside = list()
    
    for coord in zip(hphobic_df["Xcoord"],hphobic_df["Ycoord"],hphobic_df["Zcoord"]) : 
        dotsup = np.dot(np.array([coord[0], coord[1], coord[2],1]), np.array([a,b,c,7+pas]))
        dotinf = np.dot(np.array([coord[0], coord[1], coord[2],1]), np.array([a,b,c,-7+pas]))
        if dotsup>0 and dotinf<0 : 
            coord_inside.append(coord)
    coord_inside = np.array(coord_inside) 
    
    return coord_inside,pas 


def compute_hydrophobic_score(hphobic_df, coord_inside) :
    """Compute the hydrophobic score that is defined by the ratio of the number of hydrophic exposed residu to the number of total hydrophic exposed residu. 
    Parameters 
    ----------
    hphobic_df : dataframe
        Contains all hydrophic exposed residu information (chain id, residu code letter, RSA, 3D coordinates).
    coord_inside : array 
        Contains 3D coordinates of the hydrophic exposed residu information in the defined membrane.
    Returns 
    -------
    hphobic_inside_df : dataframe 
        Contains hydrophic exposed residu information (chain id, residu code letter, RSA, 3D coordinates) in the membrane. 
    hphobic_score : float 
        Contains the hightest hydrophicity score of hydrophic residu of the plane that passes the center of mass.
    """
    
    # transformation the array (coord_inside) into dataframe, and change column names by "Xcoord","Ycoord","Zcoord". 
    coord_inside_df = pd.DataFrame(coord_inside) 
    coord_inside_df.columns = ["Xcoord","Ycoord","Zcoord"]
    
    # extracts the residu code letter information from the dataframe (hphobic_df) by merging information of coordinates into the dataframe (hphobic_inside_df).
    hphobic_inside_df = hphobic_df.merge(coord_inside_df, on = ['Xcoord','Ycoord','Zcoord']) 
    
    # compute the hydrophobic score 
    hphobic_score = len(hphobic_inside_df.index)/len(hphobic_df.index) 
    
    print("HYD", hphobic_inside_df)
    return hphobic_inside_df, hphobic_score
    

def translation_mb(best_hphobic_score, hphobic_df, pas, coord_inside,a,b,c) : 
    """Translate the position of the membrane (defined by two planes), until there is no more hydrophobic residue exposed in the membrane. 
    Parameters 
    ----------
    best_hphobic_score : float 
        Supposed to contain the hightest hydrophobic score after searching of direction
    hphobic_df : dataframe 
        Contains all hydrophic exposed residu information (chain id, residu code letter, RSA, 3D coordinates)
    pas : float 
        Define the translation step (in angstrom) of the membrane, equivalent to parameter d in equation, that gives the "best_hphobic_score" after searching of direction   
    a,b,c : float 
        Defines the parameter of the equation plane : ax+by+cz+d = 0, that gives the "best_hphobic_score" that gives the "best_hphobic_score" after searching of direction
    Returns
    -------
    best_hphobic_score : float 
        Contains the hightest hydrophobic score after moving along the vector axis
    pa,pb,pc, : float, float, float
        Defines the parameters of the equation of plane : ax+by+cz+d = 0, that supposed to maximize the hydrophobicity score
    best_pas : int
        Defines the parameter d of the equation of plane : ax+by+cz+d = 0, value in angstrom. 

    """
    hphobic_score = best_hphobic_score
    best_coord_inside = coord_inside
    pa, pb, pc = a,b,c 
    best_pas = 0
    while len(coord_inside)!=0: 
        pas += 2
        coord_inside, pas = determine_hphobic_inside_mb(hphobic_df, a,b,c, pas)
        if len(coord_inside)!=0 : 
            hphobic_inside_df, hphobic_score = compute_hydrophobic_score(hphobic_df, coord_inside)
            if best_hphobic_score < hphobic_score : 
                best_hphobic_score = hphobic_score 
                best_pas = pas 
                best_coord_inside = coord_inside
                pa, pb, pc = a,b,c

                print("BEST + ", best_hphobic_score, a,b,c,best_pas)
       
    
    pas = 0
    coord_inside = best_coord_inside

    while len(coord_inside)!=0: 
        pas -= 2
        coord_inside, pas = determine_hphobic_inside_mb(hphobic_df, a,b,c, pas)
        if len(coord_inside)!=0 : 
            hphobic_inside_df, hphobic_score = compute_hydrophobic_score(hphobic_df, coord_inside)
            if best_hphobic_score < hphobic_score : 
                best_hphobic_score = hphobic_score
                best_pas = pas
                best_coord_inside = coord_inside
                pa, pb, pc = a,b,c

                print("BEST -", best_hphobic_score, pa,pb,pc,best_pas)
    
    return best_hphobic_score, best_coord_inside, pa,pb,pc, best_pas 

def plot_membrane(res_df, hphobic_df, hphobic_inside_df, a,b,c,pas,x,y,z, xcoord_inside, ycoord_inside, zcoord_inside) :
    """Give final vizualisation of the possible "best" position of membrane that maximizes the hydrophobicity score
    Parameters
    ----------
    res_df : dataframe 
        Contains protein chain A residu information (chain id, residu code letter, RSA, 3D coordinates)
    hphobic_df : dataframe 
        Contains all hydrophic exposed residu information (chain id, residu code letter, RSA, 3D coordinates)
    hphobic_inside_df : dataframe  
        Contains hydrophic exposed residu information (chain id, residu code letter, RSA, 3D coordinates) in the membrane. 
    a,b,c : float, float, float
        Defines the parameters of the equation of plane : ax+by+cz+d = 0, that gives the hightest hydrophobicity score. 
    xcoord_inside, ycoord_inside, zcoord_inside : array, array, array
        Defines the coordinates of the hydrophobe exposed residu coordinates, inside the membrane
    """
    d = pas
    # set figure parameters 
    fig= plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection='3d')
    ax.set(xlim=(-4,4), ylim=(-4,4), zlim=(-4,4))
    
    # plots the all carbon alpha using 3D coordinates from the dataframe (res_df) 
    ax.scatter(res_df["Xcoord"], res_df["Ycoord"], res_df["Zcoord"], s=1, color = "green") 

    
    # plots the carbon alpha of CA hydrophobic from chain A, using 3D coordinates from the dataframe (hphobic_df) 
    ax.scatter(hphobic_df["Xcoord"],hphobic_df["Ycoord"],hphobic_df["Zcoord"], s=8, color = "red") 

    # plots sampling points on sphere
    ax.scatter(x, y, z, color = "blue");
    # plots the vector 
    ax.quiver(0,0,0, a,b,c, length=5, normalize = True, color ="red")

    # generating the grid 
    xs = np.linspace(-5,5,2)
    ys = np.linspace(-4,4,2)
    X,Y = np.meshgrid(xs,ys)
    Z = (-d-a*X - b*Y) / c  # cartesian equation according to Z 
    Zsup = ((7+d - a*X - b*Y) / c) # position of upper plane
    Zinf = ((-7+d - a*X - b*Y) / c) # position of down plane 
    
    # plot the orthogonal plane that passes the center of mass 
    ax.plot_surface(X, Y, Z,color="red", alpha =0.5)
        
    # plot the membrane defined by two planes, with a distance of 14 angstrom 
    ax.plot_surface(X, Y, Zsup,color="blue", alpha = 0.3)
    ax.plot_surface(X, Y, Zinf,color="green", alpha = 0.3)
    
    # plot the carbone alpha inside the membrane 
    ax.scatter(xcoord_inside, ycoord_inside, zcoord_inside, color = "black", s = 16 );    

####################
### Main program ###
####################

parser = argparse.ArgumentParser()
#parser.add_argument('-f', "--file",  action = 'store_true', default = False, help = 'Get the PDB id')
parser.add_argument('-i', action='store', dest='path',
                    help='directory of pdb file', default = "/home/courage/ProjetsM2/2n90.pdb", type=str)
pdb_path = parser.parse_args() 
pdb_id = os.path.basename(pdb_path.path).strip(".pdb") # the PDB id name 

p = PDBParser()
structure = p.get_structure(pdb_id, pdb_path.path) 
model = structure[0]
dssp = DSSP(model, pdb_path.path) 

hydrophobic = ['F', 'G', 'I', 'L', 'M', 'V', 'W', 'Y'] # list of amino acid that are considered as hydrophobic residu 
hphobic_dico = defaultdict(list) 
res_dico = defaultdict(list)
res_df = pd.DataFrame()

# on PDB structure file, compute the RSA value using DSSP tool 
# browse the output "dssp" and extracts the chain id, residu id, residu name and rsa value for each residu into the dico 'hphobic_dico' when their rsa value is above 0.3. 
for key, value in zip(dssp.keys(), dssp) : 
    chainid, resid, resname, res_rsa = key[0][0], key[1][1], value[1], value[3]
    if (chainid == "A") and (resname in hydrophobic) and (res_rsa >=0.3) : 
        hphobic_dico[chainid,resid].extend([resname,res_rsa])
        
ca_coord, res_df, hphobic_df= extracts_calpha(model, res_dico,hphobic_dico)
print("RES DF", res_df)
print("HYDROPHOBIC DF", hphobic_df)


# initialise the number vector and the parameter 'pas'.
npts = 51
pas = 0

# determine the best direction 'vector' that maximizes the score of hydrophobicity 
hphobic_score, hphobic_inside_df, a, b, c, pas, coord_inside, x,y,z = determine_best_direction(npts)
print("Lorsqu'on est orthogonal au vecteur, le meilleure score d'hydrophobicité est de : {}, la position de la membrane se trouve à : {}*x+{}*y+{}*z+{}".format(hphobic_score, a,b,c,pas)) 
        
# determine the best slice along the vector that maximizes the score of hydrophobicity
hphobic_score, coord_inside, a,b,c,pas = translation_mb(hphobic_score, hphobic_df, pas, coord_inside, a,b,c)
print("Sur le longueur du vecteur, le meilleure score d'hydrophobicité est de : {}, la position de la membrane se trouve à : {}*x+{}*y+{}*z+{}".format(hphobic_score, a,b,c,pas)) 

# extracts x, y and z coordinates inside the membrane from the array 'coord_inside', in order to plot them. 
xcoord_inside = np.array([coord[0] for coord in coord_inside])
ycoord_inside = np.array([coord[1] for coord in coord_inside])
zcoord_inside = np.array([coord[2] for coord in coord_inside])

# plot the membrane and carbone alpha hydrophic exposed residu atoms of chain A protein.
plot_membrane(res_df, hphobic_df, hphobic_inside_df, a,b,c, pas, x ,y,z, xcoord_inside, ycoord_inside, zcoord_inside)
plt.show()
