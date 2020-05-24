### copyright by Mingchen Chen
### Updated on 7Dec/2019
### New features: add tertiary_frustration.dat file in AWSEM format;  Should be able to get rid of contact.map file;
### Add pymol visulaization part;
### Remove unnecessary sections for clarifications;

import random
import scipy.io
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm as cm
from scipy.stats import norm
import matplotlib.mlab as mlab
import os
from numpy import array
import sys
import warnings
from math import *
warnings.filterwarnings("ignore")


def vector(p1, p2):
    return [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]

def vabs(a):
    return sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2))

def vector_center(p1,p2):
    return [0.5*p2[0]+0.5*p1[0], 0.5*p2[1]+0.5*p1[1], 0.5*p2[2]+0.5*p1[2]]

### get the representative atom in each residue to record distances
def get_Atom_Ligand(residue):
    summ = [];
    temp_sum_update = 99999999;
    for atom1 in residue:
        temp_sum = 0.0
        for atom2 in residue:
            diff_vector  = atom1.coord - atom2.coord;
            temp_sum = np.sum(diff_vector * diff_vector) + temp_sum;
        summ.append(temp_sum);
        if temp_sum < temp_sum_update:
            temp_sum_update = temp_sum;
            atom_update = atom1;
    return atom_update

### get pairwise interactions from native structures
def calc_residue_dist(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    #print residue_one.get_resname(), residue_two.get_resname()
    se_map = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "HEM", "MSE"]
    atom_map = ['CB', 'CB','CB','CB','CB','CB','CB','CA','CB','CB','CB','CB','CB','CB','CB','CB','CB','CB','CB','CB', 'FE', 'CB'];
    atom1 = residue_one[atom_map[se_map.index(residue_one.get_resname())]];
    atom2 = residue_two[atom_map[se_map.index(residue_two.get_resname())]];
    diff_vector  = atom1.coord - atom2.coord
    return np.sqrt(np.sum(diff_vector * diff_vector))

### get pairwise interactions from native structures
def calc_residue_dist_new(residue_one, residue_two) :
    """Returns the C-alpha distance between two residues"""
    dist = 999999;
    for atom1 in residue_one:
        for atom2 in residue_two:
            diff_vector  = atom1.coord - atom2.coord;
            temp = np.sqrt(np.sum(diff_vector * diff_vector))
            if temp < dist:
                dist = temp;
    return temp

def calc_dist_matrix(ca_atoms) :
    """Returns a matrix of C-alpha distances between two chains"""
    reslen = len(ca_atoms);
    answer = [];
    answer_new = np.zeros((reslen, reslen), np.float)

    import Bio.PDB
    import numpy
    from numpy import array
    from Bio.PDB.PDBParser import PDBParser
    se_map = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE"]
    atom_map = ['CB', 'CB','CB','CB','CB','CB','CB','CA','CB','CB','CB','CB','CB','CB','CB','CB','CB','CB','CB','CB', 'CB'];
    p = PDBParser(PERMISSIVE=1)
    pdbcode = '3gso';
    s = p.get_structure(pdbcode, pdbcode+'.pdb')
    #chains = s[0].get_list()
    chains = s[0].get_list()
    for chain1 in chains:
        for res1 in chain1:
            for chain2 in chains:
                for res2 in chain2:
                    if res1.has_id('CB')==1 and res2.has_id('CB')==1 and (res1.get_resname() in se_map) and (res2.get_resname() in se_map):
                        answer.append(vabs(vector(res1['CB'].get_coord(), res2['CB'].get_coord())));
                    if res1.has_id('CB')==1 and res2.has_id('CB')==0 and res2.has_id('CA')==1 and (res1.get_resname() in se_map) and (res2.get_resname() in se_map):
                        answer.append(vabs(vector(res1['CB'].get_coord(), res2['CA'].get_coord())));
                    if res1.has_id('CB')==0 and res1.has_id('CA')==1 and res2.has_id('CB')==1 and (res1.get_resname() in se_map) and (res2.get_resname() in se_map):
                        answer.append(vabs(vector(res1['CA'].get_coord(), res2['CB'].get_coord())));
                    if res1.has_id('CB')==0 and res1.has_id('CA')==1 and res2.has_id('CB')==0 and res2.has_id('CA')==1 and (res1.get_resname() in se_map) and (res2.get_resname() in se_map):
                        answer.append(vabs(vector(res1['CA'].get_coord(), res2['CA'].get_coord())));
                    if (res1.get_resname() not in se_map) or (res2.get_resname() not in se_map):
                        result_temp = calc_residue_dist_new(res1, res2);
                        answer.append(result_temp)
    answer_new = array(answer).reshape(reslen, reslen);
    return answer_new

def get_index(pdbcode):
    import Bio.PDB
    import numpy
    from Bio.PDB.PDBParser import PDBParser
    se_map = ["ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE", "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL", "MSE"]
    se_map_b = ["A", "R", "N", "D", "C", "Q", "E", "G", "H", "I", "L", "K", "M", "F", "P", "S", "T", "W", "Y", "V", "M"]
    p = PDBParser(PERMISSIVE=1)
    s = p.get_structure(pdbcode, pdbcode+'.pdb')
    chains = s[0].get_list()
    ca_atoms = [];
    cid_list = [];
    atom_list = [];
    lig_pos = [];
    res_list_name = [];
    ichain = 0;

    #residue_one[atom_map[se_map.index(residue_one.get_resname())]];
    for chain in chains:
        for res in chain:
            cid_list.append(chain.id+str(res.get_id()[1]));
            print res.get_id()[1]
            if (res.get_resname() in se_map) and (res.get_resname() == 'GLY' or res.has_id('CB')==0 ) and res.has_id('CA')==1:
                ca_atoms.append(res['CA'].get_coord())
                atom_list.append('CA');
                lig_pos.append(0)
                res_list_name.append(se_map_b[se_map.index(res.get_resname())])
            if (res.get_resname() in se_map) and res.has_id('CB'):
                ca_atoms.append(res['CB'].get_coord())
                atom_list.append('CA');
                lig_pos.append(0)
                res_list_name.append(se_map_b[se_map.index(res.get_resname())])
            if (res.get_resname() in se_map) and res.has_id('CB')==0 and res.has_id('CA')==0:
                ca_atoms.append(res['O'].get_coord())
                atom_list.append('O');
                lig_pos.append(0)
                res_list_name.append(se_map_b[se_map.index(res.get_resname())])
            if res.get_resname() not in se_map:
                atom1 = get_Atom_Ligand(res)
                ca_atoms.append(res[atom1.get_name()].get_coord());
                atom_list.append(atom1.get_name());
                lig_pos.append(1);
                res_list_name.append("X")
    dist_matrix = calc_dist_matrix(ca_atoms);
    return dist_matrix, cid_list, atom_list, lig_pos, res_list_name, ca_atoms



def read_log(fname, cid_list, scheme):
    fin = open(fname, 'r');
    reslen = len(cid_list);
    mat = np.zeros((reslen,reslen));
    ene_res = np.zeros((reslen, ));

    for line in fin:
        strs = line.split()
        if strs[1]!='Res1' and strs[1]!='nonzero' and float(strs[4])<= 5.000:# and float(strs[16]) <:
            #print strs[1][3:], strs[2]
            ### typically, all the residues without detected vander-Walls repulsion will not be counted.
            ind1 = cid_list.index(strs[1][2:]);
            ind2 = cid_list.index(strs[2][2:]);
            if scheme == "Function1": ### leave out the rep term only
                ene = float(strs[-1])- 1.0*(float(strs[4]))- float(strs[10]) - float(strs[15]);
            if scheme == "Function2": ### leave both rep and atr terms.
                ene = float(strs[-1])- 1.0*(float(strs[4]))- 1.0*(float(strs[3])) - float(strs[10]) - float(strs[15]);
            if scheme == "Packing":
                ene = float(strs[-1]) - float(strs[10]) - float(strs[15]);
            ene_res[ind1] += 0.5*ene;
            ene_res[ind2] += 0.5*ene;
    for i in range(reslen):
        for j in range(reslen):
            mat[i][j]=ene_res[i] + ene_res[j];
            mat[j][i]=ene_res[i] + ene_res[j];
    return mat

def read_nat_log(fname, cid_list, scheme):
    fin = open(fname, 'r');
    reslen = len(cid_list);
    mat = 999*np.ones((reslen,reslen));
    ene_res = np.zeros((reslen, ));
    for line in fin:
        strs = line.split()
        if strs[1]!='Res1' and strs[1]!='nonzero' and float(strs[4])<= 5.000:
            #print strs[1][3:], strs[2]
            ### typically, all the residues without detected vander-Walls repulsion will not be counted.
            ind1 = cid_list.index(strs[1][2:]);
            ind2 = cid_list.index(strs[2][2:]);
            if scheme == "Function1": ### leave out the rep term only
                ene = float(strs[-1])- 1.0*(float(strs[4]))- float(strs[10]) - float(strs[15]);
            if scheme == "Function2": ### leave both rep and atr terms.
                ene = float(strs[-1])- 1.0*(float(strs[4]))- 1.0*(float(strs[3])) - float(strs[10]) - float(strs[15]);
            if scheme == "Packing":
                ene = float(strs[-1]) - float(strs[10]) - float(strs[15]);
            ene_res[ind1] += 0.5*ene;
            ene_res[ind2] += 0.5*ene;
    for i in range(reslen):
        for j in range(reslen):
            mat[i][j]=ene_res[i] + ene_res[j];
            mat[j][i]=ene_res[i] + ene_res[j];
    return mat

def decoy_stat(cid_list, contact_map, decoy_num, sep, scheme, lig_pos):
    reslen = len(cid_list)
    mat_all = np.zeros((reslen,reslen,decoy_num));
    ### count the number of zero matrices;
    bad_seq = 0;
    good_seq = 0;
    for i in range(decoy_num):
        temp = read_log(str(i+1)+'.log',cid_list, scheme);
        if temp.sum()==0:
            bad_seq = bad_seq + 1;
        else:
            good_seq = good_seq + 1;
            mat_all[:,:,good_seq - 1] = temp;
    print "The number of Bad Sequences is: " + str(bad_seq);
    print "The number of Good Sequences is: " + str(good_seq);

    pro_mat = [];
    lig_mat = [];
    res_mean = [0 for i in range(reslen)];
    res_std = [0 for i in range(reslen)];
    lig_mean = [];
    lig_std = [];
    ### test the number of decoys for protein-only
    for i in range(reslen):
        for j in range(i, reslen):
            if (lig_pos[i]==0 and lig_pos[j]==0) and (abs(i-j)>sep or cid_list[i][0]!=cid_list[j][0]) and contact_map[i,j] <=10.0: #and stat_std[i,j]<=5.0 and stat_mean[i,j] !=0:
                for k in range(good_seq):
                    if  mat_all[i,j,k] != 0.0:#and mat_all[i,j,k] <= 5.0:
                        pro_mat.append(mat_all[i,j,k])
    pro_mean = np.mean(pro_mat);
    pro_std = np.std(pro_mat)

    for i in range(reslen):
        if (lig_pos[i]==0):
            res_mean[i] = pro_mean;
            res_std[i] = pro_std;

    for i in range(reslen):
        if lig_pos[i]==1:
	    lig_mat=[];
            for j in range(reslen):
                if lig_pos[j]==0 and contact_map[i,j] <=10.0:
                    for k in range(good_seq):
                        if  mat_all[i,j,k] != 0.0:#and mat_all[i,j,k] <= 5.0:
                            lig_mat.append(mat_all[i,j,k])
            lig_mean.append(np.mean(lig_mat));
            lig_std.append(np.std(lig_mat));
            res_mean[i] = np.mean(lig_mat);
            res_std[i] = np.std(lig_mat);

    return pro_mean, pro_std , lig_mean, lig_std, res_mean, res_std

def hist_nat(mat, cid_list, contact_map, sep):
    new_mat = [];
    reslen = len(cid_list)
    for i in range(reslen):
        for j in range(i, reslen):
            if (abs(i-j)>sep or cid_list[i][0]!=cid_list[j][0]) and contact_map[i,j] <=10.0 and mat[i,j]!=999: #and contact_map[i,j] !=0.0: #and stat_std[i,j]<=5.0 and stat_mean[i,j] !=0:
                    new_mat.append(mat[i,j])


def frust_map(mat_nat, stat_mean, stat_std, contact_map, minvalue, maxvalue, sep, enable, cid_list, atom_list, lig_pos, res_mean, res_std, minvalue_l, maxvalue_l, ca_atoms, res_list_name):
    reslen = len(cid_list)
    frust = np.zeros((reslen,reslen));
    frust_d = np.zeros((reslen, reslen));
    fter = open('tertiary_frustration.tcl','w');
    fdat = open('tertiary_frustration.dat','w');
    fpml = open('tertiary_frustration.pml','w');
    tcl_index = 0;
    for i in range(reslen):
        for j in range(i, reslen):
            if (lig_pos[i]==0 and lig_pos[j]==0) and (abs(i-j)>sep or cid_list[i][0]!=cid_list[j][0]) and contact_map[i,j] <=10.0 and mat_nat[i,j]!=999:
                frust[i,j] = (mat_nat[i,j] - res_mean[i])/(res_std[i])
                frust[j,i] = (mat_nat[i,j] - res_mean[i])/(res_std[i])
                fdat.write(str(i)  + ' '+ str(j) + ' '+ cid_list[i][0] + ' '+ cid_list[j][0] + ' '+ str(ca_atoms[i][0]) + ' '+ str(ca_atoms[i][1]) + ' '+ str(ca_atoms[i][2]) + ' '+ str(ca_atoms[j][0]) + ' '+ str(ca_atoms[j][1]) + ' '+ str(ca_atoms[j][2]) + ' '+ str(contact_map[i,j]) + ' ' + res_list_name[i] + ' ' + res_list_name[j] + ' '+ str(mat_nat[i,j]) + ' '+ str(res_mean[i]) + ' '+ str(res_std[i]) + '\n')
            if (lig_pos[i]==1 and lig_pos[j]==0) and contact_map[i,j] <=10.0:
                frust[i,j] = (mat_nat[i,j] - res_mean[i])/(res_std[i])
                frust[j,i] = (mat_nat[i,j] - res_mean[i])/(res_std[i])
                fdat.write(str(i)  + ' '+ str(j) + ' '+ cid_list[i][0] + ' '+ cid_list[j][0] + ' '+ str(ca_atoms[i][0]) + ' '+ str(ca_atoms[i][1]) + ' '+ str(ca_atoms[i][2]) + ' '+ str(ca_atoms[j][0]) + ' '+ str(ca_atoms[j][1]) + ' '+ str(ca_atoms[j][2]) + ' '+ str(contact_map[i,j]) + ' ' + res_list_name[i] + ' ' + res_list_name[j] + ' '+ str(mat_nat[i,j]) + ' '+ str(res_mean[i]) + ' '+ str(res_std[i]) + '\n')
            if (lig_pos[i]==0 and lig_pos[j]==1) and contact_map[i,j] <=10.0:
                frust[i,j] = (mat_nat[i,j] - res_mean[j])/(res_std[j])
                frust[j,i] = (mat_nat[i,j] - res_mean[j])/(res_std[j])
                fdat.write(str(i)  + ' '+ str(j) + ' '+ cid_list[i][0] + ' '+ cid_list[j][0] + ' '+ str(ca_atoms[i][0]) + ' '+ str(ca_atoms[i][1]) + ' '+ str(ca_atoms[i][2]) + ' '+ str(ca_atoms[j][0]) + ' '+ str(ca_atoms[j][1]) + ' '+ str(ca_atoms[j][2]) + ' '+ str(contact_map[i,j]) + ' ' + res_list_name[i] + ' ' + res_list_name[j] + ' '+ str(mat_nat[i,j]) + ' '+ str(res_mean[j]) + ' '+ str(res_std[j]) + '\n')
    fdat.close()

    ### write pml script header files
    fpml.write("hide all\n")
    fpml.write("unset dynamic_measures\n")
    fpml.write("show cartoon, all\n")
    fpml.write("color grey, all\n")
    fpml.write("run draw_links.py\n")

    for i in range(reslen):
        for j in range(i, reslen):
            ### residue-residue interactions
            if (abs(i-j)>sep or cid_list[i][0]!=cid_list[j][0]) and lig_pos[i]==0 and lig_pos[j]==0 and contact_map[i,j] <=10.0 and mat_nat[i,j]!=999:
                i_resid = cid_list[i][1:];
                j_resid = cid_list[j][1:];
                i_chainid = cid_list[i][0];
                j_chainid = cid_list[j][0];
                atom_i = atom_list[i];
                atom_j = atom_list[j]
                if frust[i,j] <= minvalue:
                    frust_d[i,j] = -1; frust_d[j,i] = -1;
                    fter.write("set sel" + str(i+1) + " [atomselect top \"resid " + i_resid + " and chain " + i_chainid + " and name " + atom_i + "\"]\n");
                    fter.write("set sel" + str(j+1) + " [atomselect top \"resid " + j_resid + " and chain " + j_chainid + " and name " + atom_j + "\"]\n");
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos1\n");
                    tcl_index = tcl_index + 1;
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos2\n");
                    tcl_index = tcl_index + 1;
                    fter.write("draw color green\n");
                    fter.write("draw line $pos1 $pos2 style solid width 2\n")
                    fpml.write("draw_links resi " + i_resid + " and name " + atom_i + " and Chain " + i_chainid + ", resi " + j_resid +  " and name " + atom_j + " and Chain " + j_chainid + ", color=green, color2=green, radius=0.05, object_name=" + i_resid+":"+j_resid + "_green\n" );

                if frust[i,j] >= maxvalue:
                    frust_d[i,j] = 1; frust_d[j,i] = 1;
                    fter.write("set sel" + str(i+1) + " [atomselect top \"resid " + i_resid + " and chain " + i_chainid + " and name " + atom_i + "\"]\n");
                    fter.write("set sel" + str(j+1) + " [atomselect top \"resid " + j_resid + " and chain " + j_chainid + " and name " + atom_j + "\"]\n");
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos1\n");
                    tcl_index = tcl_index + 1;
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos2\n");
                    tcl_index = tcl_index + 1;
                    fter.write("draw color red\n");
                    fter.write("draw line $pos1 $pos2 style solid width 2\n")
                    fpml.write("draw_links resi " + i_resid + " and name " + atom_i + " and Chain " + i_chainid + ", resi " + j_resid +  " and name " + atom_j + " and Chain " + j_chainid + ", color=red, color2=red, radius=0.05, object_name=" + i_resid+":"+j_resid + "_red\n" );

                if enable !=0 and frust[i,j] > minvalue and frust[i,j] <= 0:
                    fter.write("set sel" + str(i+1) + " [atomselect top \"resid " + i_resid + " and chain " + i_chainid + " and name " + atom_i + "\"]\n");
                    fter.write("set sel" + str(j+1) + " [atomselect top \"resid " + j_resid + " and chain " + j_chainid + " and name " + atom_j + "\"]\n");
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos1\n");
                    tcl_index = tcl_index + 1;
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos2\n");
                    tcl_index = tcl_index + 1;
                    fter.write("draw color yellow\n");
                    fter.write("draw line $pos1 $pos2 style dashed width 2\n")
                if enable !=0 and frust[i,j] > 0 and frust[i,j] <= maxvalue:
                    fter.write("set sel" + str(i+1) + " [atomselect top \"resid " + i_resid + " and chain " + i_chainid + " and name " + atom_i + "\"]\n");
                    fter.write("set sel" + str(j+1) + " [atomselect top \"resid " + j_resid + " and chain " + j_chainid + " and name " + atom_j + "\"]\n");
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos1\n");
                    tcl_index = tcl_index + 1;
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos2\n");
                    tcl_index = tcl_index + 1;
                    fter.write("draw color magenta\n");
                    fter.write("draw line $pos1 $pos2 style dashed width 2\n")
            ### residue-ligand interactions
            if (abs(i-j)>sep or cid_list[i][0]!=cid_list[j][0]) and (lig_pos[i]==1 or lig_pos[j]==1) and contact_map[i,j] <=10.0 and mat_nat[i,j]!=999:
                i_resid = cid_list[i][1:];
                j_resid = cid_list[j][1:];
                i_chainid = cid_list[i][0];
                j_chainid = cid_list[j][0];
                atom_i = atom_list[i];
                atom_j = atom_list[j]
                if frust[i,j] <= minvalue_l:
                    frust_d[i,j] = -1; frust_d[j,i] = -1;
                    fter.write("set sel" + str(i+1) + " [atomselect top \"resid " + i_resid + " and chain " + i_chainid + " and name " + atom_i + "\"]\n");
                    fter.write("set sel" + str(j+1) + " [atomselect top \"resid " + j_resid + " and chain " + j_chainid + " and name " + atom_j + "\"]\n");
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos1\n");
                    tcl_index = tcl_index + 1;
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos2\n");
                    tcl_index = tcl_index + 1;
                    fter.write("draw color green\n");
                    fter.write("draw line $pos1 $pos2 style solid width 2\n")
                    fpml.write("draw_links resi " + i_resid + " and name " + atom_i + " and Chain " + i_chainid + ", resi " + j_resid +  " and name " + atom_j + " and Chain " + j_chainid + ", color=green, color2=green, radius=0.05, object_name=" + i_resid+":"+j_resid + "_green\n" );

                if frust[i,j] >= maxvalue_l:
                    frust_d[i,j] = 1; frust_d[j,i] = 1;
                    fter.write("set sel" + str(i+1) + " [atomselect top \"resid " + i_resid + " and chain " + i_chainid + " and name " + atom_i + "\"]\n");
                    fter.write("set sel" + str(j+1) + " [atomselect top \"resid " + j_resid + " and chain " + j_chainid + " and name " + atom_j + "\"]\n");
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos1\n");
                    tcl_index = tcl_index + 1;
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos2\n");
                    tcl_index = tcl_index + 1;
                    fter.write("draw color red\n");
                    fter.write("draw line $pos1 $pos2 style solid width 2\n")
                    fpml.write("draw_links resi " + i_resid + " and name " + atom_i + " and Chain " + i_chainid + ", resi " + j_resid +  " and name " + atom_j + " and Chain " + j_chainid + ", color=red, color2=red, radius=0.05, object_name=" + i_resid+":"+j_resid + "_red\n" );

                if enable !=0 and frust[i,j] > minvalue_l and frust[i,j] <= 0:
                    fter.write("set sel" + str(i+1) + " [atomselect top \"resid " + i_resid + " and chain " + i_chainid + " and name " + atom_i + "\"]\n");
                    fter.write("set sel" + str(j+1) + " [atomselect top \"resid " + j_resid + " and chain " + j_chainid + " and name " + atom_j + "\"]\n");
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos1\n");
                    tcl_index = tcl_index + 1;
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos2\n");
                    tcl_index = tcl_index + 1;
                    fter.write("draw color yellow\n");
                    fter.write("draw line $pos1 $pos2 style dashed width 2\n")
                if enable !=0 and frust[i,j] > 0 and frust[i,j] <= maxvalue_l:
                    fter.write("set sel" + str(i+1) + " [atomselect top \"resid " + i_resid + " and chain " + i_chainid + " and name " + atom_i + "\"]\n");
                    fter.write("set sel" + str(j+1) + " [atomselect top \"resid " + j_resid + " and chain " + j_chainid + " and name " + atom_j + "\"]\n");
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos1\n");
                    tcl_index = tcl_index + 1;
                    fter.write("lassign [atomselect" + str(tcl_index) + " get {x y z}] pos2\n");
                    tcl_index = tcl_index + 1;
                    fter.write("draw color magenta\n");
                    fter.write("draw line $pos1 $pos2 style dashed width 2\n")

    fter.write("mol modselect 0 top \"all\"\n");
    fter.write("mol modstyle 0 top newcartoon\n")
    fter.write("mol modcolor 0 top colorid 15\n")
    fter.write('color Display Background white\n')
    fter.close()

    fpml.write("zoom all\n");
    fpml.write("hide labels\n");
    fpml.close()

    #np.savetxt('contact.map', frust, delimiter='\t')
    return frust, frust_d


reslen = int(sys.argv[1]);
minvalue = float(sys.argv[2]);
maxvalue = float(sys.argv[3]);
decoy_num = int(sys.argv[4]);
sep = int(sys.argv[5]);
enable = int(sys.argv[6]); ### enable the display of neutral frustration
scheme = sys.argv[7]; #"Packing Frustratometer" or "Function Frustratometer"
minvalue_l = float(sys.argv[8]);
maxvalue_l = float(sys.argv[9]);

print "suggest input for packing frustratometer: python ~/Documents/Frust_Paper_21Nov/Frust_Post_v2_MultiChain_ManyBody_MultiLig_Gr_ContactList_7Dec.py 92 -0.5 1.0 200 9 0 Packing -0.5 1.0"
print "suggest input for Frunction1 frustratometer: python ~/Documents/Frust_Paper_21Nov/Frust_Post_v2_MultiChain_ManyBody_MultiLig_Gr_ContactList_7Dec.py 92 -2.5 0.5 200 9 0 Function1 -1.5 0.5"
print "the cutoff of minimal frustration is: " + str(minvalue);
print "the cutoff of high frustration is: " + str(maxvalue);
print "the number of decoys used is: " + str(decoy_num);
print "the sequence separation used is: " + str(sep);
if enable == 1:
    print "the neutral frustration is also displayed  "
if enable == 0:
    print "the neutral frustration will not be displayed "
print "the scheme of " + str(scheme) + " frustratometer is used";

contact_map, cid_list, atom_list, lig_pos, res_list_name, ca_atoms = get_index('3gso');
#print lig_pos
#print cid_list
#print atom_list
mat_nat = read_nat_log('./native.log', cid_list, scheme);
hist_nat(mat_nat, cid_list, contact_map, sep);
stat_mean, stat_std, lig_mean, lig_std, res_mean, res_std = decoy_stat(cid_list, contact_map, decoy_num, sep, scheme, lig_pos)
print stat_mean, stat_std, lig_mean, lig_std
frust, frust_d = frust_map(mat_nat, stat_mean, stat_std, contact_map, minvalue, maxvalue, sep, enable, cid_list, atom_list, lig_pos, res_mean, res_std, minvalue_l, maxvalue_l, ca_atoms, res_list_name)
