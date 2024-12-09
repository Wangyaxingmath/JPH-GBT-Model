# -*- coding: utf-8 -*-
"""
Created on Thu Aug  8 16:47:19 2024

@author: yxwang
"""
import numpy as np 
import pandas as pd
import gudhi as gd
from scipy.spatial.distance import pdist, squareform


P_atoms = ["C", "N", "O", "S", "H"]
L_atoms = ["C", "N", "O", "S", "P", "F", "Cl", "Br", "I", "H"]
M_ions = ["Zn", "Mg", "Mn", "Ca", "Na", "Fe", "Ni", "Cu", "Co", "K"]

M_index = {"Zn": 0, "ZN": 0, "Mg": 1, "Mg": 1, "Mn": 2, "MN": 2, "Ca": 3, "CA": 3, 
           "Na": 4, "NA": 4, "Fe": 5, "FE": 5, "Ni":6, "NI":6, "Cu": 7, "CU": 7,
           "Co": 8, "CO": 8}   

L_index = {"Cl": 6, "CL": 6, "Br":7, "Br":7}


"""
Extract coordination datas from  metalloprotein-ligand complexes
"""

def get_atom_index(atom, atom_type, atom_index):
    L = len(atom_type)
    if atom in atom_index:
        I = atom_index[atom]
        return I
    
    for i in range(L):
        if atom == atom_type[i]:
            return i   
    return -1

def get_distance_of_two_points(p1, p2):
    temp = pow(p1[0]-p2[0], 2) + pow(p1[1]-p2[1], 2) + pow(p1[2]-p2[2], 2)
    distance = pow(temp, 0.5)
    return distance
 
def extract_pocket_metal_ions_coordination(filename):
    f = open(filename)
    contents = f.readlines()
    f.close()
        
    M = {key: [] for key in M_ions}
    for line in contents:
        if line.startswith('HETATM'):
            atom = line[76:78].strip()
            if atom == "K":
                M["K"].append(line)
            else:
                index = get_atom_index(atom=atom, atom_type=M_ions, atom_index=M_index)
                if index == -1:
                    continue
                else:
                     M[M_ions[index]].append(line)
   
    point_cloud = {}
    for key, value in M.items():
        point = np.zeros((len(value), 3))
        
        c = 0
        for i in range(len(value)):
            x = float(value[i][30:38])
            y = float(value[i][38:46])
            z = float(value[i][46:54])
            point[c][0] = x
            point[c][1] = y
            point[c][2] = z
            c += 1

        point_cloud[key] = point
        
    return point_cloud
            
def extract_pocket_protein_coordination_within_cutoff_distance(filename1, filename2, cutoff):
    f = open(filename1)
    content1 = f.readlines()
    f.close()
    
    P = []
    for line in content1:
        if line[0:4] == "ATOM":
            P.append(line)
                
    g = open(filename2)
    content2 = g.readlines()
    g.close()
        
    n = len(content2)
    start = 0
    end = 0
    for j in range(n):
        if content2[j][0:13] == '@<TRIPOS>ATOM':
            start = j + 1
            continue
        if content2[j][0:13] == '@<TRIPOS>BOND': 
            end = j - 1
            break
            
    Ligand = []
    for k in range(start, end+1):
            Ligand.append(content2[k])
            
    protein_data = []     
    for protein in P:
        point1 = [float(protein[30:38]), float(protein[38:46]), float(protein[46:54])]
        for ligand in Ligand:
            point2 = [float(ligand[17:27]), float(ligand[27:36]), float(ligand[37:46])]
            distance = get_distance_of_two_points(p1=point1, p2=point2)
            if distance <= cutoff:
                protein_data.append(protein)
                break
            
    protein_point = {}
    for atom in P_atoms:
        protein_point[atom] = []
        
    for line in protein_data:
        if line[0:4] == "ATOM":
            atom = line[12:14]
            atom = atom.strip()
            index = get_atom_index(atom=atom, atom_type=P_atoms, atom_index=L_index)
            if index == -1:
                continue
            else:
                protein_point[P_atoms[index]].append(line)
            
    point_cloud = {}
    for key, value in protein_point.items():
        point = np.zeros((len(value), 3))
        
        c = 0
        for i in range(len(value)):
            x = float(value[i][30:38])
            y = float(value[i][38:46])
            z = float(value[i][46:54])
            point[c][0] = x
            point[c][1] = y
            point[c][2] = z
            c += 1

        point_cloud[key] = point
        
    return point_cloud

def extract_ligand_coordination(filename): 
    g = open(filename)
    contents = g.readlines()
    g.close()

    Ligand = {key: [] for key in L_atoms}
    
    n = len(contents)
    start = 0
    end = 0
    for j in range(n):
        if contents[j][0:13] == '@<TRIPOS>ATOM': 
            start = j + 1
            continue
        if contents[j][0:13] == '@<TRIPOS>BOND': 
            end = j - 1
            break

    for k in range(start, end+1):
        atom = contents[k][46:48]
        atom = atom.strip()
        index = get_atom_index(atom, atom_type=L_atoms, atom_index=L_index)
        if index == -1:
            continue
        else:
            Ligand[L_atoms[index]].append(contents[k])

    point_cloud = {}
    for key, value in Ligand.items():
        n = len(value)
        point = np.zeros((n, 3))
        c = 0
        for i in range(n):
            x = float(value[i][17:27])
            y = float(value[i][27:36])
            z = float(value[i][37:46])
            point[c][0] = x
            point[c][1] = y
            point[c][2] = z
            c += 1
        point_cloud[key] = point

    return point_cloud      


"""
Compute the join persistent homology 
"""
def get_face(simplex):
    all_face = []
    n = len(simplex)
    for i in range(n):
        face = []
        for j in range(n):
            if j != i:
                face.append(simplex[j])              
        all_face.append(face)
    return all_face

def construct_filtration_over_join(complex1):
    simplex2 = []
    fil2 = []
    t = len(complex1[0])
    exist = [1] * t
    while len(simplex2)<len(complex1[0]):
        for i in range(t-1,-1,-1):
            face = get_face(complex1[0][i])
            add = True
            if len(face)>1:
                for f in face:
                    if f not in simplex2:
                        add = False
                        break
            if add == True and exist[i]==1:
                simplex2.append(complex1[0][i])
                fil2.append(complex1[1][i])
                exist[i] = 0
                break
               
    return simplex2, fil2

def get_boundary_matrix(complex1, complex2):
    m = len(complex1)
    n = len(complex2)
    boundary_matrix = np.zeros((m+n, m+n))
    
    for i, simplex1 in enumerate(complex1):
        for j, simplex2 in enumerate(complex1):
            if len(simplex1) != len(simplex2):
                face = get_face(simplex2)
                if simplex1 in face:
                    boundary_matrix[i][j] = 1
                   
    for i, simplex1 in enumerate(complex1):
        for j, simplex2 in enumerate(complex2):
            if simplex1 == simplex2:
                boundary_matrix[i][j+m] = 1
                
    for i, simplex1 in enumerate(complex2):
        for j, simplex2 in enumerate(complex2):
            if len(simplex1) != len(simplex2):
                face = get_face(simplex2)
                if simplex1 in face:
                    boundary_matrix[i+m][j+m] = 1
    
    M = []
    for col_idx in range(boundary_matrix.shape[1]):
        col = boundary_matrix[:,col_idx]
        non_zero_indices = np.nonzero(col)[0]
        M.append(set(non_zero_indices))
    
    return M 

def show_matrix(m):
    temp = np.zeros((len(m),len(m)))
    for i in range(len(m)):
        for j in m[i]:
            temp[j,i] = 1
        
    return temp
   
"""
Matrix reduction
    input :
        boundary matrix: M
        M = [ col0, col1, col2, col3,... ]
        each column is a set, containing the position index with nonzero value
    
    output:
        reduced boundary matrix
"""

def get_max(m):
    if len(m)>0:
        return max(m)
    else:
        return -1           

def add_two_column(a, b):
    for item in a:
        if item in b:
            b.remove(item)
        else:
            b.add(item)
    return b

def matrix_reduction(sparse_matrix):
    t = len(sparse_matrix)
    
    L = [] 
    for i in range(t):
        L.append(-1)
    
    for j in range(t):
        N = get_max(sparse_matrix[j])
        while N!=-1 and L[N]!=-1:
            sparse_matrix[j] = add_two_column(sparse_matrix[L[N]], sparse_matrix[j])
            N = get_max(sparse_matrix[j])
        if N==-1:
            continue
        elif L[N]==-1:
            L[N] = j

    return sparse_matrix
         
"""
Read persistence and split into feature vector
"""

def low(j, matrix):
    col = matrix[:,j]
    col_len = len(col)
    for i in range((col_len-1) , -1, -1): 
        if col[i] == 1:
            return i
    return -1

def read_persistence_index(reduced_matrix):  
    persistence_index = []
    for j in range(reduced_matrix.shape[1]): 
        low_j = low(j, reduced_matrix)
        if low_j == -1:
            start = [j, -1]
            persistence_index.append(start) 
        else:
            feature = persistence_index.index([low_j, -1]) 
            persistence_index[feature][1] = j 
            
    return persistence_index

def read_persistence_interval(persistence_index, simplex, filtration):
    Barcode = []
    for k in persistence_index:
        start = k[0]
        end = k[1]
        d = len(simplex[start])-1
        if start <= len(filtration)//2:
            birth = filtration[start]
            death = filtration[end]
            Barcode.append((d, [birth, death]))
        
    return Barcode

def split_persistence_to_feature_vector(persistence, n, num):
    betti0 = [bar[1] for bar in persistence if bar[0]==0]
    betti1 = [[bar[1][0], 2*n-bar[1][1]] for bar in persistence if bar[0]==1]
    
    for b in betti0:
        if b[1] == 0:
            b[1] = n
        else:
            continue
    
    temp1 = np.linspace(0, n, num+1)
    inter_one = [(temp1[i], temp1[i+1]) for i in range(num)]
    c1 = [0] * num
    for i in range(num):
        for _, death in betti0:
            if death > inter_one[i][0]:
                c1[i] += 1
                
    c2 = [0] * num * 2 
    temp2 = np.linspace(0, 2*n, 2*num+1)
    inter_two = [(temp2[j], temp2[j+1]) for j in range(2*num)]
    for k in range(2*num):
        for birth, death in betti1:
            if (birth <= inter_two[k][1]) and (death >= inter_two[k][0]):
                c2[k] += 1    

    feature_vector = c1 + c2
    
    return feature_vector
    
    
def compute_interaction_feature_matrix(point1, point2, filtration_value): 
    point_cloud = {}
    point_number = {}
    for k1, v1 in point1.items():
        for k2, v2 in point2.items():
            point_number[k1+k2] = [len(v1), len(v2)]
            if (len(v1)==0) or (len(v2)==0):
                point_cloud[k1+k2] = np.array([])
            else:
                T = np.concatenate((v1, v2), axis=0)
                point_cloud[k1+k2] = T
                           
    Feature = {}
    for key, point in point_cloud.items():
        num = point_number[key]
        if (num[0]==0) or (num[1]==0):
             Feature[key] = [0] * 50 * 3     # 50 is the number of bins
        else:
            m = num[0]
            distance_matrix = squareform(pdist(point))
            distance_matrix[:m,:m] = float("inf")
            distance_matrix[m:, m:] = float("inf")
            rips_complex = gd.RipsComplex(distance_matrix=distance_matrix)
            simplex_tree = rips_complex.create_simplex_tree(max_dimension=2)
            
            simplex1 = []
            filt1 = []
            for fil_value in simplex_tree.get_filtration():
                if fil_value[1] <= filtration_value:                       
                    simplex1.append(fil_value[0])
                    filt1.append(fil_value[1])
                        
            complex1 = [simplex1, filt1]
            complex2 = construct_filtration_over_join(complex1=complex1)
            simplicial_complex = complex1[0] + complex2[0]
            filter_value = complex1[1] + complex2[1]
                
            bm = get_boundary_matrix(complex1=complex1[0], complex2=complex2[0])
            mm = matrix_reduction(bm)
            cc = show_matrix(mm) 
                
            persistence_index = read_persistence_index(reduced_matrix=cc)
            persistence = read_persistence_interval(persistence_index=persistence_index,
                                    simplex=simplicial_complex, filtration=filter_value)
            vector = split_persistence_to_feature_vector(persistence=persistence, n=filtration_value, num=50)
            
            Feature[key] = vector  
            
    Feature_vector = []
    for value in Feature.values():
        Feature_vector.extend(value)
            
    return Feature_vector 
     

filename = '/home/yxwang/Metalloprorein/MetallProtein_data/data.xlsx'
df = pd.read_excel(filename)
metalloprotein = []

n = len(df["PDB code"])
for i in range(n):
    name = df["PDB code"][i]
    metal_atom = df["metal_atom_type"][i]
    affinity_data = df["pKd pKi pIC50"][i]
    group = df["group"][i]
    file_type = df["Subset"][i]
    
    temp = [name, file_type, affinity_data, group, metal_atom]
    metalloprotein.append(temp)

root_file1 = "/home/yxwang/Metalloprorein/v2020-other-PL"
root_file2 = "/home/yxwang/Metalloprorein/refined-set"

general_file = []
refined_file = []
for line in metalloprotein:
    if line[1] == "general":
        filename1 = root_file1 + '/' + line[0] + '/' + line[0] + '_pocket.pdb'
        filename2 = root_file1 + '/' + line[0] + '/' + line[0] + '_ligand.mol2'
        temp = [filename1, filename2, line[3], line[0], line[2]]
        general_file.append(temp)
    elif line[1] == "refined":
        filename3 = root_file2 + '/' + line[0] + '/' + line[0] + '_pocket.pdb'
        filename4 = root_file2 + '/' + line[0] + '/' + line[0] + '_ligand.mol2'
        temp = [filename3, filename4, line[3], line[0], line[2]]
        refined_file.append(temp)
        
"""  
Read  train set filename and test set filename 
"""

file = general_file + refined_file 
train_filename = []
test_filename = []
for content in file:
    if content[2] == "train":
        train_filename.append(content)
    elif content[2] == "test":
        test_filename.append(content)
        
print(len(train_filename), len(test_filename))

def main(filename, cutoff_distance, bins_num):
    Feature_matrix = {}

    c = 0
    for file in filename:
        c += 1
        print(c)

        P = extract_pocket_protein_coordination_within_cutoff_distance(filename1=file[0],
                                               filename2=file[1], cutoff=cutoff_distance)    # 5 groups
        M = extract_pocket_metal_ions_coordination(filename=file[0])                         # 10 groups
        L = extract_ligand_coordination(filename=file[1])                                    # 10 groups
        
        # Interaction between metal ions and protein atoms, protein atoms and ligand atoms, metal ions and ligand  atoms
        MP = compute_interaction_feature_matrix(point1=M, point2=P, filtration_value=12.0)
        PL = compute_interaction_feature_matrix(point1=P, point2=L, filtration_value=12.0)
        ML = compute_interaction_feature_matrix(point1=M, point2=L, filtration_value=12.0)
        
        Feature_matrix[file[3]] =  MP + PL + ML + [file[4]]
         
    return Feature_matrix 


Train = main(filename=train_filename, cutoff_distance=10.0, bins_num=50)

# Resort the feature matrix into file

F = []
L = []
for key, value in Train.items():
  F.append(value[:-1])
  L.append(value[-1])


Feature_matrix = {"Feature": F, "Label": L}
df_train = pd.DataFrame(Feature_matrix)
df_train.to_csv('train.csv')








