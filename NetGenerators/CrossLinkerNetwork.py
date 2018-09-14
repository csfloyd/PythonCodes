#!/usr/bin/python3

import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import networkx as nx
from collections import OrderedDict
import csv

# read in snapshot, dissipation, tensions

args = sys.argv

#within_path = args[1]
#fname = args[2]

within_path = 'MH/Trial4/'
fname = 'MH4'
pathname = '/Users/csfloyd/Dropbox/DissipationProject/Tdir/' + within_path + 'Output/'

dissipation = open(pathname + 'dissipation.traj','r')
snapshot = open(pathname + 'snapshot.traj','r')


diss_lines = dissipation.readlines()
snap_lines = snapshot.readlines()

num_time_points = int(len(diss_lines)/2)

num_time_points = 1000


def parse_pos_line(pos_line):
    
    f_line = pos_line.split()
    num_beads = int(len(f_line)/3)
    bead_list = []
    for i in range(num_beads):
        
        bead_list.append(np.array([float(f_line[3*i]),\
                                   float(f_line[3*i+1]),\
                                   float(f_line[3*i+2])]))

    return bead_list


def check_bead_in_fil_dict_srt(bead, fil_dict):
    tol = 3
    real_dist = 10.8
    bead_com = c_o_m(bead)

    fd_sorted = sorted(fil_dict.items(), key=lambda item: \
                       np.linalg.norm(bead_com - item[1][1]))
    
    for fil in fd_sorted:
        fil_beads = fil[1][0]
        
        for bead_fil in fil_beads:
            
            dist = np.linalg.norm(bead_fil-bead)
            
            if np.abs(dist-real_dist) < tol:
                
                return int(fil[0])


          
def c_o_m(list_of_beads):
    
    c_o_m = np.array([0.0,0.0,0.0])
    
    for bead in list_of_beads:
        c_o_m += bead

    return c_o_m / len(list_of_beads) 


def safe_inverse(mat):
    (r,c) = np.shape(mat)
    new_mat = np.zeros([r,c])
    for i in range(c):
        for j in range(c):
            num = float(mat[i][j])
            if num==0:
                new_mat[i][j] = 0
            else:
                new_mat[i][j] = 1/num
    return new_mat   

def booleanize(mat):
    (r,c) = np.shape(mat)
    new_mat = np.zeros([r,c])
    for i in range(c):
        for j in range(c):
            num = float(mat[i][j])
            if num==0:
                new_mat[i][j] = 0
            else:
                new_mat[i][j] = 1
    return new_mat


# list of adjacency matrices and linker dicts
fil_matrix = None
fil_matrix_list = []
linker_list_of_dicts = []
data_list_of_dicts = []
list_of_graphs = []
i=0
for j in range(0,num_time_points):

    # get line of this timestep
    
    line=snap_lines[i].split()
    
    step_num = int(line[0])
    time_point = float(line[1])
    num_filaments = int(line[2])
    num_linkers = int(line[3])
    num_motors = int(line[4])
    num_branchers = int(line[5])
    cell_length = 2 + 2 * num_filaments + 2 * num_linkers + \
                  2 * num_motors + num_branchers

    #change this, initialize it to previous matrix
    if j==0:
        fil_matrix = np.zeros([num_filaments,num_filaments])
    

    # if there are no linkers the network is empty
    if num_linkers == 0:
        
        fil_matrix_list.append(np.zeros([num_filaments,num_filaments]))
        
    else:

        # this dict contains the id, bead position arrays, and com of the filament
        fil_dict = {}

        for n in range(0,num_filaments):
            
            fil_line = i + 1 + 2*n
            fil_num = snap_lines[fil_line].split()[1]
            pos_line = fil_line + 1
            list_of_beads = parse_pos_line(snap_lines[pos_line])
            fil_dict[fil_num] = (list_of_beads, c_o_m(list_of_beads))



        # this set contains the linker ids in this timestep
        list_of_linker_ids = set()
        
        for n in range(0,num_linkers):
            
            linker_line = i + 2*num_filaments  + 2*n + 1
            linker_id = snap_lines[linker_line].split()[1]
            list_of_linker_ids.add(linker_id)

        # thsi set contains the linker ids contained in the list of dicts previously
        ids_in_list_of_dict = set()
        
        for n in linker_list_of_dicts:
            ids_in_list_of_dict.add(list(n.keys())[0])

            


        # these sets contain the linker ids of unbound and newly bound linkers
        
        gone_ids = ids_in_list_of_dict - list_of_linker_ids
        new_ids = list_of_linker_ids - ids_in_list_of_dict

 
        # for each linker, check what to do with it then do it            
        for n in range(0,num_linkers):

            # get the linker id
            linker_line = i + 2*num_filaments  + 2*n + 1
            linker_id = snap_lines[linker_line].split()[1]

            
            # these are new, add them to the adjacency matrix and the list of dicts
            if linker_id in new_ids:
                
                pos_line = linker_line + 1
                beads = parse_pos_line(snap_lines[pos_line])

                beadL = beads[0]
                beadR = beads[1]
            
                # get the filaments this linker is bound to
                L_ind = check_bead_in_fil_dict_srt(beadL,fil_dict)
                R_ind = check_bead_in_fil_dict_srt(beadR,fil_dict)


                if L_ind and R_ind:
                    fil_matrix[L_ind][R_ind] += 1
                    fil_matrix[R_ind][L_ind] += 1

                    # add this linker to the dict
                    linker_list_of_dicts.append({linker_id : [L_ind, R_ind]})

        # these linkers need to be removed from the dict and the adjacency matrix updated

        for linker_id in gone_ids:

            temp_list = linker_list_of_dicts
            
            for linker_dict in linker_list_of_dicts:
 
                
                if linker_id == list(linker_dict.keys())[0]:
                    

                    inds = list(linker_dict.values())[0]
                    
                    if inds:
                        fil_matrix[inds[0]][inds[1]] -= 1
                        fil_matrix[inds[1]][inds[0]] -= 1
                    
                    linker_list_of_dicts.remove(linker_dict)

                    
                    
                
        # update the adjacency matrix list
        fil_matrix_list.append(np.array(fil_matrix))
        
        # analyze network on the fly
        G = nx.from_numpy_matrix(fil_matrix)
        G.edges(data=True);
        
        # create a list of the graphs for visualization
        list_of_graphs.append(G)

        
        # create inverse distance and booleanized graphs
        inv_mat = safe_inverse(fil_matrix)
        bool_mat = booleanize(fil_matrix)
        G_inv = nx.from_numpy_matrix(inv_mat)
        G_bool = nx.from_numpy_matrix(bool_mat)
        
        
        # average closeness centrality of the inverse distance graph
        try:
            closeness_cent = np.mean(list(nx.closeness_centrality(G_inv,distance = 'weight').values()))
        except:
            closeness_cent = 0
        # average eigenvector centrality of boolean graph
        try:
            e_cent = np.mean(list(nx.eigenvector_centrality(G_bool).values()))
        except:
            e_cent = 0
        
        # average clustering, weighted
        try: 
            clustering = nx.average_clustering(G_inv)
        except:
            clustering = 0
            
        # assortativity
        try:
            assort = nx.degree_assortativity_coefficient(G_bool)
        except:
            assort = 0
        
        # diameteter
        try:
            diam = nx.diameter(max(nx.connected_component_subgraphs(G_bool), key=len))
        except:
            diam = 0
        
        # density
        try:
            dens = nx.density(G_bool)
        except:
            dens = 0
        
        # number connected components

        n_c_c = nx.number_connected_components(G_bool)
        


        # number connected components
  
        n_connec = nx.average_node_connectivity(G_bool)
        
        # weighted density
        mat_sum = np.sum(fil_matrix)
        weighted_dens = mat_sum/(num_filaments*(num_filaments-1))
        
        # average node strength
        a_n_s = mat_sum / (num_filaments)
        
        data_list_of_dicts.append({'t' : j, 'averageClustering' : clustering, 'closenessCentrality' : closeness_cent, \
                                 'eigenvectorCentrality' : e_cent, 'assortativity': assort, 'diameter': diam, \
                                  'density': dens, 'weightedDensity': weighted_dens, 'meanNodeStrength': a_n_s, \
                                  'connectedComponents': n_c_c,'meanNodeConnectivity': n_connec})
        


    # update the cell length
    i += cell_length



"""
plt.figure()
sns.heatmap(fil_matrix_list[-1])

plt.figure()
sns.heatmap(fil_matrix_list[10])



G = nx.from_numpy_matrix(fil_matrix_list[-1])

G.edges(data=True);

plt.figure()

pos=nx.circular_layout(G)
nx.draw_networkx(G,pos)

all_weights = []
    #4 a. Iterate through the graph nodes to gather all the weights
for (node1,node2,data) in G.edges(data=True):
    all_weights.append(data['weight']) #we'll use this when determining edge thickness
 
    #4 b. Get unique weights
unique_weights = list(set(all_weights))
 
    #4 c. Plot the edges - one by one!
for weight in unique_weights:
        #4 d. Form a filtered list with just the weight you want to draw
    weighted_edges = [(node1,node2) for (node1,node2,edge_attr) in G.edges(data=True) if edge_attr['weight']==weight]
        #4 e. I think multiplying by [num_nodes/sum(all_weights)] makes the graphs edges look cleaner
    width = weight*num_filaments*3.0/sum(all_weights)
    nx.draw_networkx_edges(G,pos,edgelist=weighted_edges,width=width)
 
    #Plot the graph
plt.axis('off')

plt.show()

"""


# list of observables
dict_keys = list(data_list_of_dicts[0].keys())

# each observable has a list with its values in time
for i in range(len(dict_keys)):
    
    vars()["list_"+dict_keys[i]] = []
    
for dict_item in data_list_of_dicts:
    
    for i in range(len(dict_keys)):
    
        vars()["list_"+dict_keys[i]].append(dict_item[dict_keys[i]])
    


num_windows = 5
num_time_points = int(len(diss_lines)/2)
start_of_ss = 200
window_length = np.floor((num_time_points - start_of_ss)/num_windows)

# first element in window average lists is name of observables
for i in range(len(dict_keys)):
    vars()["list_"+dict_keys[i]+"av"] = [str(dict_keys[i])]

# window averaging
for n in range(num_windows):
    beginning = int(start_of_ss + (n-1)*window_length)
    end = int(start_of_ss + (n)*window_length)
    
    for i in range(len(dict_keys)):
        vars()["list_"+dict_keys[i]+"av"].append(np.mean(vars()["list_"+dict_keys[i]][beginning:end]))
        
# top row is window names
title = ["Name"]
for n in range(num_windows):
    title.append(fname+'W'+str(n+1))


# write the observables, except for time, to a .csv
path = "/Users/csfloyd/Dropbox/DissipationProject/Analysis/DissPrediction/CLNOut/"
with open(path + fname + 'XLN'+'.csv', 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(title)
    for i in range(len(dict_keys)):
        if dict_keys[i]!='t':
            writer.writerow(vars()["list_"+dict_keys[i]+"av"])



        
# close files

dissipation.close()
snapshot.close()


