#!/usr/bin/python3

import sys
import numpy as np
import networkx as nx
from collections import OrderedDict
import csv
import copy


# read in snapshot, dissipation, tensions

args = sys.argv

#within_path = args[1]
#fname = args[2]

pathname = '/Users/csfloyd/Dropbox/DissipationProject/VA3dir/Data/Trial5/Output/'

dissipation = open(pathname + 'dissipation.traj','r')
snapshot = open(pathname + 'snapshot.traj','r')

# options are ME, PE, or M
node = 'M'

# options are ['Absolute',cutoff] or ['Relative',percentage]
method = ['Relative',.3]


diss_lines = dissipation.readlines()
snap_lines = snapshot.readlines()
num_time_points = int(len(diss_lines)/2)

#num_time_points = int(len(diss_lines)/2)


def c_o_m(list_of_beads):
    
    c_o_m = np.array([0.0,0.0,0.0])
    
    for bead in list_of_beads:
        c_o_m += bead

    return c_o_m / len(list_of_beads)

def distance(pair_of_beads):

    bead_1 = pair_of_beads[0]
    bead_2 = pair_of_beads[1]

    vec = bead_2 - bead_1
    
    return np.linalg.norm(vec)



    
# list of adjacency matrices and linker dicts
matrix = None
matrix_list = []
data_list_of_dicts = []
list_of_graphs = []
i=0



for j in range(num_time_points):

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
    
    if(node == 'ME' or node == 'PE'):
        num_el = num_filaments
        off_set = 2
    else:
        num_el = num_motors
        off_set = 2 + 2*num_filaments + 2*num_linkers
        
    matrix = np.zeros([num_el,num_el])

    

    for el1 in range(num_el):
        
        element_line_1 = i + off_set + 2 * el1
        line_1 = snap_lines[element_line_1].split() 
        line_1 = [float(i) for i in line_1]

        if (node == 'ME'):
            bead_1 = np.array([line_1[0],line_1[1],line_1[2]])
        elif (node == 'PE'):
            bead_1 = np.array([line_1[-3],line_1[-2],line_1[-1]])
        else:
            bead_1 = c_o_m([np.array([line_1[0],line_1[1],line_1[2]]) ,\
                            np.array([line_1[-3],line_1[-2],line_1[-1]])])
       
        for el2 in range(num_el):

            element_line_2 = i + off_set + 2 * el2
            line_2 = snap_lines[element_line_2].split()
            line_2 = [float(i) for i in line_2]

            if (node == 'ME'):
                bead_2 = np.array([line_2[0],line_2[1],line_2[2]])
            elif (node == 'PE'):
                bead_2 = np.array([line_2[-3],line_2[-2],line_2[-1]])
            else:
                bead_2 = c_o_m([np.array([line_2[0],line_2[1],line_2[2]]) ,\
                            np.array([line_2[-3],line_2[-2],line_2[-1]])])

            matrix[el1,el2] += distance([bead_1,bead_2])

    if(method[0]=='Absolute'):
        
        nmatrix = copy.deepcopy(matrix)
        matrix[np.where(nmatrix<=method[1])]=1
        matrix[np.where(nmatrix>method[1])]=0
        
    else:
        
        ordered_dist_list = sorted(set(matrix.flatten()))                              
        l = int(len(ordered_dist_list)*method[1])
        try:
            thresh = ordered_dist_list[l]
        except:
            thresh = 0
        nmatrix = copy.deepcopy(matrix)
        matrix[np.where(nmatrix<=thresh)]=1
        matrix[np.where(nmatrix>thresh)]=0
        
    matrix = matrix - np.identity(num_el)
    matrix_list.append(matrix)

    G = nx.from_numpy_matrix(matrix)
              
        
    
    # average closeness centrality of the inverse distance graph
    try:
        closeness_cent = np.mean(list(nx.closeness_centrality(G,distance = 'weight').values()))
    except:
        closeness_cent = 0
        
    # average eigenvector centrality of boolean graph
    try:
        e_cent = np.mean(list(nx.eigenvector_centrality(G).values()))
    except:
        e_cent = 0
    
    # average clustering, weighted
    try: 
        clustering = nx.average_clustering(G)
    except:
        clustering = 0
        
    # assortativity
    try:
        assort = nx.degree_assortativity_coefficient(G)
    except:
        assort = 0
    
    # diameteter
    try:
        diam = nx.diameter(max(nx.connected_component_subgraphs(G), key=len))
    except:
        diam = 0
    
    # density
    try:
        dens = nx.density(G)
    except:
        dens = 0
    
    # number connected components
    n_c_c = nx.number_connected_components(G)
    
    # number connected components
    n_connec = nx.average_node_connectivity(G)
    
 

    
    data_list_of_dicts.append({'t' : j, 'averageClustering' : clustering, 'closenessCentrality' : closeness_cent, \
                             'eigenvectorCentrality' : e_cent, 'assortativity': assort, 'diameter': diam, \
                              'density': dens, 'connectedComponents': n_c_c,\
                              'meanNodeConnectivity': n_connec})
    


    # update the cell length
    i += cell_length

print('Done.')
