import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy import signal as sg
import networkx as nx
import collections
import copy


pathname = '/Users/csfloyd/Dropbox/DissipationProject/VA3dir/Data/Trial5/Output/'
cmgraph = open(pathname + 'CMGraph.traj','r')
cmgraph_lines = cmgraph.readlines()
num_time_points = int(len(cmgraph_lines)/3)
num_filaments = 50
dict_of_edge_lines = {}

i=1

for j in range(num_time_points):

    time_point = cmgraph_lines[i-1].split()[1]
    line = cmgraph_lines[i].split()
    num_edges = int(len(line)/3)
    edge_list_line = []

    for k in range(num_edges):
        node_1 = str(line[3*k])
        node_2 = str(line[3*k+1])
        weight = str(line[3*k+2])
        entry = node_1 + " " + node_2 + " {'weight':" + weight + "}"
        edge_list_line.append(entry)
        

    dict_of_edge_lines[str(j)] = edge_list_line   
    i = i+3


data_list_of_dicts = []

for j in range(num_time_points):
    
    G = nx.parse_edgelist(dict_of_edge_lines[str(j)])
    G_inv = copy.deepcopy(G)
    G_bool = copy.deepcopy(G)
                          
    for edge in G_inv.edges(data=True):
        u = edge[0]
        v = edge[1]
        w = G_inv[u][v]['weight']
        if w == 0:
            ret = 0
        else:
            ret = 1/w

        G_inv[u][v]['weight'] = ret

    for edge in G_bool.edges(data=True):
        u = edge[0]
        v = edge[1]
        ret = 1
        G_bool[u][v]['weight'] = ret

    # average closeness centrality of the inverse distance graph
    
    try:
        #print(list(nx.closeness_centrality(G_inv,distance = 'weight').values()))
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
    mat_sum = G.size(weight='weight')
    num_nodes = G.number_of_nodes()
    weighted_dens = mat_sum/(num_filaments*(num_filaments-1))
        
    # average node strength
    a_n_s = mat_sum / (num_filaments)
        
    data_list_of_dicts.append({'t' : j, 'averageClustering' : clustering, 'closenessCentrality' : closeness_cent, \
                                'eigenvectorCentrality' : e_cent, 'assortativity': assort, 'diameter': diam, \
                                'density': dens, 'weightedDensity': weighted_dens, 'meanNodeStrength': a_n_s, \
                                'connectedComponents': n_c_c,'meanNodeConnectivity': n_connec})

    print(j)

    
# list of observables
dict_keys = list(data_list_of_dicts[0].keys())

# each observable has a list with its values in time
for i in range(len(dict_keys)):
    
    vars()["list_"+dict_keys[i]] = []
    
for dict_item in data_list_of_dicts:
    
    for i in range(len(dict_keys)):
        
        ret = dict_item[dict_keys[i]]
        
        if(ret!=ret):
            ret = 0.0
           
        vars()["list_"+dict_keys[i]].append(ret) 
    

                          

num_windows = 3
num_time_points = int(len(diss_lines)/2)
start_of_ss = 500
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
path = "/homes/csfloyd/Dropbox/DissipationProject/Analysis/DissPrediction/XLN/"
with open(path + fname + 'XLN'+'.csv', 'w') as csvfile:
    writer = csv.writer(csvfile, delimiter=',')
    writer.writerow(title)
    for i in range(len(dict_keys)):
        if dict_keys[i]!='t':
            writer.writerow(vars()["list_"+dict_keys[i]+"av"])



        
# close files

dissipation.close()
snapshot.close()

    


