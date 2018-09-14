import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy
from scipy import signal as sg
import networkx as nx
import collections


pathname = '/Users/carlosfloyd/Dropbox/dirs/dirBS2/Output/'
cmgraph = open(pathname + 'CMGraph.traj','r')
cmgraph_lines = cmgraph.readlines()
num_time_points = int(len(cmgraph_lines)/3)
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



time_point = '1000'
G = nx.parse_edgelist(dict_of_edge_lines[time_point])

num_filaments = len(G)


plt.figure()

pos=nx.spring_layout(G)
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


degree_sequence = sorted([d for n, d in G.degree(weight=True)], reverse=True)  # degree sequence
# print "Degree sequence", degree_sequence
degreeCount = collections.Counter(degree_sequence)
deg, cnt = zip(*degreeCount.items())

fig, ax = plt.subplots()
plt.bar(deg, cnt, width=0.80, color='b')

plt.title("Degree Histogram")
plt.ylabel("Count")
plt.xlabel("Degree")
ax.set_xticks([d + 0.4 for d in deg])
ax.set_xticklabels(deg)

plt.show()
