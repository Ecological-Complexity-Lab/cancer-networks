
#! /usr/bin/Python
import MultiTensor as mt
import numpy as np
import networkx as nx
#from argparse import ArgumentParser
import sys
import itertools
import tools as tl
import time


def main1(adj,layerL,k):
	inf=10000000000000
	err_max=0.0000001

	folder="../data/"
	if(0==True):A=[ nx.MultiGraph() for l in range(layerL) ]   # list of graphs
	else:A=[ nx.MultiDiGraph() for l in range(layerL) ]   # list of graphs
	print folder, adj
        tl.read_graph(folder,adj,A)
	print "Undirected=",bool(0)
	print "Assortative=",bool(0)
	tl.print_graph_stat(A,0)

	if(0):tl.out_graph(folder,A)

	if(0==True): 
		u_list=v_list=tl.remove_zero_entries_undirected(A) 
	else:
		u_list=tl.remove_zero_entries_u(A)   # list of INT INDECES of nodes with zero out degree
		v_list=tl.remove_zero_entries_v(A)   # list of INT INDECES of nodes with zero in degree

	MT=mt.MultiTensor(  N=A[0].number_of_nodes(),
			L=layerL,K=k, 
			N_real=1,
			tolerance=0.1,
			decision=2,
			maxit=500,
			rseed=0,
			out_adjacency=bool(0),
			inf=inf,
			err_max=err_max,
			err=0.1,
			initialization=0,
			undirected=bool(0),
			folder=folder,
			end_file='.dat',
			adj=adj,
			w_file='w.dat',
			assortative=bool(0),
			name=adj[:-4]
			)

	tic = time.clock()
	N=A[0].number_of_nodes()
	B=np.empty(shape=[layerL,N,N])

	for l in range(layerL):B[l,:,:]=nx.to_numpy_matrix(A[l],weight='weight')
	
	MT.cycle_over_realizations(A,B,u_list,v_list)		

	#tl.print_graph_stat(A)	

	toc = time.clock()	  
	print "It took ",toc-tic," seconds."; 	

if __name__ == "__main__":
    layers= ['badnews','closefriend','empathetic', 'feelpositive','goodnews',
             'social advice', 'spendtime', 'support' ]
    
    pair = list(itertools.product(layers,layers))
    pair = [list(x) for x in pair]
    for item in pair:
        if item[0]==item[1]:
            pair.remove(item)
    
    three_pair = list (itertools.permutations(layers,3))
    three_pair = [list(x) for x in three_pair]
    
    
    for item in pair:
        main1(str(item[0])+"_"+str(item[1])+'.dat' , 2)


#for item in layers:
#    main1(str(item)+'.dat')

#alldorm = []
#for i in range(len(layers)):
#    a=[layers[i]]+[x for x in layers if x not in [layers[i]]]
#    alldorm.append(a)
#
#for item in alldorm:
#    name = ""
#    for i in item:
#        name = name + str(i)+"_"
#    
#    main1(name+'20.dat')
