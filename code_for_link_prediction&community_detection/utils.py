# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 15:54:37 2020

@author: hexie
"""
import random
import numpy as np
import pandas as pd
import os 
import networkx as nx
from sklearn.metrics import roc_curve
from sklearn.metrics import roc_auc_score
from sklearn.metrics import average_precision_score
from sklearn.metrics import recall_score
from matplotlib import pyplot
import itertools
#os.chdir(r"C:\Users\hexie\Desktop\MURI_LINK\MultiTensor-master\MultiTensor-master\data") 
#path = "\nas\longleaf\home\hexie\MURI\MultiTensor-master\data\\"
path = r"../data/"
#directory1 =  "C:\Users\hexie\Desktop\MURI_LINK\data\Stanford_dorm\Stanford_dorm"

def create_edgepair(data_dict):
    edgepair=[]
    for item in data_dict.values():
        tmp = zip ([int(x) for x in item[:,0]],[int(x) for x in (item[:,1])])
        if not edgepair:
            edgepair = tmp
        else:
            edgepair = edgepair + tmp
        # number of edges in all layers added together,4145
        # edge_count = len(edgepair)
        # total edge pair count (different layers count as one),1759
        edgepair= list(set(edgepair))
    return edgepair

#GIVEN A DATA DICT, EDGE PAIR AND EDGE PAIR COUNT, CREATE OUTPUT NUMPY ARRAY
def create_output(data_dict, edgepair, edge_count, order): 
    output= np.array(edgepair)
    list_E= ['E']*edge_count
    list_E = np.array(list_E)
    list_E = list_E.reshape(-1,1)
    output = np.concatenate((list_E,output),axis=1)
    
    for o in order:
        item = data_dict[o]
#    
#    for item in data_dict.values():
        z = np.zeros((edge_count,1), dtype = int)
        tmp = zip ([int(x) for x in item[:,0]],[int(x) for x in (item[:,1])])
        for it in tmp:
            if it in edgepair:
                z[edgepair.index(it)] = 1
        output = np.concatenate((output, z), axis=1)
    return output

#GIVEN A DATA DICT GET NODE IDX
def get_idx(data_dict):
    # Get idx of everynode,194
    idx= []
    for item in data_dict.values():
        idx = np.concatenate((idx, item[:,0]), axis=None)
        idx = np.concatenate((idx, item[:,1]), axis=None)
    idx=np.unique(idx)
    idx = [int(x) for x in idx]
    return idx

def get_idx_from_edge(mylist):
    idx= [item for t in mylist for item in t] 
    idx=np.unique(idx)
    idx = [int(x) for x in idx]
    return idx
    

def get_idmap(nodes_name):
    
    idx= list(range(len(nodes_name)))
    idmap = dict(zip(nodes_name, idx))
    
    return idmap
    


# GET CERTAIN LAYER AS SUBDATA FROM DATA
def get_specified_layer(data_dict, layer_list):
    out = {x: data_dict[x] for x in layer_list}
    return out


def remove_edges(graph):
    g = graph.copy()
    nodes_to_remove = []
    for n in g.nodes():
        if g.in_degree(n) == 0 and g.out_degree(n) == 0:
            nodes_to_remove.append(n)
    g.remove_nodes_from(nodes_to_remove)
    return g

def remain_layer(idx,info,lay):
    idx = [int(i) for i in idx]
    inf= info.copy()
    q = inf[lay]
    count = 0
    qq = []
    deleted = []
    for item in q:    
        if (int(item[0]) in idx) and (int(item[1]) in idx):
            qq.append(item)
        else:
            deleted.append(item[0])
            deleted.append(item[1])
            count = count+1
    #print("counted delete these edges",count )
    inf[lay] = np.array(qq)
    deleted = np.unique(deleted)
    return inf, deleted


# GET THE FINAL OUTPUT FILE TO BECOME A GRAPH
def get_output_matrix(single_layer,k,r,prefix):
    #path="\nas\longleaf\home\hexie\MURI\MultiTensor-master\data\\"    

    
    name = str(single_layer)
    
    with open(path+name+ prefix + "_" + str(k) + "_" + str(r)+'_u_K'+str(k)+'.dat') as f:
        lines = (line for line in f if not line.startswith('#'))
        uk5 = np.loadtxt(lines, delimiter=' ', skiprows = 0 )
    uk5= uk5[:,1:]
    uk5 = np.array(uk5)
    
    with open(path+name+ prefix + "_" + str(k) + "_" + str(r)+'_v_K'+str(k)+'.dat') as f:
        lines = (line for line in f if not line.startswith('#'))
        vk5 = np.loadtxt(lines, delimiter=' ', skiprows = 0 )
    vk5= vk5[:,1:]
    vk5 = np.array(vk5)
    
    w5={}
    try:
        wk5 = pd.read_csv(path+name+ prefix + "_" + str(k) + "_" + str(r) +'_w_K'+str(k)+'.dat', sep = " ")
        count = 0
        for i in range(wk5.shape[0]):
            if wk5.iloc[i][0] == "a=":
                w5[count]=wk5.iloc[i+1:i+k+1,0:k ]
                count = count+1
                i = i+k+1
                
    except:
        w5 = {}
        wk5 = pd.read_csv(path+name+ prefix + "_" + str(k) + "_" + str(r) +'_w_K'+str(k)+'.dat',index_col=0 )
        
        wk5 = wk5.index.values
        wk5 = [x for x in wk5 if "a=" not in x]
        
        wk5 = np.array(wk5)
        w51 = []
        for lst in wk5:
            w51.append([float(s) for s in lst.split(' ')])
        w51 = np.array(w51)
        
        count = 0
        for i in range(w51.shape[0]):
            tmp = w51[i:i+k, :]
            w5[count] = tmp
            i = i+k+1
            count = count+1
        
    
    #M = np.matmul(np.matmul(np.array(uk5).astype(float).T, np.array(vk5).astype(float)).astype(float), np.array(w5[0]).astype(float)) 
    M = np.matmul(np.matmul(np.array(uk5), np.array(w5[0],
                            dtype='float')), np.array(vk5).T)

    return M
    
def get_roc(name, p,q):
    # calculate scores
    ns_auc = roc_auc_score(p,q)
    # summarize scores
    #print(name, 'ROC AUC=%.3f' % (ns_auc))
    # calculate roc curves
    fpr, tpr, thresholds = roc_curve(p, q)
    # plot the roc curve for the model
    pyplot.plot(fpr, tpr, marker='.', label=name + str (ns_auc))
    # axis labels
    pyplot.xlabel('False Positive Rate')
    pyplot.ylabel('True Positive Rate')
    # show the legend
    pyplot.legend()
    # show the plot
    #pyplot.savefig(str(single_layer) + "_solo.png")
    pyplot.show()
    return ns_auc



def read_single_true(data, single,k,r,prefix):
    #path="C:\Users\hexie\Desktop\MURI_LINK\MultiTensor-master\MultiTensor-master\data\\"
    name = str(single)
    print("layer is :", name)
    layer_info = get_specified_layer(data, [single])
    edge_in_layer = create_edgepair(layer_info)
    node_in_layer = get_idx(layer_info)
    print("it originally has:",len(edge_in_layer), "edges")
    
    M = get_output_matrix(name,k,r, prefix)
    removed_edges = np.load(path+ name+ prefix + "_" + str(k) + "_" + str(r) +".npy")
    removed_edges = zip ([int(x) for x in removed_edges[:,0]],[int(x) for x in (removed_edges[:,1])])
    
    
    remain_edges = [x for x in edge_in_layer if x not in removed_edges]
    actual_removed_edges = [x for x in edge_in_layer if x in removed_edges]
    print("after deletion it has:",len(remain_edges), "edges")
    
    idmap = get_idmap(node_in_layer)
    
    #####original graph info###########
    G= nx.DiGraph()
    G.add_nodes_from(node_in_layer)
    G.add_edges_from(edge_in_layer)
    A = nx.adjacency_matrix(G)    
    A= A.todense()
    
    
    G2= nx.DiGraph()
    G2.add_nodes_from(node_in_layer)
    G2.add_edges_from(remain_edges)
    G2 = remove_edges(G2)
    
    removed_nodes =[x for x in node_in_layer if x not in G2.nodes()]
    
    for i in removed_nodes:
        j = idmap[i]
        r,c = M.shape
        r = int(r)
        col = np.zeros(r)
        col = col.reshape(r,1)
        M = np.hstack((M[:,:j], col, M[:,j:]))
        r,c = M.shape
        c = int(c)
        row = np.zeros(c)
        row = row.reshape(1,c)
        M = np.vstack((M[:j,:], row, M[j:,:]))
    
    removed_index = [(idmap[x[0]], idmap[x[1]]) for x in removed_edges]
    idx_i, idx_j = zip(*removed_index)
    
    
    
    p = A[idx_i, idx_j]
    p = p.tolist()
    p = p[0]
    q = M[idx_i, idx_j]
    q = list(q)
    
    ns_auc = roc_auc_score(p,q)
    #q1 = np.where(q > 0, 1, 0)
    ns_auc1 = average_precision_score(p,q)
    fpr, tpr, thresholds = roc_curve(p, q)
    return fpr, tpr, ns_auc
    

def read_double_true(data, double,k,r, prefix):
    #path="C:\Users\hexie\Desktop\MURI_LINK\MultiTensor-master\MultiTensor-master\data\\"
    name = str(double[0])+"_"+str(double[1])
    
    layer_info = get_specified_layer(data, [double[0]])
    edge_in_layer = create_edgepair(layer_info)
    node_in_layer = get_idx(layer_info)
    M = get_output_matrix(name,k,r, prefix)
    removed_edges = np.load(path+name+ prefix + "_" + str(k) + "_" + str(r) +".npy")
    removed_edges = zip ([int(x) for x in removed_edges[:,0]],[int(x) for x in (removed_edges[:,1])])
    remain_edges = [x for x in edge_in_layer if x not in removed_edges]
    actual_removed_edges = [x for x in edge_in_layer if x in removed_edges]
    rest = get_specified_layer(data, [double[1]])
    
    
    node_in_rest = get_idx(rest) 
    node_remain = get_idx_from_edge(remain_edges)
    node_in_all = np.unique(node_remain + node_in_rest)
    
    idmap_predict = get_idmap(node_in_layer)
    idmap_all = get_idmap(node_in_all)
    
    
    #####original graph info###########
    G= nx.DiGraph()
    G.add_nodes_from(node_in_layer)
    G.add_edges_from(edge_in_layer)
    A = nx.adjacency_matrix(G)    
    A= A.todense()
    
    
    G2= nx.DiGraph()
    G2.add_nodes_from(node_in_layer)
    G2.add_edges_from(remain_edges)
    G2 = remove_edges(G2)
    
    removed_nodes =[x for x in node_in_layer if x not in node_in_all]
    nodes_not_in_predict = [x for x in node_in_all if x not in node_in_layer]
    
    
    need_deletion = [idmap_all[x] for x in nodes_not_in_predict]
    
    M = np.delete(M, need_deletion, 0)
    M = np.delete(M, need_deletion, 1)  # delete j th column of M

    for i in removed_nodes:
        j = idmap_predict[i]
        r,c = M.shape
        r = int(r)
        col = np.zeros(r)
        col = col.reshape(r,1)
        M = np.hstack((M[:,:j], col, M[:,j:]))
        r,c = M.shape
        c = int(c)
        row = np.zeros(c)
        row = row.reshape(1,c)
        M = np.vstack((M[:j,:], row, M[j:,:]))
    
    removed_index = [(idmap_predict[x[0]], idmap_predict[x[1]]) for x in removed_edges]
    idx_i, idx_j = zip(*removed_index)
        
    p = A[idx_i, idx_j]
    p = p.tolist()
    p = p[0]
    q = M[idx_i, idx_j]
    q = list(q)

    ns_auc = roc_auc_score(p,q)
    #q1 = np.where(q > 0, 1, 0)
    ns_auc1 = average_precision_score(p,q)
    fpr, tpr, thresholds = roc_curve(p, q)
    return fpr, tpr, ns_auc


def write_single(data, item, k,r, prefix):
    name = str(item)

    layer_info = get_specified_layer(data, [item])
    edge_in_layer = create_edgepair(layer_info)
    node_in_layer = get_idx(layer_info)
    print("it originally has:",len(edge_in_layer))
    
    ## take 20% of that edge pair randomly and set any edges appeared in that 20% to be 0
    all_edge_pair = list(itertools.product(node_in_layer,node_in_layer))
    all_edge_pair = [x for x in all_edge_pair if x[0]!=x[1]]
    remove_count = int(round(len(all_edge_pair)*0.5))
    removed_edges = random.sample(all_edge_pair, remove_count)
    removed_edges = [(int(x[0]), int(x[1])) for x in removed_edges]
    
    remain_edges = [x for x in edge_in_layer if x not in removed_edges]
    actual_removed_edges = [x for x in edge_in_layer if x in removed_edges]
    print("after deletion it has:",len(remain_edges))
    out = create_output(layer_info,remain_edges, len(remain_edges),[item])
    out = pd.DataFrame(out)
    
    #path="C:\Users\hexie\Desktop\MURI_LINK\MultiTensor-master\MultiTensor-master\data\\"
    np.save(path+name+ prefix  + "_" + str(k) + "_" + str(r) + '.npy', removed_edges)
    out.to_csv(path+name+ prefix + "_" + str(k) + "_" + str(r) +'.dat', sep =' ',header=False, index=False)


def write_double(data, item,k,r, prefix):
    name = str(item[0])+"_"+str(item[1])
    layer_info = get_specified_layer(data, [item[0]])
    edge_in_layer = create_edgepair(layer_info)
    node_in_layer = get_idx(layer_info)
    print("it originally has:",len(edge_in_layer))
    
    ## take 20% of that edge pair randomly and set any edges appeared in that 20% to be 0
    all_edge_pair = list(itertools.product(node_in_layer,node_in_layer))
    all_edge_pair = [x for x in all_edge_pair if x[0]!=x[1]]
    remove_count = int(round(len(all_edge_pair)*0.5))
    removed_edges = random.sample(all_edge_pair, remove_count)
    removed_edges = [(int(x[0]), int(x[1])) for x in removed_edges]
    
    remain_edges = [x for x in edge_in_layer if x not in removed_edges]
    actual_removed_edges = [x for x in edge_in_layer if x in removed_edges]
    print("after deletion it has:",len(remain_edges))
    
    
    double = get_specified_layer(data, [item[1]])
    
    double_edges = create_edgepair(double)
    double_edges = double_edges+remain_edges
    both = get_specified_layer(data, item)
    
    double_edges = list(set(double_edges))
    
    out = create_output(both, double_edges, len(double_edges), item)
    out = pd.DataFrame(out)
    
    
    #path="C:\Users\hexie\Desktop\MURI_LINK\MultiTensor-master\MultiTensor-master\data\\"
    np.save(path+name+ prefix + "_" + str(k) + "_" + str(r) +'.npy', removed_edges)
    out.to_csv(path+name+ prefix + "_" + str(k) + "_" + str(r) +'.dat', sep =' ',header=False, index=False)
    

def write_multiple(data, item,k,r, prefix):
    
    name = ""
    for i in item:
        name = name + str(i)+"_"
    
    
    
    
    layer_info = get_specified_layer(data, [item[0]])
    edge_in_layer = create_edgepair(layer_info)
    node_in_layer = get_idx(layer_info)
    print("it originally has:",len(edge_in_layer))
    
    ## take 20% of that edge pair randomly and set any edges appeared in that 20% to be 0
    all_edge_pair = list(itertools.product(node_in_layer,node_in_layer))
    all_edge_pair = [x for x in all_edge_pair if x[0]!=x[1]]
    remove_count = int(round(len(all_edge_pair)*0.5))
    removed_edges = random.sample(all_edge_pair, remove_count)
    removed_edges = [(int(x[0]), int(x[1])) for x in removed_edges]
    
    remain_edges = [x for x in edge_in_layer if x not in removed_edges]
    actual_removed_edges = [x for x in edge_in_layer if x in removed_edges]
    print("after deletion it has:",len(remain_edges))
    
    
    multiple = get_specified_layer(data, item[1:])
    
    multiple_edges = create_edgepair(multiple)
    
    print(len(get_idx(multiple)))
    print(len(get_idx_from_edge(multiple_edges)))
    
    multiple_edges = multiple_edges+remain_edges
    mul = get_specified_layer(data, item)
    
    multiple_edges = list(set(multiple_edges))
    

    out = create_output(mul, multiple_edges, len(multiple_edges), item)
    out = pd.DataFrame(out)
    
    
    #path="C:\Users\hexie\Desktop\MURI_LINK\MultiTensor-master\MultiTensor-master\data\\"
    np.save(path+name + prefix + "_" + str(k) + "_" + str(r) +'.npy', removed_edges)
    out.to_csv(path+name+ prefix + "_" + str(k) + "_" + str(r) +'.dat', sep =' ',header=False, index=False)
    
    return name


def write_multiple_community(data, item):
    
    name = "community_"
    for i in item:
        name = name + str(i)+"_"
    
 
    multiple = get_specified_layer(data, item)
    multiple_edges = create_edgepair(multiple)
    multiple_edges = list(set(multiple_edges))
    

    out = create_output(multiple, multiple_edges, len(multiple_edges), item)
    out = pd.DataFrame(out)
    
    
    #path="C:\Users\hexie\Desktop\MURI_LINK\MultiTensor-master\MultiTensor-master\data\\"
    #np.save(path+name+'.npy', removed_edges)
    out.to_csv(path+name+'.dat', sep =' ',header=False, index=False)
    
    return name


    
    
def read_multiple_true(data, item,k,r,prefix):
    #path="C:\Users\hexie\Desktop\MURI_LINK\MultiTensor-master\MultiTensor-master\data\\"
    name = ""
    for i in item:
        name = name + str(i)+"_"
    
    
    layer_info = get_specified_layer(data, [item[0]])
    edge_in_layer = create_edgepair(layer_info)
    node_in_layer = get_idx(layer_info)
    M = get_output_matrix(name,k,r, prefix)
    removed_edges = np.load(path+name+ prefix + "_" + str(k) + "_" + str(r) +".npy")
    removed_edges = zip ([int(x) for x in removed_edges[:,0]],[int(x) for x in (removed_edges[:,1])])
    remain_edges = [x for x in edge_in_layer if x not in removed_edges]
    actual_removed_edges = [x for x in edge_in_layer if x in removed_edges]
    rest = get_specified_layer(data, item[1:])
    
    print(len(node_in_layer))
    node_in_rest = get_idx(rest) 
    node_remain = get_idx_from_edge(remain_edges)
    node_in_all = np.unique(node_remain + node_in_rest)
    
    print(len(node_in_all))
    
    
    idmap_predict = get_idmap(node_in_layer)
    idmap_all = get_idmap(node_in_all)
    
    
    #####original graph info###########
    G= nx.DiGraph()
    G.add_nodes_from(node_in_layer)
    G.add_edges_from(edge_in_layer)
    A = nx.adjacency_matrix(G)    
    A= A.todense()

    
    G2= nx.DiGraph()
    G2.add_nodes_from(node_in_layer)
    G2.add_edges_from(remain_edges)
    G2 = remove_edges(G2)
    
    
    removed_nodes =[x for x in node_in_layer if x not in node_in_all]
    nodes_not_in_predict = [x for x in node_in_all if x not in node_in_layer]
    
    print(A.shape)
    print(M.shape)
    
    print(removed_nodes)
    print(nodes_not_in_predict)
    
    
    
    need_deletion = [idmap_all[x] for x in nodes_not_in_predict]
    
    M = np.delete(M, need_deletion, 0)
    M = np.delete(M, need_deletion, 1)  # delete j th column of M

    for i in removed_nodes:
        j = idmap_predict[i]
        r,c = M.shape
        r = int(r)
        col = np.zeros(r)
        col = col.reshape(r,1)
        M = np.hstack((M[:,:j], col, M[:,j:]))
        r,c = M.shape
        c = int(c)
        row = np.zeros(c)
        row = row.reshape(1,c)
        M = np.vstack((M[:j,:], row, M[j:,:]))
    
    print(A.shape)
    print(M.shape)
    
    removed_index = [(idmap_predict[x[0]], idmap_predict[x[1]]) for x in removed_edges]
    idx_i, idx_j = zip(*removed_index)
    
    print(removed_nodes)
    
    p = A[idx_i, idx_j]
    p = p.tolist()
    p = p[0]
    q = M[idx_i, idx_j]
    q = list(q)

    ns_auc = roc_auc_score(p,q)
    #q1 = np.where(q > 0, 1, 0)
    ns_auc1 = average_precision_score(p,q)
    fpr, tpr, thresholds = roc_curve(p, q)
    return fpr, tpr, ns_auc
