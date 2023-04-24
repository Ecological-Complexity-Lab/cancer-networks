from contextlib import closing
import pickle
import numpy as np
import pandas as pd
from utils import *
from main1 import main1
import itertools
from multiprocessing import Pool, freeze_support

GENES = ["breast", "colon", "head", "kidneyc", "kidneyp", "liver", "lunga", "lungs", "prostate", "stomach", "thyroid","uterine"]
N = 5
KK = range(2,4)

weird_list = []

for i in range(N):
    for k in KK:
        weird_list.append([k,i])


def load_cancer_data():

    layers = ["breast", "colon", "head", "kidneyc", "kidneyp", "liver", "lunga", "lungs", 
              "prostate", "stomach", "thyroid","uterine"]
    data_path = r"../data/"
    data ={}
    
    for lay in layers:
        data[lay] = np.loadtxt(data_path + lay +".txt")
    
    pair = list(itertools.product(layers,layers))
    pair = [list(x) for x in pair]
    for item in pair:
        if item[0]==item[1]:
            pair.remove(item)
            
    allcomb = []
    for i in range(len(layers)):
        a=[layers[i]]+[x for x in layers if x not in [layers[i]]]
        allcomb.append(a)
            
    return data, layers, pair, allcomb

def crazy_pl(k,r):
    data1, layers1, pair1, allcomb1 = load_cancer_data()
    prefix = "_50"

    single_auc = []
    double_auc = []
    all_auc = []

    for la in layers1:
        write_single(data1, la,k,r, prefix)
        main1(la+ prefix + "_" + str(k) + "_" + str(r)+".dat",1,k)
        fpr, tpr, tmp = read_single_true(data1, la,k,r, prefix)
        single_auc.append(tmp)
    for pa in pair1:
        name = str(pa[0])+"_"+str(pa[1])
        write_double(data1, pa,k,r, prefix)
        main1(name+ prefix + "_" + str(k) + "_" + str(r)+".dat",2, k)
        fpr, tpr, tmp = read_double_true(data1, pa, k,r, prefix)
        double_auc.append(tmp)
    for al in allcomb1:
        name = ""
        for j in al :
            name = name + str(j)+"_"
        write_multiple(data1, al,k,r, prefix)
        main1(name+ prefix + "_" + str(k) + "_" + str(r)+".dat",12,k)
        fpr, tpr, tmp = read_multiple_true(data1, al, k,r, prefix)
        all_auc.append(tmp)

    np.savetxt("single_"+str(k)+ prefix + "_" + str(r)+".txt", single_auc)
    np.savetxt("double_"+str(k)+ prefix + "_" + str(r)+".txt", double_auc)
    np.savetxt("all_"+str(k)+ prefix + "_" + str(r)+".txt", all_auc)

def func_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return crazy_pl(*a_b)

def main():
    pool = Pool()
    N = 5
    KK = range(2,3)

    weird_list = []
    for i in range(N):
        for k in KK:
            weird_list.append([k,i])
    pool.map(func_star, weird_list)

freeze_support()
main()
