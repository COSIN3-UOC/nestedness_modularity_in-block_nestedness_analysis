#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar  1 14:55:14 2018

@author: mariapalazzi
"""
import numpy as np
import pandas as pd
from nestedness_metrics_other_functions import from_edges_to_matrix
import extremal_bi
import extremal_uni
import multiprocessing as multi
import glob, os
import sys
#%%
def structural_network_analysis(fname):
    '''
    function to perform structural anaylsis in binary unipartite and bipartite networks
    by means of nestedness (as defined in ASR et al, PRE 2018), in-block nested and modularity.
        
    The optimization of modularity and in-block nestedness is done employing the extremal optimization
    algorithm.

    Inputs:
    ----------
       fname: list
       
           fname[0]: name of the network file to read
           fname[1]: boolean to indicate if "filename" is a bipartite (True) or 
           unipartite (False) network
           fname[2]: boolean indicating the format of the data file. Three-column
           or edge list (True) or matrix format (False)
         
    
    '''
    name="results_"+str(os.path.basename(fname[0]).split('.csv')[0])+".npz"
    if len(glob.glob(name))==0:
        if fname[2]==True: #wheter the data is edge list or adajcency matrix
            M=from_edges_to_matrix(fname[0],fname[1])
        else:
            M = np.loadtxt(fname[0],dtype='int',delimiter=',')
            
        '''starting the structural analysis '''
        print('starting the structural analysis', str(os.path.basename(fname[0]).split('.csv')[0]))
        
        
        if fname[1]==True: #if the network is bipartite or not
            cols_degr=M.sum(axis=0)
            row_degr=M.sum(axis=1)
            R,C=M.shape #rows and cols
            #Nestednes
            # In-block nestedness with B=1
            Cn_=[np.repeat(1, R),np.repeat(1, C)]
            max_blockN=max(max(Cn_[0]),max(Cn_[1]))+1
            lambdasN=extremal_bi.call_lambda_i(M,cols_degr,row_degr,Cn_[1],Cn_[0],max_blockN,True)
            nestedness_=extremal_bi.calculate_Fitness(M,cols_degr,row_degr,lambdasN[0],lambdasN[1],True)
        
            #Modularity Extremal
            C_=extremal_bi.recursive_step(M,cols_degr,row_degr,.7,3,False)
            max_blockQ=max(max(C_[0]),max(C_[1]))+1
            lambdasQ=extremal_bi.call_lambda_i(M,cols_degr,row_degr,C_[1],C_[0],max_blockQ,False)
            Q_=extremal_bi.calculate_Fitness(M,cols_degr,row_degr,lambdasQ[0],lambdasQ[1],False)
    
            # Inblock nestedness extremal
            Ci_=extremal_bi.recursive_step(M,cols_degr,row_degr,.7,3,True)
            max_blockI=max(max(Ci_[0]),max(Ci_[1]))+1
            lambdasI=extremal_bi.call_lambda_i(M,cols_degr,row_degr,Ci_[1],Ci_[0],max_blockI,True)
            I_=extremal_bi.calculate_Fitness(M,cols_degr,row_degr,lambdasI[0],lambdasI[1],True)
            
        else:
            cols_degr=M.sum(axis=0)
            row_degr=M.sum(axis=1)
            R,C=M.shape #rows and cols
            #Nestednes
            # IBN with B=1
            Cn_=np.repeat(1, C).tolist()
            max_blockN=max(Cn_)+1
            lambdasN=extremal_uni.call_lambda_i(M,cols_degr,Cn_,max_blockN,True)
            nestedness_=extremal_uni.calculate_Fitness(M,cols_degr,lambdasN,True) #in-block nestedness value
            
            # Modularity
            C_=extremal_uni.recursive_step(M,cols_degr,.7,3,False) # vector with labels of partitions
            max_blockQ=max(C_)+1
            lambdasQ=extremal_uni.call_lambda_i(M,cols_degr,C_,max_blockQ,False)
            Q_=extremal_uni.calculate_Fitness(M,cols_degr,lambdasQ,False) # modularity value

            # Inblock nestedness
            Ci_=extremal_uni.recursive_step(M,cols_degr,.7,3,True) # vector with labels of partitions
            max_blockI=max(Ci_)+1
            lambdasI=extremal_uni.call_lambda_i(M,cols_degr,Ci_,max_blockI,True)
            I_=extremal_uni.calculate_Fitness(M,cols_degr,lambdasI,True) #in-block nestedness value
            
            
        ''' Saving results of analysis'''
        print('saving results for', str(os.path.basename(fname[0]).split('.csv')[0]))
        dfq=pd.DataFrame({'rows': pd.Series(C_[0]), 'cols': pd.Series(C_[1])})
        dfi=pd.DataFrame({'rows': pd.Series(Ci_[0]), 'cols': pd.Series(Ci_[1])})
        dfq.to_csv("modularity_partitions_"+str(os.path.basename(fname[0]).split('.csv')[0])+".csv",index=False,float_format='%.0f')
        dfi.to_csv("in-block_partitions_"+str(os.path.basename(fname[0]).split('.csv')[0])+".csv",index=False,float_format='%.0f')
        np.savez_compressed("results_"+str(os.path.basename(fname[0]).split('.csv')[0])+".npz", N=nestedness_,Q=Q_,I=I_)
#%%
def str_to_bool(s):
    if s == 'True':
         return True
    elif s == 'False':
         return False
    else:
         raise ValueError
#%%
def arguments_list_to_pool(argv0,argv1,argv3):
    ''' 
    Function that generate the list of lists with the arguments needed to perfomr structural 
    analysis of different networks in parallel processes.
    
    inputs:
    ----------
    args[0]: 
        directory where the network file to read are
   args[1]: 
       boolean to indicate if "filename" is a bipartite (True) or 
       unipartite (False) network
   args[2]: 
       boolean indicating the format of the data file. Three-column
       or edge list (True) or matrix format (False)
    '''  
    path = str(argv0)
    filenames = sorted(glob.glob(path+"*.csv"))
    bipartite=bool(argv1)
    edges=bool(argv3)
    parameters=list()
    [parameters.append([filenames[i],bipartite,edges]) for i in range(len(filenames))]
    return parameters
#%%
''' 
To perform parallel analysis:
This function will split the list of containing the B parameter, the format of the data files 
and the network type (uni or bi) and call the main function to generate several networks and perfom the structural analysis

nc= number of simultaneous processes 
'''
if __name__ == '__main__':
#    print('parameters')
    parameters=arguments_list_to_pool(sys.argv[1],str_to_bool(sys.argv[2]), str_to_bool(sys.argv[3]))
    n_cpus = multi.cpu_count()
    if n_cpus > 3:
        nc = n_cpus - 1
    else:
        nc = 1
    pool=multi.Pool(processes=nc)
    pool.map(structural_network_analysis,parameters)
    pool.terminate()
    del(pool)
    
    filenames = sorted(glob.glob("results_*.npz"))
    N=[]
    Q=[]
    I=[]
    fi=[]
    for f in filenames:
        ff=(os.path.basename(f).split('.npz')[0])
        fi.append((os.path.basename(ff).split('results_')[1]))
        data=np.load(f)
        N.append((data['N']))
        Q.append((data['Q']))
        I.append((data['I']))
    
    df=pd.DataFrame()
    df['name']=fi
    df['N']=N
    df['Q']=Q
    df['I']=I
    df.to_csv("data_structures_NQI_results.csv",index=False,sep=',')
    for f in filenames:
        os.remove(f)
