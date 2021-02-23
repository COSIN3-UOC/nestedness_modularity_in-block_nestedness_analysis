#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 16:16:57 2017

@author: maria
"""
#%%
import numpy as np
#%%
def glob_nestnulb(M):
    '''
    function to calculate the nestedness fitness N, a modified version
    of Ñ, corrected by a null model
    Metric developed by ASR et al, PRE 2018.

    Inputs:n
    ----------
        M: array
            A matrix to which I want to calculate the N
    
    output:
    ----------
    N: number
        The N score for the whole matrix
    '''
    rw,cl=M.shape
    colN=np.zeros((cl,cl))
    rowN=np.zeros((rw,rw))
    cols_degr = M.sum(axis=0) # dregree of the cols nodes
    rows_degr = M.sum(axis=1) # degree of the rows nodes
    
    #Find N col score
    for i in range(cl): # at a left position with respect to column j
      	for j in range(cl):
              if M[i,j]==1:
                  if (cols_degr[i]>=cols_degr[j]) & (cols_degr[j]>0): # heaviside
                      if (cols_degr[i]==cols_degr[j]):
                          colN[i,j]=(np.sum((M[:,i]*M[:,j]),dtype=float)-((cols_degr[i]*cols_degr[j])/rw))/(2*cols_degr[j]) #paired overlap
                      else:
                          colN[i,j]=(np.sum((M[:,i]*M[:,j]),dtype=float)-((cols_degr[i]*cols_degr[j])/rw))/cols_degr[j]
        
    N_COL = (np.sum(colN,dtype=float)/(cl-1))
    
    for i in range(rw): #at an upper position with respect to row j
        for j in range(rw):
            if M[i,j]==1:
                if (rows_degr[i]>=rows_degr[j]) & (rows_degr[j]>0): # Heaviside
                    if (rows_degr[i]==rows_degr[j]):
                        rowN[i,j]=(np.sum((M[i,:]*M[j,:]),dtype=float)-((rows_degr[i]*rows_degr[j])/cl))/(2*rows_degr[j]) #paired overlap
                    else:
                        rowN[i,j]=(np.sum((M[i,:]*M[j,:]),dtype=float)-((rows_degr[i]*rows_degr[j])/cl))/rows_degr[j] #paired overlap

    N_ROW = (np.sum(rowN,dtype=float)/(rw-1))
    
    #Find N
    N=(N_COL+N_ROW)*(2./(rw+cl))
    return N
#%%
def glob_nestnul(M):
    '''
    function to calculate the nestedness fitness N, a modified version
    of Ñ, corrected by a null model
    Metric developed by ASR et al, PRE 2018.

    Inputs:n
    ----------
        M: array
            A matrix to which I want to calculate the N
    
    output:
    ----------
    N: number
        The N score for the whole matrix
    '''
    rw,cl=M.shape
    colN=np.zeros((cl,cl))
    cols_degr = M.sum(axis=0) # dregree of the cols nodes
    
    #Find N col score
    for i in range(cl): # at a left position with respect to column j
      	for j in range(cl):
              if (cols_degr[i]>=cols_degr[j]) & (cols_degr[j]>0): # heaviside
                  if (cols_degr[i]==cols_degr[j]):
                      colN[i,j]=(np.sum((M[:,i]*M[:,j]),dtype=float)-((cols_degr[i]*cols_degr[j])/rw))/(2*cols_degr[j]) #paired overlap
                  else:
                      colN[i,j]=(np.sum((M[:,i]*M[:,j]),dtype=float)-((cols_degr[i]*cols_degr[j])/rw))/cols_degr[j]
    
    N_COL = (np.sum(colN,dtype=float)/(cl-1))
        
    #Find N
    N=(N_COL)*(2./(cl))
    return N
#%%
def from_edges_to_matrix(fname,bipartite=True):
    '''
    function to generate the network adjacency matrix given the edge list

    Inputs:n
    ----------
        fname: data file
            The data file containing the edge lists. accepts .csv, .edges, .txt
    
    output:
    ----------
    matrix: array
        The numpy array containing the adjacency matrix
    '''
    aa=np.loadtxt(fname,dtype='int')
    if bipartite==True:
        nodes_cols = int(max(aa[j,1] for j in range(aa.shape[0]))+1)
        nodes_rows= int(max(aa[j,0] for j in range(aa.shape[0]))+1)
        matrix=np.zeros((nodes_rows,nodes_cols),dtype='int')
        for j in range(aa.shape[0]):
            matrix[aa[j,0],aa[j,1]] = 1
    else:
        nodes = int(max(max(aa[j,0], aa[j,1]) for j in range(aa.shape[0]))+1)
        matrix=np.zeros((nodes,nodes),dtype='int')
        for j in range(aa.shape[0]):
            matrix[aa[j,0],aa[j,1]] = 1
        np.fill_diagonal(matrix, 0)
        matrix=np.triu(matrix,k=1)+(np.triu(matrix,k=1)).T
    return matrix