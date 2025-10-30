#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 20 16:16:57 2017

@author: maria
"""
#%%
import numpy as np
#%%
#%%
def NODF(M):
    '''
    function to calculate the nestedness by overlap and decreasing fill (NODF).
    Metric developed by Almeida-Neto et al., 2008.

    Inputs:
    ----------
        M: array
            An matrix to which I want to calculate the NODF
    
    output:
    ----------
    NODF: number
        The NODF score for the whole matrix
    '''
    rw,cl=M.shape
    colN=np.zeros((cl,cl))
    rowN=np.zeros((rw,rw))
    
    #Find NODF column score
    for i in range(cl-1): # at a left position with respect to column j
      	for j in range(i+1,cl):
              #if (i!=j):
            if (np.sum(M[:,i])>np.sum(M[:,j]))&(np.sum(M[:,j])>0): # DF =! to zero, then NP =! to zero
                colN[i,j]=(M[:,i]*M[:,j]).sum()/(np.sum(M[:,j]))
    
#    NODF_COL = (2*np.sum(colN)/(cl*(cl-1)))*100
    
    #Find NODF row score
    for i in range(rw-1): #at an upper position with respect to row j
        for j in range(i+1,rw):
            #if (i!=j):
            if (np.sum(M[i,:])>np.sum(M[j,:]))&(np.sum(M[j,:])>0): # DF =! to zero, then NP =! to zero
                rowN[i,j]=(M[i,:]*M[j,:]).sum()/(np.sum(M[j,:]))
    
#    NODF_ROW = (2*np.sum(rowN)/(rw*(rw-1)))*100
    
    #Find NODF
    NODF=(2*(np.sum(rowN)+np.sum(colN))/(cl*(cl-1) + rw*(rw-1) ))
    return  NODF

#%%
def spectral_radius(M):    
    """
    Spectral Radius described in Staniczenko et al., 2013
    INPUT: 
        M: array
            the bipartite biadjacency matrix
    OUTPUT: 
    - spectral radius.
  
    """
    #build the adajacency matrix
    r,cl=M.shape
    L=M.sum()
    ntotal=r+cl
    theta_ik=np.zeros((ntotal,ntotal))
    theta_ik[0:r,r::]=M;
    theta_ik[r::,0:r]=M.T;  
    max_eig = max(np.abs(np.linalg.eig(theta_ik)[0].real))
    return max_eig/np.sqrt(L)
#%%
# def glob_nestnulb(M):
#     '''
#     function to calculate the nestedness fitness N, a modified version
#     of Ñ, corrected by a null model
#     Metric developed by ASR et al, PRE 2018.

#     Inputs:n
#     ----------
#         M: array
#             A matrix to which I want to calculate the N
    
#     output:
#     ----------
#     N: number
#         The N score for the whole matrix
#     '''
#     rw,cl=M.shape
#     colN=np.zeros((cl,cl))
#     rowN=np.zeros((rw,rw))
#     cols_degr = M.sum(axis=0) # degree of the cols nodes
#     rows_degr = M.sum(axis=1) # degree of the rows nodes
    
#     #Find N col score
#     for i in range(cl): # at a left position with respect to column j
#       	for j in range(cl):
#               if M[i,j]==1:
#                   if (cols_degr[i]>=cols_degr[j]) & (cols_degr[j]>0): # heaviside
#                       if (cols_degr[i]==cols_degr[j]):
#                           colN[i,j]=(np.sum((M[:,i]*M[:,j]),dtype=float)-((cols_degr[i]*cols_degr[j])/rw))/(2*cols_degr[j]) #paired overlap
#                       else:
#                           colN[i,j]=(np.sum((M[:,i]*M[:,j]),dtype=float)-((cols_degr[i]*cols_degr[j])/rw))/cols_degr[j]
        
#     N_COL = (np.sum(colN,dtype=float)/(cl-1))
    
#     for i in range(rw): # at an upper position with respect to row j
#         for j in range(rw):
#             if M[i,j]==1:
#                 if (rows_degr[i]>=rows_degr[j]) & (rows_degr[j]>0): # Heaviside
#                     if (rows_degr[i]==rows_degr[j]):
#                         rowN[i,j]=(np.sum((M[i,:]*M[j,:]),dtype=float)-((rows_degr[i]*rows_degr[j])/cl))/(2*rows_degr[j]) #paired overlap
#                     else:
#                         rowN[i,j]=(np.sum((M[i,:]*M[j,:]),dtype=float)-((rows_degr[i]*rows_degr[j])/cl))/rows_degr[j] #paired overlap

#     N_ROW = (np.sum(rowN,dtype=float)/(rw-1))
    
#     #Find N
#     N=(N_COL+N_ROW)*(2./(rw+cl))
#     return N
#%%

def mathcalN_bipartite(M, equal_weight=1.0):
    """
    Global nestedness fitness N (null-model corrected) for a bipartite adjacency matrix.

    Parameters
    ----------
        M : array-like (Nr x Nc), binary 0/1
        equal_weight : float in [0,1]
            Contribution weight for equal-degree pairs.

    Returns
    -------
        N : float
            The N score for the whole matrix
    """
    A = np.asarray(M, dtype=float)
    Nr, Nc = A.shape
    if Nr < 2 and Nc < 2:
        return 0.0

    # degrees
    k_r = A.sum(axis=1)        # rows    
    k_c = A.sum(axis=0)        # cols
    
    # overlaps and null-model expectations
    O_r = A @ A.T                              # (Nr x Nr), common columns
    E_r = np.outer(k_r, k_r) / (Nc)
    
    O_c = A.T @ A                              # (Nc x Nc), common rows
    E_c = np.outer(k_c, k_c) / (Nr)
    
    # ===== weights (relaxed Heaviside), exclude self-pairs =====
    # rows
    kr_i = k_r[:, None]
    kr_j = k_r[None, :]
    denom_r = kr_j
    w_r = (kr_i > kr_j).astype(float)
    w_r += equal_weight * ((kr_i == kr_j) & (kr_i > 0))   # equals, but only if degree>0
    np.fill_diagonal(w_r, 0.0)                            # exclude i==j
    w_r = np.triu(w_r)
    w_r *= (denom_r > 0)                                  # no contribution if denom==0
    
    # cols
    kc_i = k_c[:, None]
    kc_j = k_c[None, :]
    denom_c = kc_j
    w_c = (kc_i > kc_j).astype(float)
    w_c += equal_weight * ((kc_i == kc_j) & (kc_i > 0))
    np.fill_diagonal(w_c, 0.0)
    w_c = np.triu(w_c)
    w_c *= (denom_c > 0)
    
    # ===== safe division: divide only where denom>0 to avoid inf; elsewhere 0 =====
    R = np.zeros_like(O_r, dtype=float)
    np.divide(O_r - E_r, denom_r, out=R, where=(denom_r > 0))
    R *= w_r

    C = np.zeros_like(O_c, dtype=float)
    np.divide(O_c - E_c, denom_c, out=C, where=(denom_c > 0))
    C *= w_c

    N_row = R.sum() / (Nr - 1)
    N_col = C.sum() / (Nc - 1)
    return float((2.0 / (Nr + Nc)) * (N_row + N_col))


def mathcalN_unipartite(M):
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
        N: float
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

