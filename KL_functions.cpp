//
//  EO_functions.cpp
//  
//
//  Created by María Palazzi Nieves on 13/11/18.
//
//mex -largeArrayDims IBN_bi.cpp

//#include <mex.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <istream>
#include <string>
#include <sstream>
#include <fstream>
#include <cmath>
#include <random>
#include <valarray>
#include <cstdio>
#include <numeric>
//#include "KL.hpp"
#include <cfloat>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/stl_bind.h>
namespace py = pybind11;
#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) < (Y)) ? (Y) : (X))
#define REMOVE true
#define ADD false
#define IBN true
#define MOD false

using namespace std;

vector<int> indexes_sorted_vector(vector<double> & x){
    
    std::vector<int> vs(x.size());
    std::iota(vs.begin(), vs.end(), 0);
    auto comparator = [&x](int a, int b){ return x[a] < x[b]; };
    std::sort(vs.begin(), vs.end(), comparator);
    return vs;
}

void lambdas_inblock(
                     vector<vector<int> > & input_matrix,
                     vector<int> & k_cols,
                     vector<int> & k_rows,
                     vector<int> label_cols,
                     vector<int> label_rows,
                     vector<double>& lambda_cols,
                     vector<double>& lambda_rows,
                     int max_number_blocks){
    
    int N_rows=input_matrix.size();
    int N_cols=input_matrix[0].size();
    lambda_cols= vector<double>(N_cols,0);
    lambda_rows = vector<double>(N_rows,0);
    vector<double> lambda_cols_aux = vector<double>(N_cols,0);
    vector<double> lambda_rows_aux = vector<double>(N_rows,0);
    
    vector<vector<int> > list_of_blocks_cols(max_number_blocks+1,vector<int>());
    vector<vector<int> > list_of_blocks_rows(max_number_blocks+1,vector<int>());
    
    for (size_t j = 0; j < label_cols.size(); ++j){
        list_of_blocks_cols[label_cols[j]].push_back(j);
    }
    
    //rows
    for (size_t j = 0; j < label_rows.size(); ++j){
        list_of_blocks_rows[label_rows[j]].push_back(j);
    }
    
    //computing the column nodes contribution (pair overlap)
    for (size_t l=0; l<list_of_blocks_cols.size(); l++){
        int size_blocks_cols=(int)list_of_blocks_cols[l].size();
        
        
        for (int i=0; i<size_blocks_cols; i++){
            int ii = list_of_blocks_cols[l][i];
            
            
            for (int j=0; j<size_blocks_cols; j++){
                
                int jj = list_of_blocks_cols[l][j];
                double PO_col_i=0;
                double normalization = ((double)(k_cols[jj]))*((double)(size_blocks_cols-1.0));
                
                if ((k_cols[ii]>k_cols[jj])){
                    
                    for (size_t k=0; k<list_of_blocks_rows[l].size(); k++){
                        
                        int kk = list_of_blocks_rows[l][k];
                        
                        if ((input_matrix[kk][ii]*input_matrix[kk][jj])==1){
                            PO_col_i++;
                        }
                    }
                    
                    double null_model = ((double)(k_cols[ii]*k_cols[jj]))/((double)(N_rows));
                    
                    PO_col_i=(PO_col_i-null_model)/normalization;
                    if(isnan(PO_col_i)){
                        PO_col_i = 0;
                    }
                    //                    sprintf(outputString,"fprintf('overlap col node %i,%i = %f  \\n');",ii,jj,PO_col_i-null_model);
                    //                    mexEvalString(outputString);
                    
                    lambda_cols[ii]+=PO_col_i;
                    
                }
            }
        }
    }
    
    
    //computing the rows pair overlap
    for (size_t l=0; l<list_of_blocks_rows.size(); l++){
        
        int size_blocks_rows=(int)list_of_blocks_rows[l].size();
        for (int i=0; i<(size_blocks_rows); i++){
            
            int ii = list_of_blocks_rows[l][i];
            for (int j=0; j<size_blocks_rows; j++){
                
                int jj =list_of_blocks_rows[l][j];
                double PO_row_i=0;
                double normalization = ((double)(k_rows[jj]))*((double)(size_blocks_rows-1.0));
                
                if ((k_rows[ii]>k_rows[jj])){
                    
                    for (size_t k=0; k<list_of_blocks_cols[l].size(); k++){
                        int kk = list_of_blocks_cols[l][k];
                        
                        if (input_matrix[ii][kk] == 1 && input_matrix[jj][kk]==1){
                            PO_row_i = PO_row_i + 1;
                        }
                    }
                    double null_model = ((double)(k_rows[ii]*k_rows[jj]))/((double)(N_cols));
                    
                    PO_row_i=(PO_row_i-null_model)/normalization;
                    if(isnan(PO_row_i)){
                        PO_row_i = 0;
                    }
                    lambda_rows[ii]+=PO_row_i;
                }
            }
        }
    }
    
}

void lambdas_inblock_change_node_partition_row(
                                               int nodeId,
                                               int paritionId,
                                               vector<vector<int> > & input_matrix,
                                               vector<int> & k_cols,
                                               vector<int> & k_rows,
                                               vector<int> label_cols,
                                               vector<int> label_rows,
                                               vector<double> & lambda_cols,
                                               vector<double> & lambda_rows,
                                               int is_removal,
                                               int max_number_blocks){
    
    int N_cols=input_matrix[0].size();
    
    vector<vector<int> > list_of_blocks_cols(max_number_blocks+1,vector<int>());
    vector<vector<int> > list_of_blocks_rows(max_number_blocks+1,vector<int>());
    
    //creating a vector that separate the nodes according to the block the belong to columns
    for (size_t j = 0; j < label_cols.size(); ++j){
        list_of_blocks_cols[label_cols[j]].push_back(j);
    }
    
    //rows
    for (size_t j = 0; j < label_rows.size(); ++j){
        list_of_blocks_rows[label_rows[j]].push_back(j);
    }
    
    //computing the column nodes contribution (pair overlap)
    int l = paritionId;
    int size_blocks_rows= (double)list_of_blocks_rows[l].size();
    int ii = nodeId;
    if(is_removal){
        lambda_rows[ii] = 0;
    }else{
        //recompute the contribution of the node to the other community
        
        for (int j=0; j<size_blocks_rows; j++){
            int jj =list_of_blocks_rows[l][j];
            double PO_row_i=0;
            double normalization = ((double)(k_rows[jj]))*((double)((size_blocks_rows-1)+1));
            
            if ( (k_rows[ii]>k_rows[jj])){
                
                for (size_t k=0; k<list_of_blocks_cols[l].size(); k++){
                    int kk = list_of_blocks_cols[l][k];
                    if ((input_matrix[ii][kk]*input_matrix[jj][kk])==1){
                        PO_row_i++;
                    }
                }
                
                double null_model = ((double)(k_rows[ii]*k_rows[jj]))/((double)(N_cols));
                PO_row_i=(PO_row_i-null_model)/normalization;
                if(isnan(PO_row_i)){
                    PO_row_i = 0;
                }
                lambda_rows[ii]+=PO_row_i;
            }
        }
    }
    
    //update contribution for the same dimension where we update the node
    for (int j=0; j<size_blocks_rows; j++){
        int jj =list_of_blocks_rows[l][j];
        double PO_row_i=0;
        double normalization_old = ((double)((size_blocks_rows-1)));
        double normalization_new;
        if(is_removal){
            normalization_new = ((double)((size_blocks_rows-1)-1));
        }else{
            normalization_new = ((double)((size_blocks_rows-1)+1));
        }
        
        if ((k_rows[ii]<k_rows[jj])){
            //remove old constant
            lambda_rows[jj] = normalization_old*lambda_rows[jj];
            for (size_t k=0; k<list_of_blocks_cols[l].size(); k++){
                int kk = list_of_blocks_cols[l][k];
                if ((input_matrix[ii][kk]*input_matrix[jj][kk])==1){
                    PO_row_i++;
                }
            }
            
            double null_model = ((double)(k_rows[ii]*k_rows[jj]))/((double)(N_cols));
            PO_row_i=(PO_row_i-null_model)/((double)(k_rows[ii]));
            
            if(isnan(PO_row_i)){
                PO_row_i = 0;
            }
            
            if(is_removal){
                lambda_rows[jj] = (lambda_rows[jj]-PO_row_i)/normalization_new;
            }else{
                lambda_rows[jj] = (lambda_rows[jj]+PO_row_i)/normalization_new;
            }
            
            if(isnan(lambda_rows[jj])){
                lambda_rows[jj] = 0;
            }
        }
        
        if ( (k_rows[ii]>=k_rows[jj])){
            //remove old constant
            lambda_rows[jj] = normalization_old*lambda_rows[jj];
            lambda_rows[jj] = lambda_rows[jj]/normalization_new;
            
            if(isnan(lambda_rows[jj])){
                lambda_rows[jj] = 0;
            }
        }
    }
    
    //update contribution for the other dimension where we update the node
    l = paritionId;
    double size_blocks_cols=(double)list_of_blocks_cols[l].size();
    for (int i=0; i<(size_blocks_cols); i++){
        
        int ii = list_of_blocks_cols[l][i];
        for (int j=0; j<size_blocks_cols; j++){
            
            int jj =list_of_blocks_cols[l][j];
            
            double normalization = ((double)(k_cols[jj]))*((double)(size_blocks_cols-1.0));
            
            if ((k_cols[ii]>k_cols[jj]) ){
                double PO_col_i=0;
                
                int kk = nodeId;
                if ((input_matrix[kk][ii]*input_matrix[kk][jj])==1){
                    PO_col_i++;
                }
                
                double delta_value = PO_col_i/normalization;
                
                if(isnan(delta_value)){
                    delta_value = 0;
                }
                
                if(is_removal){
                    lambda_cols[ii] = lambda_cols[ii] - delta_value;
                }else{
                    lambda_cols[ii] = lambda_cols[ii] + delta_value;
                }
                
            }
        }
    }
}


void lambdas_inblock_change_node_partition_col(
                                               int nodeId,
                                               int paritionId,
                                               vector<vector<int> > & input_matrix,
                                               vector<int> & k_cols,
                                               vector<int> & k_rows,
                                               vector<int> label_cols,
                                               vector<int> label_rows,
                                               vector<double>& lambda_cols,
                                               vector<double>& lambda_rows,
                                               int is_removal,
                                               int max_number_blocks){
    
    int N_rows=input_matrix.size();
    
    vector<vector<int> > list_of_blocks_cols(max_number_blocks+1,vector<int>());
    vector<vector<int> > list_of_blocks_rows(max_number_blocks+1,vector<int>());
    
    //creating a vector that separate the nodes according to the block the belong to columns
    for (size_t j = 0; j < label_cols.size(); ++j){
        list_of_blocks_cols[label_cols[j]].push_back(j);
    }
    
    //rows
    for (size_t j = 0; j < label_rows.size(); ++j){
        list_of_blocks_rows[label_rows[j]].push_back(j);
    }
    
    //computing the column nodes contribution (pair overlap)
    int l = paritionId;
    int size_blocks_cols=(double)list_of_blocks_cols[l].size();
    int ii = nodeId;
    
    if(is_removal){
        lambda_cols[ii] = 0;
    }else{
        //recompute the contribution of the node to the other community
        for (int j=0; j<size_blocks_cols; j++){
            
            int jj =list_of_blocks_cols[l][j];
            double PO_col_i=0;
            double normalization = ((double)(k_cols[jj]))*((double)((size_blocks_cols-1)+1));
            
            if ((k_cols[ii]>k_cols[jj])){
                
                for (size_t k=0; k<list_of_blocks_rows[l].size(); k++){
                    int kk = list_of_blocks_rows[l][k];
                    if ((input_matrix[kk][ii]*input_matrix[kk][jj])==1){
                        PO_col_i++;
                    }
                }
                
                double null_model = ((double)(k_cols[ii]*k_cols[jj]))/((double)(N_rows));
                PO_col_i=(PO_col_i-null_model)/normalization;
                if(isnan(PO_col_i)){
                    PO_col_i = 0;
                }
                
                lambda_cols[ii]+=PO_col_i;
            }
        }
    }
    
    //update contribution for the same dimension where we update the node
    for (int j=0; j<size_blocks_cols; j++){
        int jj =list_of_blocks_cols[l][j];
        double PO_col_i=0;
        double normalization_old = ((double)((size_blocks_cols-1)));
        double normalization_new;
        if(is_removal){
            normalization_new = ((double)((size_blocks_cols-1)-1));
        }else{
            normalization_new = ((double)((size_blocks_cols-1)+1));
        }
        
        if ((k_cols[ii]<k_cols[jj])){
            //remove old constant
            lambda_cols[jj] = normalization_old*lambda_cols[jj];
            
            for (size_t k=0; k<list_of_blocks_rows[l].size(); k++){
                int kk = list_of_blocks_rows[l][k];
                if ((input_matrix[kk][ii]*input_matrix[kk][jj])==1){
                    PO_col_i++;
                }
            }
            
            double null_model = ((double)(k_cols[ii]*k_cols[jj]))/((double)(N_rows));
            PO_col_i=(PO_col_i-null_model)/((double)(k_cols[ii]));
            
            if(isnan(PO_col_i)){
                PO_col_i = 0;
            }
            
            if(is_removal){
                lambda_cols[jj] = (lambda_cols[jj]-PO_col_i)/normalization_new;
            }else{
                lambda_cols[jj] = (lambda_cols[jj]+PO_col_i)/normalization_new;
            }
            
            if(isnan(lambda_cols[jj])){
                lambda_cols[jj] = 0;
            }
        }
        
        if ( (k_cols[ii]>=k_cols[jj]) ){
            //remove old constant
            lambda_cols[jj] = normalization_old*lambda_cols[jj];
            lambda_cols[jj] = lambda_cols[jj]/normalization_new;
            
            if(isnan(lambda_cols[jj])){
                lambda_cols[jj] = 0;
            }
        }
    }
    
    
    //update contribution for the other dimension where we update the node
    l = paritionId;
    double size_blocks_rows=(double)list_of_blocks_rows[l].size();
    for (int i=0; i<(size_blocks_rows); i++){
        
        int ii = list_of_blocks_rows[l][i];
        for (int j=0; j<size_blocks_rows; j++){
            
            int jj =list_of_blocks_rows[l][j];
            
            double normalization = ((double)(k_rows[jj]))*((double)(size_blocks_rows-1.0));
            
            if ((k_rows[ii]>k_rows[jj]) ){
                double PO_row_i=0;
                
                int kk = nodeId;
                if ((input_matrix[ii][kk]*input_matrix[jj][kk])==1){
                    PO_row_i++;
                }
                
                double delta_value = PO_row_i/normalization;
                
                if(isnan(delta_value)){
                    delta_value = 0;
                }
                if(is_removal){
                    lambda_rows[ii] = lambda_rows[ii] - delta_value;
                }else{
                    lambda_rows[ii] = lambda_rows[ii] + delta_value;
                }
                
            }
        }
    }
}


void lambdas_modularity(
                        vector<vector<int> > & input_matrix,
                        vector<int> & k_cols,
                        vector<int> & k_rows,
                        vector<int> label_cols,
                        vector<int> label_rows,
                        vector<double>& lambda_cols,
                        vector<double>& lambda_rows,
                        int max_number_blocks){
    
    int N_rows=input_matrix.size();
    int N_cols=input_matrix[0].size();
    lambda_rows = vector<double>(N_rows,0);
    lambda_cols=vector<double>(N_cols,0);
    vector<double> kappa_rows(N_rows,0);
    vector<double> kappa_cols(N_cols,0);
    
    //obtaining current number of blocks
    vector<double> links_blocks_rows(max_number_blocks,0);
    vector<double> links_blocks_cols(max_number_blocks,0);
    double total_links=0.;
    
    //getting the total links
    total_links = accumulate(k_rows.begin(), k_rows.end(), 0.0);
    
    // getting the kappa's
    for (unsigned int i = 0; i < label_rows.size(); ++i) {
        for (unsigned int j = 0; j < label_cols.size(); ++j) {
            if  ((input_matrix[i][j]==1) && (label_rows[i]==label_cols[j])){
                kappa_rows[i]+=1./(double)k_rows[i];
                kappa_cols[j]+=1./(double)k_cols[j];
            }
        }
    }
    
    // getting the a_r(i)'s
    for (unsigned int i = 0; i < label_rows.size(); ++i) {
        links_blocks_rows[label_rows[i]]+=k_rows[i]/total_links;
    }
    
    for (unsigned int i = 0; i < label_cols.size(); ++i) {
        links_blocks_cols[label_cols[i]]+=k_cols[i]/total_links;
    }
    //getting the lambda fitness
    for (unsigned int i = 0; i < label_rows.size(); ++i) {
        lambda_rows[i]=kappa_rows[i] - links_blocks_rows[label_rows[i]];
    }
    for (unsigned int i = 0; i < label_cols.size(); ++i) {
        lambda_cols[i]=kappa_cols[i] - links_blocks_cols[label_cols[i]];
    }
}

void lambdas_modularity_change_node_row(
                                        unsigned int nodeId,
                                        int paritionId,
                                        vector<vector<int> > & input_matrix,
                                        vector<int> & k_cols,
                                        vector<int> & k_rows,
                                        vector<int> label_cols,
                                        vector<int> label_rows,
                                        vector<double> & lambda_cols,
                                        vector<double> & lambda_rows,
                                        int operation,
                                        int max_number_blocks){
    
    int N_rows=input_matrix.size();
    int N_cols=input_matrix[0].size();
    vector<double> kappa_rows_diff(N_rows,0);
    vector<double> kappa_cols_diff(N_cols,0);
    
    //obtaining current number of blocks
    vector<double> links_blocks_rows_diff(max_number_blocks,0);
    vector<double> links_blocks_cols_diff(max_number_blocks,0);
    vector<double> links_blocks_rows(max_number_blocks,0);
    vector<double> links_blocks_cols(max_number_blocks,0);
    double total_links=0;
    
    total_links = accumulate(k_rows.begin(), k_rows.end(), 0.0);
    
    unsigned int i = nodeId;
    for (unsigned int j = 0; j < label_cols.size(); ++j) {
        if  ((input_matrix[i][j]==1) && (paritionId==label_cols[j])){
            kappa_rows_diff[i]+=1./(double)k_rows[i];
            kappa_cols_diff[j]+=1./(double)k_cols[j];
        }
    }
    
    // getting the a_r(i)'s
    if(operation == ADD){
        for (unsigned int i = 0; i < label_rows.size(); ++i) {
            if(paritionId == label_rows[i] && i != nodeId){
                links_blocks_rows[label_rows[i]]+=k_rows[i]/total_links;
            }
        }
    }
    links_blocks_rows_diff[paritionId] = k_rows[nodeId]/total_links;
    
    for (unsigned int i = 0; i < label_rows.size(); i++) {
        if(i==nodeId){
            if(operation == ADD){
                lambda_rows[i]= lambda_rows[i]
                + kappa_rows_diff[i] - links_blocks_rows[paritionId] - links_blocks_rows_diff[paritionId];
            }else{
                lambda_rows[i] = 0;//lambda_rows[i] - kappa_rows_diff[i] + links_blocks_rows_diff[paritionId];
            }
        }else{
            if (label_rows[i] == paritionId){
                if(operation == ADD){
                    lambda_rows[i]= lambda_rows[i] - links_blocks_rows_diff[label_rows[i]];
                }else{
                    lambda_rows[i]= lambda_rows[i] + links_blocks_rows_diff[label_rows[i]];
                }
            }else{
                //lambda_rows[i] = lambda_rows[i] - kappa_rows_diff[i];
            }
        }
    }
    
    for (unsigned int j = 0; j < label_cols.size(); j++) {
        if(operation == ADD){
            lambda_cols[j] = lambda_cols[j] + kappa_cols_diff[j];
        }else{
            lambda_cols[j] = lambda_cols[j] - kappa_cols_diff[j];
        }
    }
}

void lambdas_modularity_change_node_col(
                                        unsigned int nodeId,
                                        int paritionId,
                                        vector<vector<int> > & input_matrix,
                                        vector<int> & k_cols,
                                        vector<int> & k_rows,
                                        vector<int> label_cols,
                                        vector<int> label_rows,
                                        vector<double>& lambda_cols,
                                        vector<double>& lambda_rows,
                                        int operation,
                                        int max_number_blocks){
    
    int N_rows=input_matrix.size();
    int N_cols=input_matrix[0].size();
    vector<double> kappa_rows_diff(N_rows,0);
    vector<double> kappa_cols_diff(N_cols,0);
    //vector<int> k_rows(N_rows,0);
    //vector<int> k_cols(N_cols,0);
    
    //obtaining current number of blocks
    vector<double> links_blocks_rows_diff(max_number_blocks,0);
    vector<double> links_blocks_cols_diff(max_number_blocks,0);
    vector<double> links_blocks_rows(max_number_blocks,0);
    vector<double> links_blocks_cols(max_number_blocks,0);
    double total_links=0.;
    
    total_links = accumulate(k_rows.begin(), k_rows.end(), 0.0);
    
    for (unsigned int i = 0; i < label_rows.size(); ++i) {
        unsigned int j = nodeId;
        if  ((input_matrix[i][j]==1) && (label_rows[i]==paritionId)){
            kappa_rows_diff[i]+=1./(double)k_rows[i];
            kappa_cols_diff[j]+=1./(double)k_cols[j];
        }
    }
    
    // getting the a_r(i)'s
    if(operation == ADD){
        for (unsigned int i = 0; i < label_cols.size(); ++i) {
            if(paritionId == label_cols[i] && i != nodeId){
                links_blocks_cols[label_cols[i]]+=k_cols[i]/total_links;
            }
        }
    }
    links_blocks_cols_diff[paritionId] = k_cols[nodeId]/total_links;
    
    //getting the lambda fitness
    for (unsigned int i = 0; i < label_cols.size(); i++) {
        if(i==nodeId){
            if(operation == ADD){
                lambda_cols[i]= lambda_cols[i]
                + kappa_cols_diff[i] - links_blocks_cols[paritionId] - links_blocks_cols_diff[paritionId];
            }else{
                lambda_cols[i] = 0;//lambda_rows[i] - kappa_rows_diff[i] + links_blocks_rows_diff[paritionId];
            }
        }else{
            if (label_cols[i] == paritionId){
                if(operation == ADD){
                    lambda_cols[i]= lambda_cols[i] - links_blocks_cols_diff[label_cols[i]];
                }else{
                    lambda_cols[i]= lambda_cols[i] + links_blocks_cols_diff[label_cols[i]];
                }
            }else{
                //lambda_rows[i] = lambda_rows[i] - kappa_rows_diff[i];
            }
        }
    }
    
    for (unsigned int j = 0; j < label_rows.size(); j++) {
        if(operation == ADD){
            lambda_rows[j] = lambda_rows[j] + kappa_rows_diff[j];
        }else{
            lambda_rows[j] = lambda_rows[j] - kappa_rows_diff[j];
        }
    }
}

double calculate_inblock_nestedness(
                                    vector<vector<int> > & input_matrix,
                                    vector<double> lambda_cols,
                                    vector<double> lambda_rows){
    
    double I;
    int N_rows=input_matrix.size();
    int N_cols=input_matrix[0].size();
    double i_col=0;
    double i_row=0;
    
    i_col = accumulate(lambda_cols.begin(), lambda_cols.end(), 0.0);
    
    i_row = accumulate(lambda_rows.begin(), lambda_rows.end(), 0.0);
    
    I=(2.0/((double)N_cols+(double)N_rows))*(i_col+i_row);
    return I;
}

double calculate_modularity(
                            vector<vector<int> > & input_matrix,
                            vector<int> & k_cols,
                            vector<int> & k_rows,
                            vector<double> lambda_cols,
                            vector<double> lambda_rows){
    
    double Q;
    double q_rows=0;
    double q_cols=0;
    int N_rows=input_matrix.size();
    int N_cols=input_matrix[0].size();
    double total_links=0;
    
    total_links = accumulate(k_rows.begin(), k_rows.end(), 0.0);
    
    for (int i=0; i<N_cols; i++){
        q_cols+=k_cols[i]*lambda_cols[i];
    }
    for (int j=0; j<N_rows; j++){
        q_rows+=k_rows[j]*lambda_rows[j];
    }
    Q=(q_cols+q_rows)/(2*total_links);
    return Q;
    
}

void lambda_i(
              vector<vector<int> > & input_matrix,
              vector<int> & k_cols,
              vector<int> & k_rows,
              vector<int> label_cols,
              vector<int> label_rows,
              vector<double> & lambda_cols,
              vector<double> & lambda_rows,
              int max_number_blocks,
              bool ibn){
    
    if (ibn==true){
        /* Function that calculate the fitness contributions of nodes (lambda_i) for the in-block nestedness*/
        lambdas_inblock(
                        input_matrix,
                        k_cols,
                        k_rows,
                        label_cols,
                        label_rows,
                        lambda_cols,
                        lambda_rows,
                        max_number_blocks);
    }else {
        /* Function that calculate the fitness contributions of nodes (lambda_i) for the modularity*/
        lambdas_modularity(
                           input_matrix,
                           k_cols,
                           k_rows,
                           label_cols,
                           label_rows,
                           lambda_cols,
                           lambda_rows,
                           max_number_blocks);
    }
}

vector<vector<double> > call_lambda_i(
                                      vector<vector<int> > & input_matrix,
                                      vector<int> & k_cols,
                                      vector<int> & k_rows,
                                      vector<int> label_cols,
                                      vector<int> label_rows,
                                      int max_number_blocks,
                                      bool ibn){
    vector<vector<double> > out;
    vector<double> lambda_cols(label_cols.size(),0);
    vector<double> lambda_rows(label_rows.size(),0);
    lambda_i(
             input_matrix,
             k_cols,
             k_rows,
             label_cols,
             label_rows,
             lambda_cols,
             lambda_rows,
             max_number_blocks,
             ibn);
    out.push_back(lambda_cols);
    out.push_back(lambda_rows);
    return out;
}

double calculate_Fitness(
                         vector<vector<int> > & input_matrix,
                         vector<int>  k_cols,
                         vector<int>  k_rows,
                         vector<double> lambda_cols,
                         vector<double> lambda_rows,
                         bool ibn){
    
    /*  Calculate the EO metric of the whole network */
    double metric;
    if (ibn==true) {
        // in block nestedness
        metric=calculate_inblock_nestedness(input_matrix, lambda_cols,lambda_rows);
    }else{ // modularity
        metric=calculate_modularity(
                                    input_matrix,
                                    k_cols,
                                    k_rows,
                                    lambda_cols,
                                    lambda_rows);
    }
    return metric;
}


/*Obtain the node that is going to be moved from its block, based on
 probability distribution. tau-EO like in Duch et al 2005.*/
int low_fitness_node(
                     vector<double> lambda_cols,
                     vector<double> lambda_rows) {
    
    int low_node;
    int N=lambda_rows.size()+lambda_cols.size();
    double tau = 1.+1./log(N); //tau-exponent
    vector<double> probabilities(N,0);
    double p_norm=0;
    vector<double> lambda_stacked=lambda_rows;
    vector<int> lambda_sorted(N);
    
    // concatenating lambda_row and lambda_col into one vector
    lambda_stacked.insert( lambda_stacked.end(), lambda_cols.begin(), lambda_cols.end() );
    
    // generate distribution of probabilities
    for (int i=0; i<N ; i++){
        probabilities[i]= pow(i+1,-tau);
        p_norm+=probabilities[i];
    }
    for (int j=0; j<N; j++){
        probabilities[j]=probabilities[j]/p_norm;
    }
    
    discrete_distribution<int> distribution(probabilities.begin(), probabilities.end());
    random_device rd;
    mt19937 gen(rd());
    
    // sorting the lambda_stacked vector
    lambda_sorted = indexes_sorted_vector(lambda_stacked);
    low_node=lambda_sorted[distribution(gen)];
    
    return low_node;
}

void update_partitions_labels(
                              vector<int> & total_in_label_cols,
                              vector<int> & total_in_label_rows,
                              vector<int> label_partition_cols,
                              vector<int> label_partition_rows,
                              int current_block){
    
    int j=0;
    //updating the label vector with the new partitions
    for (unsigned int i=0; i<total_in_label_cols.size(); i++){
        if (total_in_label_cols[i]==current_block){
            total_in_label_cols[i]=label_partition_cols[j];
            j++;
        }
    }
    int k=0;
    for (unsigned int l=0; l<total_in_label_rows.size(); l++){
        if (total_in_label_rows[l]==current_block){
            total_in_label_rows[l]=label_partition_rows[k];
            k++;
        }
    }
}

void load_matrix(istream* is,vector< vector<int> >* matrix,\
                 const string& delim = " \t"){
    string      line;
    string      strnum;
    
    // clear first
    matrix->clear();
    
    // parse line by line
    while (getline(*is, line))
    {
        matrix->push_back(vector<int>());
        
        for (string::const_iterator i = line.begin(); i != line.end(); ++ i)
        {
            // If i is not a delim, then append it to strnum
            if (delim.find(*i) == string::npos)
            {
                strnum += *i;
                if (i + 1 != line.end()) // If it's the last char, do not continue
                    continue;
            }
            
            // if strnum is still empty, it means the previous char is also a
            // delim (several delims appear together). Ignore this char.
            if (strnum.empty())
                continue;
            
            // If we reach here, we got a number. Convert it to double.
            int       number;
            
            istringstream(strnum) >> number;
            matrix->back().push_back(number);
            
            strnum.clear();
        }
    }
    //log_file("File read: pass");
}


struct data_out_KL {
    double QI_kl;
    vector<vector<int> > partitions_result;
};

struct KL_bestMovement_struct {
    int best_node_ind = -1;
    int best_comm_index = -1;
    double max_found = -DBL_MAX;
};

void find_Qmax(KL_bestMovement_struct & KL_bestMovement,
               vector<vector<double>> Q, int num_nodes, int num_part){
    
    KL_bestMovement.best_node_ind = -1;
    KL_bestMovement.best_comm_index = -1;
    KL_bestMovement.max_found = -DBL_MAX;
    
    for(int i=0; i<num_nodes; i++){
        for(int c=0; c<num_part; c++){
            if( Q[i][c] > KL_bestMovement.max_found ){
                KL_bestMovement.best_node_ind = i;
                KL_bestMovement.best_comm_index = c;
                KL_bestMovement.max_found = Q[i][c];
            }
        }
    }
}

KL_bestMovement_struct KL_IBN_optimization_testMovement_rows(
                                                             vector<vector<int>> & M,
                                                             vector<int> & k_rows,
                                                             vector<int> & k_cols,
                                                             vector<int> partition_rows,
                                                             vector<int> partition_cols,
                                                             vector<vector<double>> & Q_rows,
                                                             int max_number_blocks){
    
    int N_rows=M.size();
    
    vector<vector<double> > lambdas = call_lambda_i(
                                                    M,
                                                    k_cols,
                                                    k_rows,
                                                    partition_cols,
                                                    partition_rows,
                                                    max_number_blocks,
                                                    IBN);
    vector<double> lambdas_rows = lambdas[1];
    vector<double> lambdas_cols = lambdas[0];
    
    vector<int> partition_rows_temp = partition_rows;
    vector<int> partition_cols_temp = partition_cols;
    vector<double> lambdas_rows_temp;
    vector<double> lambdas_cols_temp;
    KL_bestMovement_struct KL_bestMovement;
    KL_bestMovement.best_node_ind = -1;
    KL_bestMovement.best_comm_index = -1;
    KL_bestMovement.max_found = -DBL_MAX;
    //compute increase in fitness for changing node
    for(int node_ind=0; node_ind < N_rows; node_ind++){
        for(int comm_index = 0; comm_index < max_number_blocks; comm_index++){
            
            if(partition_rows_temp[node_ind] != comm_index){
                //fast way
                partition_rows_temp = partition_rows;
                partition_cols_temp = partition_cols;
                lambdas_rows_temp = lambdas_rows;
                lambdas_cols_temp = lambdas_cols;
                
                lambdas_inblock_change_node_partition_row(
                                                          node_ind,
                                                          partition_rows_temp[node_ind],
                                                          M,
                                                          k_cols,
                                                          k_rows,
                                                          partition_cols_temp,
                                                          partition_rows_temp,
                                                          lambdas_cols_temp,
                                                          lambdas_rows_temp,
                                                          REMOVE,
                                                          max_number_blocks);
                
                lambdas_inblock_change_node_partition_row(
                                                          node_ind,
                                                          comm_index,
                                                          M,
                                                          k_cols,
                                                          k_rows,
                                                          partition_cols_temp,
                                                          partition_rows_temp,
                                                          lambdas_cols_temp,
                                                          lambdas_rows_temp,
                                                          ADD,
                                                          max_number_blocks);
                
                Q_rows[node_ind][comm_index] = calculate_Fitness(
                                                                 M,
                                                                 k_cols,
                                                                 k_rows,
                                                                 lambdas_cols_temp,
                                                                 lambdas_rows_temp,
                                                                 IBN);
                
            }else{
                Q_rows[node_ind][comm_index] = -DBL_MAX;
            }
            
            if(Q_rows[node_ind][comm_index]>KL_bestMovement.max_found){
                KL_bestMovement.best_node_ind = node_ind;
                KL_bestMovement.best_comm_index = comm_index;
                KL_bestMovement.max_found = Q_rows[node_ind][comm_index];
            }
        }
        //or new community???
        //fast way
        partition_rows_temp = partition_rows;
        partition_cols_temp = partition_cols;
        lambdas_rows_temp = lambdas_rows;
        lambdas_cols_temp = lambdas_cols;
        
        lambdas_inblock_change_node_partition_row(
                                                  node_ind,
                                                  partition_rows_temp[node_ind],
                                                  M,
                                                  k_cols,
                                                  k_rows,
                                                  partition_cols_temp,
                                                  partition_rows_temp,
                                                  lambdas_cols_temp,
                                                  lambdas_rows_temp,
                                                  REMOVE,
                                                  max_number_blocks);
        
        lambdas_inblock_change_node_partition_row(
                                                  node_ind,
                                                  max_number_blocks,
                                                  M,
                                                  k_cols,
                                                  k_rows,
                                                  partition_cols_temp,
                                                  partition_rows_temp,
                                                  lambdas_cols_temp,
                                                  lambdas_rows_temp,
                                                  ADD,
                                                  max_number_blocks+1);
        
        Q_rows[node_ind][max_number_blocks] = calculate_Fitness(
                                                                M,
                                                                k_cols,
                                                                k_rows,
                                                                lambdas_cols_temp,
                                                                lambdas_rows_temp,
                                                                IBN);
        
        if(Q_rows[node_ind][max_number_blocks]>KL_bestMovement.max_found){
            KL_bestMovement.best_node_ind = node_ind;
            KL_bestMovement.best_comm_index = max_number_blocks;
            KL_bestMovement.max_found = Q_rows[node_ind][max_number_blocks];
        }
    }
    
    return KL_bestMovement;
}

KL_bestMovement_struct KL_IBN_optimization_testMovement_cols(
                                                             vector<vector<int>> & M,
                                                             vector<int> & k_rows,
                                                             vector<int> & k_cols,
                                                             vector<int> partition_rows,
                                                             vector<int> partition_cols,
                                                             vector<vector<double>> & Q_cols,
                                                             int max_number_blocks){
    
    int N_cols=M[0].size();
    
    vector<vector<double> > lambdas = call_lambda_i(
                                                    M,
                                                    k_cols,
                                                    k_rows,
                                                    partition_cols,
                                                    partition_rows,
                                                    max_number_blocks,
                                                    IBN);
    vector<double> lambdas_rows = lambdas[1];
    vector<double> lambdas_cols = lambdas[0];
    
    vector<int> partition_rows_temp = partition_rows;
    vector<int> partition_cols_temp = partition_cols;
    vector<double> lambdas_rows_temp;
    vector<double> lambdas_cols_temp;
    KL_bestMovement_struct KL_bestMovement;
    KL_bestMovement.best_node_ind = -1;
    KL_bestMovement.best_comm_index = -1;
    KL_bestMovement.max_found = -DBL_MAX;
    //compute increase in fitness for changing node
    for(int node_ind=0; node_ind < N_cols; node_ind++){
        for(int comm_index = 0; comm_index < max_number_blocks; comm_index++){
            
            if(partition_cols_temp[node_ind] != comm_index){
                
                //fast way
                partition_rows_temp = partition_rows;
                partition_cols_temp = partition_cols;
                lambdas_rows_temp = lambdas_rows;
                lambdas_cols_temp = lambdas_cols;
                
                lambdas_inblock_change_node_partition_col(
                                                          node_ind,
                                                          partition_cols_temp[node_ind],
                                                          M,
                                                          k_cols,
                                                          k_rows,
                                                          partition_cols_temp,
                                                          partition_rows_temp,
                                                          lambdas_cols_temp,
                                                          lambdas_rows_temp,
                                                          REMOVE,
                                                          max_number_blocks);
                
                lambdas_inblock_change_node_partition_col(
                                                          node_ind,
                                                          comm_index,
                                                          M,
                                                          k_cols,
                                                          k_rows,
                                                          partition_cols_temp,
                                                          partition_rows_temp,
                                                          lambdas_cols_temp,
                                                          lambdas_rows_temp,
                                                          ADD,
                                                          max_number_blocks);
                
                Q_cols[node_ind][comm_index] = calculate_Fitness(
                                                                 M,
                                                                 k_cols,
                                                                 k_rows,
                                                                 lambdas_cols_temp,
                                                                 lambdas_rows_temp,
                                                                 IBN);
                
            }else{
                Q_cols[node_ind][comm_index] = -DBL_MAX;
            }
            
            if(Q_cols[node_ind][comm_index]>KL_bestMovement.max_found){
                KL_bestMovement.best_node_ind = node_ind;
                KL_bestMovement.best_comm_index = comm_index;
                KL_bestMovement.max_found = Q_cols[node_ind][comm_index];
            }
        }
        //or new community???
        //fast way
        partition_rows_temp = partition_rows;
        partition_cols_temp = partition_cols;
        lambdas_rows_temp = lambdas_rows;
        lambdas_cols_temp = lambdas_cols;
        
        lambdas_inblock_change_node_partition_col(
                                                  node_ind,
                                                  partition_cols_temp[node_ind],
                                                  M,
                                                  k_cols,
                                                  k_rows,
                                                  partition_cols_temp,
                                                  partition_rows_temp,
                                                  lambdas_cols_temp,
                                                  lambdas_rows_temp,
                                                  REMOVE,
                                                  max_number_blocks);
        
        lambdas_inblock_change_node_partition_col(
                                                  node_ind,
                                                  max_number_blocks,
                                                  M,
                                                  k_cols,
                                                  k_rows,
                                                  partition_cols_temp,
                                                  partition_rows_temp,
                                                  lambdas_cols_temp,
                                                  lambdas_rows_temp,
                                                  ADD,
                                                  max_number_blocks+1);
        
        Q_cols[node_ind][max_number_blocks] = calculate_Fitness(
                                                                M,
                                                                k_cols,
                                                                k_rows,
                                                                lambdas_cols_temp,
                                                                lambdas_rows_temp,
                                                                IBN);
        
        if(Q_cols[node_ind][max_number_blocks]>KL_bestMovement.max_found){
            KL_bestMovement.best_node_ind = node_ind;
            KL_bestMovement.best_comm_index = max_number_blocks;
            KL_bestMovement.max_found = Q_cols[node_ind][max_number_blocks];
        }
    }
    
    return KL_bestMovement;
}

void KL_IBN_optimization_testMovement_rows_updateClass(
                                                       int communityId,
                                                       double offsetTo_noneUpdateNodes,
                                                       vector<double> lambdas_rows,
                                                       vector<double> lambdas_cols,
                                                       vector<vector<int>> & M,
                                                       vector<int> & k_rows,
                                                       vector<int> & k_cols,
                                                       vector<int> partition_rows,
                                                       vector<int> partition_cols,
                                                       vector<vector<double>> & Q_rows,
                                                       int max_number_blocks){
    
    int N_rows=M.size();
    
    vector<int> partition_rows_temp = partition_rows;
    vector<int> partition_cols_temp = partition_cols;
    
    vector<double> lambdas_rows_temp = lambdas_rows;
    vector<double> lambdas_cols_temp = lambdas_cols;
    //compute increase in fitness for changing node
    for(int node_ind=0; node_ind < N_rows; node_ind++){
        for(int comm_index = 0; comm_index < max_number_blocks; comm_index++){
            //fast way
            partition_rows_temp = partition_rows;
            partition_cols_temp = partition_cols;
            lambdas_rows_temp = lambdas_rows;
            lambdas_cols_temp = lambdas_cols;
            
            //if the node belongs to the community to update we must recompute increase to
            //move it to all other communities
            if(partition_rows_temp[node_ind] == communityId || comm_index == communityId){
                if(partition_rows_temp[node_ind] != comm_index){
                    lambdas_inblock_change_node_partition_row(
                                                              node_ind,
                                                              partition_rows_temp[node_ind],
                                                              M,
                                                              k_cols,
                                                              k_rows,
                                                              partition_cols_temp,
                                                              partition_rows_temp,
                                                              lambdas_cols_temp,
                                                              lambdas_rows_temp,
                                                              REMOVE,
                                                              max_number_blocks);
                    
                    lambdas_inblock_change_node_partition_row(
                                                              node_ind,
                                                              comm_index,
                                                              M,
                                                              k_cols,
                                                              k_rows,
                                                              partition_cols_temp,
                                                              partition_rows_temp,
                                                              lambdas_cols_temp,
                                                              lambdas_rows_temp,
                                                              ADD,
                                                              max_number_blocks);
                    
                    Q_rows[node_ind][comm_index] = calculate_Fitness(
                                                                     M,
                                                                     k_cols,
                                                                     k_rows,
                                                                     lambdas_cols_temp,
                                                                     lambdas_rows_temp,
                                                                     IBN);
                    
                }else{
                    Q_rows[node_ind][comm_index] = -DBL_MAX;
                }
            }else{
                Q_rows[node_ind][comm_index] = Q_rows[node_ind][comm_index] + offsetTo_noneUpdateNodes;
            }
        }
        //fast way
        partition_rows_temp = partition_rows;
        partition_cols_temp = partition_cols;
        lambdas_rows_temp = lambdas_rows;
        lambdas_cols_temp = lambdas_cols;
        
        lambdas_inblock_change_node_partition_row(
                                                  node_ind,
                                                  partition_rows_temp[node_ind],
                                                  M,
                                                  k_cols,
                                                  k_rows,
                                                  partition_cols_temp,
                                                  partition_rows_temp,
                                                  lambdas_cols_temp,
                                                  lambdas_rows_temp,
                                                  REMOVE,
                                                  max_number_blocks);
        
        lambdas_inblock_change_node_partition_row(
                                                  node_ind,
                                                  max_number_blocks,
                                                  M,
                                                  k_cols,
                                                  k_rows,
                                                  partition_cols_temp,
                                                  partition_rows_temp,
                                                  lambdas_cols_temp,
                                                  lambdas_rows_temp,
                                                  ADD,
                                                  max_number_blocks+1);
        
        Q_rows[node_ind][max_number_blocks] = calculate_Fitness(
                                                                M,
                                                                k_cols,
                                                                k_rows,
                                                                lambdas_cols_temp,
                                                                lambdas_rows_temp,
                                                                IBN);
    }
}

void KL_IBN_optimization_testMovement_cols_updateClass(
                                                       int communityId,
                                                       double offsetTo_noneUpdateNodes,
                                                       vector<double> lambdas_rows,
                                                       vector<double> lambdas_cols,
                                                       vector<vector<int>> & M,
                                                       vector<int> & k_rows,
                                                       vector<int> & k_cols,
                                                       vector<int> partition_rows,
                                                       vector<int> partition_cols,
                                                       vector<vector<double>> & Q_cols,
                                                       int max_number_blocks){
    
    int N_cols=M[0].size();
    
    vector<int> partition_rows_temp = partition_rows;
    vector<int> partition_cols_temp = partition_cols;
    vector<double> lambdas_rows_temp = lambdas_rows;
    vector<double> lambdas_cols_temp = lambdas_cols;
    
    //compute increase in fitness for changing node
    for(int node_ind=0; node_ind < N_cols; node_ind++){
        for(int comm_index = 0; comm_index < max_number_blocks; comm_index++){
            
            partition_rows_temp = partition_rows;
            partition_cols_temp = partition_cols;
            lambdas_rows_temp = lambdas_rows;
            lambdas_cols_temp = lambdas_cols;
            
            //if the node belongs to the community to update we must recompute increase to
            //move it to all other communities
            if(partition_cols_temp[node_ind] == communityId || comm_index == communityId){
                if(partition_cols_temp[node_ind] != comm_index){
                    
                    //fast way
                    lambdas_inblock_change_node_partition_col(
                                                              node_ind,
                                                              partition_cols_temp[node_ind],
                                                              M,
                                                              k_cols,
                                                              k_rows,
                                                              partition_cols_temp,
                                                              partition_rows_temp,
                                                              lambdas_cols_temp,
                                                              lambdas_rows_temp,
                                                              REMOVE,
                                                              max_number_blocks);
                    
                    lambdas_inblock_change_node_partition_col(
                                                              node_ind,
                                                              comm_index,
                                                              M,
                                                              k_cols,
                                                              k_rows,
                                                              partition_cols_temp,
                                                              partition_rows_temp,
                                                              lambdas_cols_temp,
                                                              lambdas_rows_temp,
                                                              ADD,
                                                              max_number_blocks);
                    
                    Q_cols[node_ind][comm_index] = calculate_Fitness(
                                                                     M,
                                                                     k_cols,
                                                                     k_rows,
                                                                     lambdas_cols_temp,
                                                                     lambdas_rows_temp,
                                                                     IBN);
                    
                }else{
                    Q_cols[node_ind][comm_index] = -DBL_MAX;
                }
            }else{
                Q_cols[node_ind][comm_index] = Q_cols[node_ind][comm_index] + offsetTo_noneUpdateNodes;
            }
        }
        //or new community???
        partition_rows_temp = partition_rows;
        partition_cols_temp = partition_cols;
        lambdas_rows_temp = lambdas_rows;
        lambdas_cols_temp = lambdas_cols;
        
        lambdas_inblock_change_node_partition_col(
                                                  node_ind,
                                                  partition_cols_temp[node_ind],
                                                  M,
                                                  k_cols,
                                                  k_rows,
                                                  partition_cols_temp,
                                                  partition_rows_temp,
                                                  lambdas_cols_temp,
                                                  lambdas_rows_temp,
                                                  REMOVE,
                                                  max_number_blocks);
        
        lambdas_inblock_change_node_partition_col(
                                                  node_ind,
                                                  max_number_blocks,
                                                  M,
                                                  k_cols,
                                                  k_rows,
                                                  partition_cols_temp,
                                                  partition_rows_temp,
                                                  lambdas_cols_temp,
                                                  lambdas_rows_temp,
                                                  ADD,
                                                  max_number_blocks+1);
        
        Q_cols[node_ind][max_number_blocks] = calculate_Fitness(
                                                                M,
                                                                k_cols,
                                                                k_rows,
                                                                lambdas_cols_temp,
                                                                lambdas_rows_temp,
                                                                IBN);
    }
}



data_out_KL KL_optimization_IBN( vector<vector<int>> & M, vector<int> partition_rows, vector<int> partition_cols){
    data_out_KL var;
    //vector<vector<int> > M;
    //ifstream is(Filename);
    //load_matrix(&is, &M);
   // log_prefix_global = log_prefix;
    //load_matrix_sparseMatlab(is, M);
    
    //network degress by rows and columns
    int N_cols=M[0].size();
    int N_rows=M.size();
    vector<int> k_rows(N_rows,0);
    vector<int> k_cols(N_cols,0);
    for (int i = 0; i < N_rows; ++i) {
        for (int j = 0; j < N_cols; ++j) {
            k_rows[i]+=M[i][j];
            k_cols[j]+=M[i][j];
        }
    }
    
    //!!!!!int count = 0;
    // bool done = false;
    double current_IBN = 0;
    double current_IBN_old = 0;
    int oldCommunity = -1;
    int newCommunity = -1;
    
    double fitness_diff;
    
    //obtain required variables
    int max_number_blocks_rows=*max_element(partition_rows.begin(), partition_rows.end());
    int max_number_blocks_cols=*max_element(partition_cols.begin(), partition_cols.end());
    int max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows) + 1;
    
    //obtain initial fitness
    vector<vector<double> > lambdas = call_lambda_i(
                                                    M,
                                                    k_cols,
                                                    k_rows,
                                                    partition_cols,
                                                    partition_rows,
                                                    max_number_blocks,
                                                    IBN);
    
    current_IBN = calculate_Fitness(
                                    M,
                                    k_cols,
                                    k_rows,
                                    lambdas[0],
                                    lambdas[1],
                                    IBN);
    
    //sprintf(outputString,"fprintf('%s: Initial IBN = %f\\n');",log_prefix_global,current_IBN);
    //mexEvalString(outputString);
    
    vector<vector<double>> Q_rows;
    vector<vector<double>> Q_cols;
    
    Q_rows = vector<vector<double>> (N_rows, vector<double> (max_number_blocks+1,-DBL_MAX));
    Q_cols = vector<vector<double>> (N_cols, vector<double> (max_number_blocks+1,-DBL_MAX));
    
    KL_bestMovement_struct KL_bestMovement_rows = KL_IBN_optimization_testMovement_rows(
                                                                                        M,k_rows,k_cols,partition_rows,partition_cols,Q_rows,max_number_blocks);
    
    KL_bestMovement_struct KL_bestMovement_cols = KL_IBN_optimization_testMovement_cols(
                                                                                        M,k_rows,k_cols,partition_rows,partition_cols,Q_cols,max_number_blocks);
    
    current_IBN_old = current_IBN;
    double global_MAX_IBN = current_IBN;
    
    int countj = 0;
    while(  (global_MAX_IBN<KL_bestMovement_rows.max_found || global_MAX_IBN<KL_bestMovement_cols.max_found) && countj < 2*(N_rows+N_cols) ){
        countj++;
        // if (countj == (N_rows+N_cols))
        // {
        //     printf("Condition reached! %d KL IBN loops\n", countj);
        // }

        //update the required node partition
        if(KL_bestMovement_rows.max_found>=KL_bestMovement_cols.max_found){
            oldCommunity = partition_rows[KL_bestMovement_rows.best_node_ind];
            partition_rows[KL_bestMovement_rows.best_node_ind] = KL_bestMovement_rows.best_comm_index;
            newCommunity = partition_rows[KL_bestMovement_rows.best_node_ind];
            //sprintf(outputString,"fprintf('%s: count = %i,\\t rows %f from %i to %i \\t');",
              //      log_prefix_global,count,KL_bestMovement_rows.max_found,oldCommunity,newCommunity);
            //mexEvalString(outputString);
        }else{
            oldCommunity = partition_cols[KL_bestMovement_cols.best_node_ind];
            partition_cols[KL_bestMovement_cols.best_node_ind] = KL_bestMovement_cols.best_comm_index;
            newCommunity = partition_cols[KL_bestMovement_cols.best_node_ind];
            
            //sprintf(outputString,"fprintf('%s: count = %i,\\t col %f from %i to %i \\t');",
              //      log_prefix_global,count,KL_bestMovement_cols.max_found,oldCommunity,newCommunity);
            //mexEvalString(outputString);
        }
        
        //obtain required variables
        max_number_blocks_rows=*max_element(partition_rows.begin(), partition_rows.end());
        max_number_blocks_cols=*max_element(partition_cols.begin(), partition_cols.end());
        max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows) + 1;
        
        //obtain initial fitness
        vector<double> lambdas_rows;
        vector<double> lambdas_cols;
        lambdas = call_lambda_i(
                                M,
                                k_cols,
                                k_rows,
                                partition_cols,
                                partition_rows,
                                max_number_blocks,
                                IBN);
        
       // sprintf(outputString,"fprintf('old IBN =%f, ');",current_IBN);
        //mexEvalString(outputString);
        current_IBN = calculate_Fitness(
                                        M,
                                        k_cols,
                                        k_rows,
                                        lambdas[0],
                                        lambdas[1],
                                        IBN);

        if(global_MAX_IBN < current_IBN){
            global_MAX_IBN = current_IBN;
        }
        
        //sprintf(outputString,"fprintf('new IBN = %f\\n');",current_IBN);
        //mexEvalString(outputString);
        
        //if we increase the number of classes we need also to increase the Q_old's vectors
        if(Q_rows[0].size() < (static_cast<size_t>(max_number_blocks+1))){
          //  sprintf(outputString,"fprintf('increase number of communities in Q vectors\\n');");
            //mexEvalString(outputString);
            //increase the size of the vectors
            for(int i=0; i<N_rows; i++){
                Q_rows[i].resize(Q_rows[i].size()+1,0);
            }
            for(int i=0; i<N_cols; i++){
                Q_cols[i].resize(Q_cols[i].size()+1,0);
            }
        }
        
        //fast code
        lambdas_rows = lambdas[1];
        lambdas_cols = lambdas[0];
        fitness_diff = current_IBN - current_IBN_old;
        KL_IBN_optimization_testMovement_rows_updateClass(
                                                          oldCommunity,fitness_diff,
                                                          lambdas_rows,lambdas_cols,M,k_rows,k_cols,
                                                          partition_rows,partition_cols,Q_rows,
                                                          max_number_blocks);
        KL_IBN_optimization_testMovement_rows_updateClass(
                                                          newCommunity,0,
                                                          lambdas_rows,lambdas_cols,M,k_rows,k_cols,
                                                          partition_rows,partition_cols,Q_rows,
                                                          max_number_blocks);
        find_Qmax(KL_bestMovement_rows,
                  Q_rows,N_rows,max_number_blocks+1);
        
        KL_IBN_optimization_testMovement_cols_updateClass(
                                                          oldCommunity,fitness_diff,
                                                          lambdas_rows,lambdas_cols,M,k_rows,k_cols,
                                                          partition_rows,partition_cols,Q_cols,
                                                          max_number_blocks);
        KL_IBN_optimization_testMovement_cols_updateClass(
                                                          newCommunity,0,
                                                          lambdas_rows,lambdas_cols,M,k_rows,k_cols,
                                                          partition_rows,partition_cols,Q_cols,
                                                          max_number_blocks);
        find_Qmax(KL_bestMovement_cols,
                  Q_cols,N_cols,max_number_blocks+1);
        current_IBN_old = current_IBN;
        //!!!!count ++;
    }
    //printf("N_rows = %d\t N_cols = %d\t Number of loops (IBN): %d\n", N_rows, N_cols, countj);
    
    var.partitions_result.push_back(partition_rows);
    var.partitions_result.push_back(partition_cols);
    var.QI_kl = current_IBN;
    return var;
}

KL_bestMovement_struct KL_MOD_optimization_testMovement_cols(
                                                             vector<vector<int>> & M,
                                                             vector<int> & k_rows,
                                                             vector<int> & k_cols,
                                                             vector<int> partition_rows,
                                                             vector<int> partition_cols,
                                                             vector<vector<double>> & Q_cols,
                                                             int max_number_blocks){
    
    int N_cols=M[0].size();
    //int N_rows=M.size();
    
    //compute the original lambdas
    vector<vector<double> > lambdas = call_lambda_i(
                                                    M,
                                                    k_cols,
                                                    k_rows,
                                                    partition_cols,
                                                    partition_rows,
                                                    max_number_blocks,
                                                    MOD);
    vector<double> lambdas_cols = lambdas[0];
    vector<double> lambdas_rows = lambdas[1];
    vector<double> lambdas_cols_temp;
    vector<double> lambdas_rows_temp;
    
    vector<int> partition_rows_temp = partition_rows;
    vector<int> partition_cols_temp = partition_cols;
    KL_bestMovement_struct KL_bestMovement;
    KL_bestMovement.best_node_ind = -1;
    KL_bestMovement.best_comm_index = -1;
    KL_bestMovement.max_found = -DBL_MAX;
    //compute increase in fitness for changing node
    for(int node_ind=0; node_ind < N_cols; node_ind++){
        for(int comm_index = 0; comm_index < max_number_blocks; comm_index++){
            
            if(partition_cols_temp[node_ind] != comm_index){
                
                
                //fast way
                partition_rows_temp = partition_rows;
                partition_cols_temp = partition_cols;
                lambdas_rows_temp = lambdas_rows;
                lambdas_cols_temp = lambdas_cols;
                
                lambdas_modularity_change_node_col(
                                                   node_ind,
                                                   partition_cols_temp[node_ind],
                                                   M,
                                                   k_cols,
                                                   k_rows,
                                                   partition_cols_temp,
                                                   partition_rows_temp,
                                                   lambdas_cols_temp,
                                                   lambdas_rows_temp,
                                                   REMOVE,
                                                   max_number_blocks);
                
                lambdas_modularity_change_node_col(
                                                   node_ind,
                                                   comm_index,
                                                   M,
                                                   k_cols,
                                                   k_rows,
                                                   partition_cols_temp,
                                                   partition_rows_temp,
                                                   lambdas_cols_temp,
                                                   lambdas_rows_temp,
                                                   ADD,
                                                   max_number_blocks);
                
                Q_cols[node_ind][comm_index] = calculate_Fitness(
                                                                 M,
                                                                 k_cols,
                                                                 k_rows,
                                                                 lambdas_cols_temp,
                                                                 lambdas_rows_temp,
                                                                 MOD);
                
            }else{
                Q_cols[node_ind][comm_index] = -DBL_MAX;
            }
            
            if(Q_cols[node_ind][comm_index]>KL_bestMovement.max_found){
                KL_bestMovement.best_node_ind = node_ind;
                KL_bestMovement.best_comm_index = comm_index;
                KL_bestMovement.max_found = Q_cols[node_ind][comm_index];
            }
        }
        //or new community???
        //fast way
        partition_rows_temp = partition_rows;
        partition_cols_temp = partition_cols;
        lambdas_rows_temp = lambdas_rows;
        lambdas_cols_temp = lambdas_cols;
        
        lambdas_modularity_change_node_col(
                                           node_ind,
                                           partition_cols_temp[node_ind],
                                           M,
                                           k_cols,
                                           k_rows,
                                           partition_cols_temp,
                                           partition_rows_temp,
                                           lambdas_cols_temp,
                                           lambdas_rows_temp,
                                           REMOVE,
                                           max_number_blocks);
        
        lambdas_modularity_change_node_col(
                                           node_ind,
                                           max_number_blocks,
                                           M,
                                           k_cols,
                                           k_rows,
                                           partition_cols_temp,
                                           partition_rows_temp,
                                           lambdas_cols_temp,
                                           lambdas_rows_temp,
                                           ADD,
                                           max_number_blocks+1);
        
        Q_cols[node_ind][max_number_blocks] = calculate_Fitness(
                                                                M,
                                                                k_cols,
                                                                k_rows,
                                                                lambdas_cols_temp,
                                                                lambdas_rows_temp,
                                                                MOD);
        
        if(Q_cols[node_ind][max_number_blocks]>KL_bestMovement.max_found){
            KL_bestMovement.best_node_ind = node_ind;
            KL_bestMovement.best_comm_index = max_number_blocks;
            KL_bestMovement.max_found = Q_cols[node_ind][max_number_blocks];
        }
    }
    
    return KL_bestMovement;
}

KL_bestMovement_struct KL_MOD_optimization_testMovement_rows(
                                                             vector<vector<int>> & M,
                                                             vector<int> & k_rows,
                                                             vector<int> & k_cols,
                                                             vector<int> partition_rows,
                                                             vector<int> partition_cols,
                                                             vector<vector<double>> & Q_rows,
                                                             int max_number_blocks){
    
    //int N_cols=M[0].size();
    int N_rows=M.size();
    
    //compute the original lambdas
    vector<vector<double> > lambdas = call_lambda_i(
                                                    M,
                                                    k_cols,
                                                    k_rows,
                                                    partition_cols,
                                                    partition_rows,
                                                    max_number_blocks,
                                                    MOD);
    vector<double> lambdas_cols = lambdas[0];
    vector<double> lambdas_rows = lambdas[1];
    vector<double> lambdas_cols_temp;
    vector<double> lambdas_rows_temp;
    
    vector<int> partition_rows_temp = partition_rows;
    vector<int> partition_cols_temp = partition_cols;
    KL_bestMovement_struct KL_bestMovement;
    KL_bestMovement.best_node_ind = -1;
    KL_bestMovement.best_comm_index = -1;
    KL_bestMovement.max_found = -DBL_MAX;
    //compute increase in fitness for changing node
    for(int node_ind=0; node_ind < N_rows; node_ind++){
        for(int comm_index = 0; comm_index < max_number_blocks; comm_index++){
            
            //sprintf(outputString,"fprintf('p1 - %i,%i\\n');",node_ind,comm_index);
            //mexEvalString(outputString);
            
            if(partition_rows_temp[node_ind] != comm_index){
                
                //fast way
                partition_rows_temp = partition_rows;
                partition_cols_temp = partition_cols;
                lambdas_rows_temp = lambdas_rows;
                lambdas_cols_temp = lambdas_cols;
                
                lambdas_modularity_change_node_row(
                                                   node_ind,
                                                   partition_rows_temp[node_ind],
                                                   M,
                                                   k_cols,
                                                   k_rows,
                                                   partition_cols_temp,
                                                   partition_rows_temp,
                                                   lambdas_cols_temp,
                                                   lambdas_rows_temp,
                                                   REMOVE,
                                                   max_number_blocks);
                
                lambdas_modularity_change_node_row(
                                                   node_ind,
                                                   comm_index,
                                                   M,
                                                   k_cols,
                                                   k_rows,
                                                   partition_cols_temp,
                                                   partition_rows_temp,
                                                   lambdas_cols_temp,
                                                   lambdas_rows_temp,
                                                   ADD,
                                                   max_number_blocks);
                
                Q_rows[node_ind][comm_index] = calculate_Fitness(
                                                                 M,
                                                                 k_cols,
                                                                 k_rows,
                                                                 lambdas_cols_temp,
                                                                 lambdas_rows_temp,
                                                                 MOD);
                
            }else{
                Q_rows[node_ind][comm_index] = -DBL_MAX;
            }
            
            //sprintf(outputString,"fprintf('\\t p2\\n');");
            //mexEvalString(outputString);
            
            if(Q_rows[node_ind][comm_index]>KL_bestMovement.max_found){
                KL_bestMovement.best_node_ind = node_ind;
                KL_bestMovement.best_comm_index = comm_index;
                KL_bestMovement.max_found = Q_rows[node_ind][comm_index];
            }
            //sprintf(outputString,"fprintf('\\t p3\\n');");
            //mexEvalString(outputString);
            
        }
        
        //sprintf(outputString,"fprintf('\\t p4\\n');");
        //mexEvalString(outputString);
        
        //or new community???
        //fast way
        partition_rows_temp = partition_rows;
        partition_cols_temp = partition_cols;
        lambdas_rows_temp = lambdas_rows;
        lambdas_cols_temp = lambdas_cols;
        
        lambdas_modularity_change_node_row(
                                           node_ind,
                                           partition_rows_temp[node_ind],
                                           M,
                                           k_cols,
                                           k_rows,
                                           partition_cols_temp,
                                           partition_rows_temp,
                                           lambdas_cols_temp,
                                           lambdas_rows_temp,
                                           REMOVE,
                                           max_number_blocks);
        
        lambdas_modularity_change_node_row(
                                           node_ind,
                                           max_number_blocks,
                                           M,
                                           k_cols,
                                           k_rows,
                                           partition_cols_temp,
                                           partition_rows_temp,
                                           lambdas_cols_temp,
                                           lambdas_rows_temp,
                                           ADD,
                                           max_number_blocks+1);
        
        Q_rows[node_ind][max_number_blocks] = calculate_Fitness(
                                                                M,
                                                                k_cols,
                                                                k_rows,
                                                                lambdas_cols_temp,
                                                                lambdas_rows_temp,
                                                                MOD);
        
        if(Q_rows[node_ind][max_number_blocks]>KL_bestMovement.max_found){
            KL_bestMovement.best_node_ind = node_ind;
            KL_bestMovement.best_comm_index = max_number_blocks;
            KL_bestMovement.max_found = Q_rows[node_ind][max_number_blocks];
        }
    }
    
    return KL_bestMovement;
}

void KL_MOD_optimization_testMovement_rows_updateClass(
                                                       int communityId,
                                                       double offsetTo_noneUpdateNodes,
                                                       vector<double> lambdas_rows,
                                                       vector<double> lambdas_cols,
                                                       vector<vector<int>> & M,
                                                       vector<int> & k_rows,
                                                       vector<int> & k_cols,
                                                       vector<int> partition_rows,
                                                       vector<int> partition_cols,
                                                       vector<vector<double>> & Q_rows,
                                                       int max_number_blocks){
    
    //int N_cols=M[0].size();
    int N_rows=M.size();
    
    vector<int> partition_rows_temp = partition_rows;
    vector<int> partition_cols_temp = partition_cols;
    
    vector<double> lambdas_cols_temp = lambdas_cols;
    vector<double> lambdas_rows_temp = lambdas_rows;
    
    //compute increase in fitness for changing node
    for(int node_ind=0; node_ind < N_rows; node_ind++){
        for(int comm_index = 0; comm_index < max_number_blocks; comm_index++){
            //fast way
            partition_rows_temp = partition_rows;
            partition_cols_temp = partition_cols;
            lambdas_rows_temp = lambdas_rows;
            lambdas_cols_temp = lambdas_cols;
            
            //if the node belongs to the community to update we must recompute increase to
            //move it to all other communities
            if(partition_rows_temp[node_ind] == communityId || comm_index == communityId){
                if(partition_rows_temp[node_ind] != comm_index){
                    
                    lambdas_modularity_change_node_row(
                                                       node_ind,
                                                       partition_rows_temp[node_ind],
                                                       M,
                                                       k_cols,
                                                       k_rows,
                                                       partition_cols_temp,
                                                       partition_rows_temp,
                                                       lambdas_cols_temp,
                                                       lambdas_rows_temp,
                                                       REMOVE,
                                                       max_number_blocks);
                    
                    lambdas_modularity_change_node_row(
                                                       node_ind,
                                                       comm_index,
                                                       M,
                                                       k_cols,
                                                       k_rows,
                                                       partition_cols_temp,
                                                       partition_rows_temp,
                                                       lambdas_cols_temp,
                                                       lambdas_rows_temp,
                                                       ADD,
                                                       max_number_blocks);
                    
                    Q_rows[node_ind][comm_index] = calculate_Fitness(
                                                                     M,
                                                                     k_cols,
                                                                     k_rows,
                                                                     lambdas_cols_temp,
                                                                     lambdas_rows_temp,
                                                                     MOD);
                    
                }else{
                    Q_rows[node_ind][comm_index] = -DBL_MAX;
                }
            }else{
                Q_rows[node_ind][comm_index] = Q_rows[node_ind][comm_index] + offsetTo_noneUpdateNodes;
            }
        }
        //fast way
        partition_rows_temp = partition_rows;
        partition_cols_temp = partition_cols;
        lambdas_rows_temp = lambdas_rows;
        lambdas_cols_temp = lambdas_cols;
        
        lambdas_modularity_change_node_row(
                                           node_ind,
                                           partition_rows_temp[node_ind],
                                           M,
                                           k_cols,
                                           k_rows,
                                           partition_cols_temp,
                                           partition_rows_temp,
                                           lambdas_cols_temp,
                                           lambdas_rows_temp,
                                           REMOVE,
                                           max_number_blocks);
        
        lambdas_modularity_change_node_row(
                                           node_ind,
                                           max_number_blocks,
                                           M,
                                           k_cols,
                                           k_rows,
                                           partition_cols_temp,
                                           partition_rows_temp,
                                           lambdas_cols_temp,
                                           lambdas_rows_temp,
                                           ADD,
                                           max_number_blocks+1);
        
        Q_rows[node_ind][max_number_blocks] = calculate_Fitness(
                                                                M,
                                                                k_cols,
                                                                k_rows,
                                                                lambdas_cols_temp,
                                                                lambdas_rows_temp,
                                                                MOD);
    }
}

void KL_MOD_optimization_testMovement_cols_updateClass(
                                                       int communityId,
                                                       double offsetTo_noneUpdateNodes,
                                                       vector<double> lambdas_rows,
                                                       vector<double> lambdas_cols,
                                                       vector<vector<int>> & M,
                                                       vector<int> & k_rows,
                                                       vector<int> & k_cols,
                                                       vector<int> partition_rows,
                                                       vector<int> partition_cols,
                                                       vector<vector<double>> & Q_cols,
                                                       int max_number_blocks){
    
    int N_cols=M[0].size();
   // int N_rows=M.size();
    
    vector<int> partition_rows_temp = partition_rows;
    vector<int> partition_cols_temp = partition_cols;
    
    vector<double> lambdas_cols_temp = lambdas_cols;
    vector<double> lambdas_rows_temp = lambdas_rows;
    
    //compute increase in fitness for changing node
    for(int node_ind=0; node_ind < N_cols; node_ind++){
        for(int comm_index = 0; comm_index < max_number_blocks; comm_index++){
            
            partition_rows_temp = partition_rows;
            partition_cols_temp = partition_cols;
            lambdas_rows_temp = lambdas_rows;
            lambdas_cols_temp = lambdas_cols;
            
            if(partition_cols_temp[node_ind] == communityId || comm_index == communityId){
                if(partition_cols_temp[node_ind] != comm_index){
                    
                    lambdas_modularity_change_node_col(
                                                       node_ind,
                                                       partition_cols_temp[node_ind],
                                                       M,
                                                       k_cols,
                                                       k_rows,
                                                       partition_cols_temp,
                                                       partition_rows_temp,
                                                       lambdas_cols_temp,
                                                       lambdas_rows_temp,
                                                       REMOVE,
                                                       max_number_blocks);
                    
                    lambdas_modularity_change_node_col(
                                                       node_ind,
                                                       comm_index,
                                                       M,
                                                       k_cols,
                                                       k_rows,
                                                       partition_cols_temp,
                                                       partition_rows_temp,
                                                       lambdas_cols_temp,
                                                       lambdas_rows_temp,
                                                       ADD,
                                                       max_number_blocks);
                    
                    Q_cols[node_ind][comm_index] = calculate_Fitness(
                                                                     M,
                                                                     k_cols,
                                                                     k_rows,
                                                                     lambdas_cols_temp,
                                                                     lambdas_rows_temp,
                                                                     MOD);
                    
                }else{
                    Q_cols[node_ind][comm_index] = -DBL_MAX;
                }
            }else{
                Q_cols[node_ind][comm_index] = Q_cols[node_ind][comm_index] + offsetTo_noneUpdateNodes;
            }
        }
        //or new community???
        partition_rows_temp = partition_rows;
        partition_cols_temp = partition_cols;
        lambdas_rows_temp = lambdas_rows;
        lambdas_cols_temp = lambdas_cols;
        
        lambdas_modularity_change_node_col(
                                           node_ind,
                                           partition_cols_temp[node_ind],
                                           M,
                                           k_cols,
                                           k_rows,
                                           partition_cols_temp,
                                           partition_rows_temp,
                                           lambdas_cols_temp,
                                           lambdas_rows_temp,
                                           REMOVE,
                                           max_number_blocks);
        
        lambdas_modularity_change_node_col(
                                           node_ind,
                                           max_number_blocks,
                                           M,
                                           k_cols,
                                           k_rows,
                                           partition_cols_temp,
                                           partition_rows_temp,
                                           lambdas_cols_temp,
                                           lambdas_rows_temp,
                                           ADD,
                                           max_number_blocks+1);
        
        Q_cols[node_ind][max_number_blocks] = calculate_Fitness(
                                                                M,
                                                                k_cols,
                                                                k_rows,
                                                                lambdas_cols_temp,
                                                                lambdas_rows_temp,
                                                                MOD);
    }
}


data_out_KL KL_optimization_MOD( vector<vector<int>> & M, vector<int> partition_rows, vector<int> partition_cols){
    data_out_KL var;
    //vector<vector<int> > M;
   // ifstream is(Filename);
   // load_matrix(&is, &M);
    // log_prefix_global = log_prefix;
    // load_matrix_sparseMatlab(is, M);
    
    //network degress by rows and columns
    int N_cols=M[0].size();
    int N_rows=M.size();
    vector<int> k_rows(N_rows,0);
    vector<int> k_cols(N_cols,0);
    for (int i = 0; i < N_rows; ++i) {
        for (int j = 0; j < N_cols; ++j) {
            k_rows[i]+=M[i][j];
            k_cols[j]+=M[i][j];
        }
    }
    
    //!!!!int count = 0;
    // bool done = false;
    double current_MOD = 0;
    double current_MOD_old = 0;
    int oldCommunity = -1;
    int newCommunity = -1;
    
    vector<double> lambdas_rows;
    vector<double> lambdas_cols;
    double fitness_diff;
    
    //obtain required variables
    int max_number_blocks_cols=*max_element(partition_rows.begin(), partition_rows.end());
    int max_number_blocks_rows=*max_element(partition_cols.begin(), partition_cols.end());
    int max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows) + 1;
    
    //obtain initial fitness
    vector<vector<double> > lambdas = call_lambda_i(
                                                    M,
                                                    k_cols,
                                                    k_rows,
                                                    partition_cols,
                                                    partition_rows,
                                                    max_number_blocks,
                                                    MOD);
    
    current_MOD = calculate_Fitness(
                                    M,
                                    k_cols,
                                    k_rows,
                                    lambdas[0],
                                    lambdas[1],
                                    MOD);
    
    //sprintf(outputString,"fprintf('%s: Initial MOD = %f\\n');",log_prefix_global,current_MOD);
    //mexEvalString(outputString);
    
    vector<vector<double>> Q_rows;
    vector<vector<double>> Q_cols;
    
    Q_rows = vector<vector<double>> (N_rows, vector<double> (max_number_blocks+1,-DBL_MAX));
    Q_cols = vector<vector<double>> (N_cols, vector<double> (max_number_blocks+1,-DBL_MAX));
    
    KL_bestMovement_struct KL_bestMovement_rows = KL_MOD_optimization_testMovement_rows(
                                                                                        M,k_rows,k_cols,partition_rows,partition_cols,Q_rows,max_number_blocks);
    
    KL_bestMovement_struct KL_bestMovement_cols = KL_MOD_optimization_testMovement_cols(
                                                                                        M,k_rows,k_cols,partition_rows,partition_cols,Q_cols,max_number_blocks);
    
    current_MOD_old = current_MOD;
    double global_MAX_MOD = current_MOD;
    
    int countj = 0;
    while( (global_MAX_MOD<KL_bestMovement_rows.max_found || global_MAX_MOD<KL_bestMovement_cols.max_found) && countj < 2*(N_rows+N_cols) ){
        countj++; // allows only 2*(N_rows+N_cols) attempts to improve
        // if (countj == (N_rows+N_cols))
        // {
        //     printf("Condition reached! %d KL MOD loops\n", countj);
        // }

        //update the required node partition
        if(KL_bestMovement_rows.max_found>=KL_bestMovement_cols.max_found){
            oldCommunity = partition_rows[KL_bestMovement_rows.best_node_ind];
            partition_rows[KL_bestMovement_rows.best_node_ind] = KL_bestMovement_rows.best_comm_index;
            newCommunity = partition_rows[KL_bestMovement_rows.best_node_ind];
            
            // sprintf(outputString,"fprintf('%s: count = %i,\\t rows %i from %i to %i \\t');",
            //       log_prefix_global,count,KL_bestMovement_rows.best_node_ind,oldCommunity,newCommunity);
            //mexEvalString(outputString);
        }else{
            oldCommunity = partition_cols[KL_bestMovement_cols.best_node_ind];
            partition_cols[KL_bestMovement_cols.best_node_ind] = KL_bestMovement_cols.best_comm_index;
            newCommunity = partition_cols[KL_bestMovement_cols.best_node_ind];
            
            // sprintf(outputString,"fprintf('%s: count = %i,\\t col %i from %i to %i \\t');",
            //       log_prefix_global,count,KL_bestMovement_cols.best_node_ind,oldCommunity,newCommunity);
            //mexEvalString(outputString);
        }
        
        //obtain required variables
        max_number_blocks_cols=*max_element(partition_rows.begin(), partition_rows.end());
        max_number_blocks_rows=*max_element(partition_cols.begin(), partition_cols.end());
        max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows) + 1;
        
        //obtain initial fitness
        lambdas = call_lambda_i(
                                M,
                                k_cols,
                                k_rows,
                                partition_cols,
                                partition_rows,
                                max_number_blocks,
                                MOD);
        
        // sprintf(outputString,"fprintf('old MOD =%f, ');",current_MOD);
        //mexEvalString(outputString);
        
        current_MOD = calculate_Fitness(
                                        M,
                                        k_cols,
                                        k_rows,
                                        lambdas[0],
                                        lambdas[1],
                                        MOD);
        
        if(global_MAX_MOD < current_MOD){
            global_MAX_MOD = current_MOD;
        }
        
        //  sprintf(outputString,"fprintf('new MOD = %f\\n');",current_MOD);
        //mexEvalString(outputString);
        
        //if we increase the number of classes we need also to increase the Q_old's vectors
        if(Q_rows[0].size() < (static_cast<size_t>(max_number_blocks+1))){
            // sprintf(outputString,"fprintf('increase number of communities in Q vectors\\n');");
            //mexEvalString(outputString);
            //increase the size of the vectors
            for(int i=0; i<N_rows; i++){
                Q_rows[i].resize(Q_rows[i].size()+1,0);
            }
            for(int i=0; i<N_cols; i++){
                Q_cols[i].resize(Q_cols[i].size()+1,0);
            }
        }
        
        //fast code
        lambdas_rows = lambdas[1];
        lambdas_cols = lambdas[0];
        fitness_diff = current_MOD - current_MOD_old;
        
        KL_MOD_optimization_testMovement_rows_updateClass(
                                                          oldCommunity,fitness_diff,
                                                          lambdas_rows,lambdas_cols,M,k_rows,k_cols,
                                                          partition_rows,partition_cols,Q_rows,
                                                          max_number_blocks);
        KL_MOD_optimization_testMovement_rows_updateClass(
                                                          newCommunity,0,
                                                          lambdas_rows,lambdas_cols,M,k_rows,k_cols,
                                                          partition_rows,partition_cols,Q_rows,
                                                          max_number_blocks);
        find_Qmax(KL_bestMovement_rows,
                  Q_rows,N_rows,max_number_blocks+1);
        
        KL_MOD_optimization_testMovement_cols_updateClass(
                                                          oldCommunity,fitness_diff,
                                                          lambdas_rows,lambdas_cols,M,k_rows,k_cols,
                                                          partition_rows,partition_cols,Q_cols,
                                                          max_number_blocks);
        KL_MOD_optimization_testMovement_cols_updateClass(
                                                          newCommunity,0,
                                                          lambdas_rows,lambdas_cols,M,k_rows,k_cols,
                                                          partition_rows,partition_cols,Q_cols,
                                                          max_number_blocks);
        find_Qmax(KL_bestMovement_cols,
                  Q_cols,N_cols,max_number_blocks+1);
        
        current_MOD_old = current_MOD;
        //!!!!count ++;
    }
    
    //printf("N_rows = %d\t N_cols = %d\t Number of loops (MOD): %d\n", N_rows, N_cols, countj);

    var.partitions_result.push_back(partition_rows);
    var.partitions_result.push_back(partition_cols);
    var.QI_kl = current_MOD;
    return var;
}

PYBIND11_MODULE(kernighan_lin_bi,m) {
    //py::module m("example", "Generating primes in c++ with python bindings using pybind11");
    py::class_<data_out_KL> dp(m, "data_out_KL");
    dp.def(py::init<>());
    m.def("call_lambda_i", &call_lambda_i, "A function that calculates the fitness contribution of each node");
    m.def("calculate_Fitness", &calculate_Fitness, "A function that calculate the Ieo and Qeo of the whole network");
    dp.def_readwrite("KL_Fitness", &data_out_KL::QI_kl);
    dp.def_readwrite("KL_partitions", &data_out_KL::partitions_result);
    m.def("KL_optimization_MOD", &KL_optimization_MOD, "A function that perform KL algoritihm to modularity");
    m.def("KL_optimization_IBN", &KL_optimization_IBN, "A function that perform KL algoritihm to in-block nestedness");
    //vector<vector<int> >
    //m.def("KL_metric", py::overload_cast<double>(&data_out_KL::QI_kl));
   // m.def("KL_partitions",(&data_out_KL::partitions_result));
    // return m.ptr();
}
// compilation on linux g++ -O3 -Wall -shared -std=c++11 -fPIC `python -m pybind11 --includes` KL_functions.cpp -o kernighan_lin_bi.so

//compilation on mac: g++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python -m pybind11 --includes` KL_functions.cpp -o kernighan_lin_bi.so
