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
#include <cfloat>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/stl_bind.h>
//#include "EO.hpp"

using namespace std;

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) < (Y)) ? (Y) : (X))
#define REMOVE true
#define ADD false
#define IBN true
#define MOD false

vector<int> indexes_sorted_vector(vector<double> & x){
    
    std::vector<int> vs(x.size());
    std::iota(vs.begin(), vs.end(), 0);
    auto comparator = [&x](int a, int b){ return x[a] < x[b]; };
    std::sort(vs.begin(), vs.end(), comparator);
    return vs;
}
void shuffle_partition(vector<int> &labels){
    
    int half_size = round(labels.size()/2);
    for (int i=0; i<labels.size();i++){
        if (i < half_size){
            labels[i]=1;
        }else{
            labels[i]=2;
        }
    }
    
    /* Shuffling the labels*/
    random_device rd; //non fixed seed
    mt19937 gen(rd());
    shuffle ( labels.begin(), labels.end(), gen);
    //log_file("\tleaving bipartition function");
}

/* Function to perform the bipartition of a given array.
 It receives an NxM matrix and two empty vectors where the label nodes will be storage*/
void bipartition(
                 vector<vector<int> > & input_matrix,
                 vector<int> & label_col,
                 vector<int> & label_row,
                 int nextBlockId){
    
    if(label_col.size()>2){
        int rnd_label = 0;
        int num_label1 = 0;
        int num_label2 = 0;
        for (int i=0 ; i<label_col.size() ; i++){
            rnd_label = rand() % 2;
            label_col[i] = rand() % 2 + nextBlockId;
            
            if(rnd_label == 1){
                num_label1 ++;
            }else{
                num_label2 ++;
            }
        }
        if(num_label1 == 0 || num_label2 == 0){
            
            shuffle_partition(label_col);
        }
    }else{
        label_col[0] = nextBlockId;
        label_col[1] = nextBlockId + 1;
    }
    
    if(label_row.size()>2){
        int rnd_label = 0;
        int num_label1 = 0;
        int num_label2 = 0;
        for (int i=0 ; i<label_row.size() ; i++){
            rnd_label = rand() % 2;
            label_row[i] = rand() % 2 + nextBlockId;
            
            if(rnd_label == 1){
                num_label1 ++;
            }else{
                num_label2 ++;
            }
        }
        if(num_label1 == 0 || num_label2 == 0){
            shuffle_partition(label_row);
        }
    }else{
        label_row[0] = nextBlockId;
        label_row[1] = nextBlockId + 1;
    }
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
    
    for (int j = 0; j < label_cols.size(); ++j){
        list_of_blocks_cols[label_cols[j]].push_back(j);
    }
    
    //rows
    for (int j = 0; j < label_rows.size(); ++j){
        list_of_blocks_rows[label_rows[j]].push_back(j);
    }
    
    //computing the column nodes contribution (pair overlap)
    for (int l=0; l<list_of_blocks_cols.size(); l++){
        int size_blocks_cols=(int)list_of_blocks_cols[l].size();
        
        
        for (int i=0; i<size_blocks_cols; i++){
            int ii = list_of_blocks_cols[l][i];
            
            
            for (int j=0; j<size_blocks_cols; j++){
                
                int jj = list_of_blocks_cols[l][j];
                double PO_col_i=0;
                double normalization = ((double)(k_cols[jj]))*((double)(size_blocks_cols-1.0));
                
                if ((k_cols[ii]>k_cols[jj]) & (k_cols[jj]>0) & (normalization!=0.0)){
                    
                    for (int k=0; k<list_of_blocks_rows[l].size(); k++){
                        
                        int kk = list_of_blocks_rows[l][k];
                        
                        if ((input_matrix[kk][ii]*input_matrix[kk][jj])==1){
                            PO_col_i++;
                        }
                    }
                    
                			double null_model = ((double)(k_cols[ii]*k_cols[jj]))/((double)(N_rows));
                			PO_col_i=(PO_col_i-null_model)/normalization;
                			//lambda_cols[ii]+=PO_col_i;
                		
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
    for (int l=0; l<list_of_blocks_rows.size(); l++){
        
        int size_blocks_rows=(int)list_of_blocks_rows[l].size();
        for (int i=0; i<(size_blocks_rows); i++){
            
            int ii = list_of_blocks_rows[l][i];
            for (int j=0; j<size_blocks_rows; j++){
                
                int jj =list_of_blocks_rows[l][j];
                double PO_row_i=0;
                double normalization = ((double)(k_rows[jj]))*((double)(size_blocks_rows-1.0));
                
                if ((k_rows[ii]>k_rows[jj]) & (k_rows[jj]>0) & (normalization!=0.0) ){
                    
                    for (int k=0; k<list_of_blocks_cols[l].size(); k++){
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
    for (int j = 0; j < label_cols.size(); ++j){
        list_of_blocks_cols[label_cols[j]].push_back(j);
    }
    
    //rows
    for (int j = 0; j < label_rows.size(); ++j){
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
            
            if ( (k_rows[ii]>k_rows[jj]) & (normalization!=0.0)){
                
                for (int k=0; k<list_of_blocks_cols[l].size(); k++){
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
            for (int k=0; k<list_of_blocks_cols[l].size(); k++){
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
        
        if ( (k_rows[ii]>=k_rows[jj]) ){
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
            
            if ((k_cols[ii]>k_cols[jj])  ){
                double PO_col_i=0;
                double delta_value;
                int kk = nodeId;
                if ((input_matrix[kk][ii]*input_matrix[kk][jj])==1){
                    PO_col_i++;
                }
                
                delta_value = PO_col_i/normalization;
                
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
    for (int j = 0; j < label_cols.size(); ++j){
        list_of_blocks_cols[label_cols[j]].push_back(j);
    }
    
    //rows
    for (int j = 0; j < label_rows.size(); ++j){
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
                
                for (int k=0; k<list_of_blocks_rows[l].size(); k++){
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
            
            for (int k=0; k<list_of_blocks_rows[l].size(); k++){
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
                double delta_value;
                
                int kk = nodeId;
                if ((input_matrix[ii][kk]*input_matrix[jj][kk])==1){
                    PO_row_i++;
                }
        
                    delta_value = PO_row_i/normalization;
                
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
            if  ((input_matrix[i][j]==1) & (label_rows[i]==label_cols[j])){
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

void first_step(
                vector<vector<int> > &input_matrix,
                vector<int> & k_cols,
                vector<int> & k_rows,
                vector<int> label_cols,
                vector<int> label_rows,
                int Malpha,
                vector<int> & label_cols_new,
                vector<int> & label_rows_new,
                int blockId1,
                int blockId2,
                int max_number_blocks,
                bool ibn){
    /*Given a matrix arbitrarily partitioned in two, move the nodes with "lowest" fitness
     until optimal Ieo is reached. The stopping criteria for the optimal Ieo is given by
     Malpha.*/
  
    vector<double> lambda_cols(label_cols.size(),0);
    vector<double> lambda_rows(label_rows.size(),0);
    vector<double> lambda_cols_copy(label_cols.size(),0);
    vector<double> lambda_rows_copy(label_rows.size(),0);
    
    vector<int> label_temporal_cols;
    vector<int> label_temporal_rows;
    double Io=-1.0;
    double In=-1.0;
    int alpha=0;
    int max_alpha=0;
    int low_node=0;
    
    //compute the first lambda_i
    label_temporal_cols = label_cols;
    label_temporal_rows = label_rows;
    
    //sprintf(outputString,"fprintf('inside current_block->first_step->alpha %i, Malpha %i, max_blocks %i\\n');",
    //        alpha,Malpha,max_number_blocks);
    //mexEvalString(outputString);
    
    lambda_i(
             input_matrix,
             k_cols,
             k_rows,
             label_temporal_cols,
             label_temporal_rows,
             lambda_cols,
             lambda_rows,
             max_number_blocks,
             ibn);
    
    while (alpha<Malpha){
        
        label_temporal_cols = label_cols;
        label_temporal_rows = label_rows;
        
        lambda_cols_copy = lambda_cols;
        lambda_rows_copy = lambda_rows;
        
        //sprintf(outputString,"fprintf('inside current_block->first_step->alpha %i, Malpha %i \\n');",alpha,Malpha);
        //mexEvalString(outputString);
        
        //sprintf(outputString,"fprintf('inside current_block->first_step->lambda_i \\n');");
        //mexEvalString(outputString);
        
        //sprintf(outputString,"fprintf('inside current_block->first_step->eval low_fitness_node\\n');");
        //mexEvalString(outputString);
        
        low_node=low_fitness_node(lambda_cols,lambda_rows);
        
        //sprintf(outputString,"fprintf('inside current_block->first_step->low_fitness_node %i \\n');",low_node);
        //mexEvalString(outputString);
        
        //TODO: do this update of lambda_i's for IBN in a nicer way
        if(ibn){
            // changing the node label
            // test instead of the if, replace the operation by using booleans
            //consider low_node ranges from 0 to (N+M)-1
            //sprintf(outputString,"fprintf('updating lambdas\\n');");
            //mexEvalString(outputString);
            
            if (low_node<label_temporal_rows.size()){
                //sprintf(outputString,"fprintf('begin updating lambdas 1\\n');");
                //mexEvalString(outputString);
                if (label_temporal_rows[low_node]==blockId1){
                    
                    lambdas_inblock_change_node_partition_row(
                                                              low_node,blockId1,input_matrix,k_cols,k_rows,
                                                              label_temporal_cols,label_temporal_rows,lambda_cols,lambda_rows,true,max_number_blocks);
                    lambdas_inblock_change_node_partition_row(
                                                              low_node,blockId2,input_matrix,k_cols,k_rows,
                                                              label_temporal_cols,label_temporal_rows,lambda_cols,lambda_rows,false,max_number_blocks);
                    
                    label_temporal_rows[low_node]=blockId2;
                    
                }else {
                    
                    lambdas_inblock_change_node_partition_row(
                                                              low_node,blockId2,input_matrix,k_cols,k_rows,
                                                              label_temporal_cols,label_temporal_rows,lambda_cols,lambda_rows,true,max_number_blocks);
                    lambdas_inblock_change_node_partition_row(
                                                              low_node,blockId1,input_matrix,k_cols,k_rows,
                                                              label_temporal_cols,label_temporal_rows,lambda_cols,lambda_rows,false,max_number_blocks);
                    
                    label_temporal_rows[low_node]=blockId1;
                    
                }
                //sprintf(outputString,"fprintf('end updating lambdas 1\\n');");
                //mexEvalString(outputString);
            }else{
                int low_node_col = low_node-((int)label_temporal_rows.size());
                
                //sprintf(outputString,"fprintf('end updating lambdas 1 %i = %i - %i \\n');",low_node_col,low_node,(int)label_temporal_rows.size());
                //mexEvalString(outputString);
                if (label_temporal_cols[low_node_col]==blockId1){
                    
                    lambdas_inblock_change_node_partition_col(
                                                              low_node_col,blockId1,input_matrix,k_cols,k_rows,
                                                              label_temporal_cols,label_temporal_rows,lambda_cols,lambda_rows,true,max_number_blocks);
                    lambdas_inblock_change_node_partition_col(
                                                              low_node_col,blockId2,input_matrix,k_cols,k_rows,
                                                              label_temporal_cols,label_temporal_rows,lambda_cols,lambda_rows,false,max_number_blocks);
                    
                    label_temporal_cols[low_node_col]=blockId2;
                    
                }else {
                    
                    lambdas_inblock_change_node_partition_col(
                                                              low_node_col,blockId2,input_matrix,k_cols,k_rows,
                                                              label_temporal_cols,label_temporal_rows,lambda_cols,lambda_rows,true,max_number_blocks);
                    lambdas_inblock_change_node_partition_col(
                                                              low_node_col,blockId1,input_matrix,k_cols,k_rows,
                                                              label_temporal_cols,label_temporal_rows,lambda_cols,lambda_rows,false,max_number_blocks);
                    
                    label_temporal_cols[low_node_col]=blockId1;
                    
                }
                //sprintf(outputString,"fprintf('end updating lambdas 2\\n');");
                //mexEvalString(outputString);
            }
            
        }else{
            // changing the node label
            // test instead of the if, replace the operation by using booleans
            //consider low_node ranges from 0 to (N+M)-1
            if (low_node<label_temporal_rows.size()){
                if (label_temporal_rows[low_node]==blockId1){
                    label_temporal_rows[low_node]=blockId2;
                }else {
                    label_temporal_rows[low_node]=blockId1;
                }
            }else{
                if (label_temporal_cols[low_node-label_temporal_rows.size()]==blockId1){
                    label_temporal_cols[low_node-label_temporal_rows.size()]=blockId2;
                }else {
                    label_temporal_cols[low_node-label_temporal_rows.size()]=blockId1;
                }
            }
            
            lambda_i(
                     input_matrix,
                     k_cols,
                     k_rows,
                     label_temporal_cols,
                     label_temporal_rows,
                     lambda_cols,
                     lambda_rows,
                     max_number_blocks,
                     ibn);
        }
        
        //sprintf(outputString,"fprintf('inside current_block->first_step->lambda_i_2 \\n');");
        //mexEvalString(outputString);
        
        In = calculate_Fitness(
                               input_matrix,
                               k_cols,
                               k_rows,
                               lambda_cols,
                               lambda_rows,
                               ibn);
        
        //sprintf(outputString,"fprintf('inside current_block->first_step->calculate_fitness \\n');");
        //mexEvalString(outputString);
        
        int suma_col_1=0,suma_col_2=0,suma_row_1=0,suma_row_2=0;
        for(int i=0; i<label_temporal_cols.size();i++){
            if (label_temporal_cols[i]==blockId1)
                suma_col_1++;
            else
                suma_col_2++;
        }
        for(int k=0; k<label_temporal_rows.size();k++){
            if (label_temporal_rows[k]==blockId1)
                suma_row_1++;
            else
                suma_row_2++;
        }
        
        
        // if an increase is found, save the changes

        if (In>Io){
            Io=In;
            label_cols = label_temporal_cols;
            label_rows = label_temporal_rows;
            
            alpha=0;
        }else{
            lambda_cols = lambda_cols_copy; //this and the restoration of the copy of label cols has the same
            lambda_rows = lambda_rows_copy; //objective, just thath the variable names operate differnt, the labels_temporal_cols logic
            //shold be improved
            alpha++;
        }
        if(max_alpha<alpha){
            max_alpha = alpha;
           // sprintf(outputString,"fprintf('%s: max_alpha %i of Malpha %i \\n');",log_prefix_global,max_alpha,Malpha);
           // mexEvalString(outputString);
        }
    }
    label_cols_new = label_cols;
    label_rows_new = label_rows;
}

vector<vector<int> > recursive_step(
                                    vector<vector<int> > & adjacency_matrix,
                                    vector<int> & k_cols,
                                    vector<int> & k_rows,
                                    double alpha_scale_factor,
                                    int repetitions,
                                    bool ibn){
  
    vector<vector<int> >  out;
    int N_cols=adjacency_matrix[0].size();
    int N_rows=adjacency_matrix.size();
    //int Malpha=int((N_cols+N_rows)*alpha_scale_factor);
    
    vector<int> col_labels_out(N_cols,1);
    vector<int> row_labels_out(N_rows,1);
    vector<int > labels_temp_col(N_cols,1);
    vector<int > labels_temp_row(N_rows,1);
    vector<double > lambda_cols(N_cols,0);
    vector<double > lambda_rows(N_rows,0);
    
    double If=-DBL_MAX;
    double I;
    int current_block=0;
    for (int i=0; i<repetitions; i++){
        
        //sprintf(outputString,"fprintf('Begin recursive_step function and calling other functions %i repetitions\\n');",i);
        //mexEvalString(outputString);
        //sprintf(outputString,"fprintf('%s: Initializing 1st block\\n');", log_prefix_global);
        //mexEvalString(outputString);
        
        double In=-DBL_MAX;
        
        vector<int> labels_final_col(N_cols,1);
        vector<int> labels_final_row(N_rows,1);
        
        int max_number_blocks_cols=*max_element(labels_final_col.begin(), labels_final_col.end());
        int max_number_blocks_rows=*max_element(labels_final_row.begin(), labels_final_row.end());
        int max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows) + 1;
        
        lambda_i(adjacency_matrix,
                 k_cols,
                 k_rows,
                 labels_final_col,
                 labels_final_row,
                 lambda_cols,
                 lambda_rows,
                 max_number_blocks,
                 ibn);
        
        I = calculate_Fitness(
                              adjacency_matrix,
                              k_cols,
                              k_rows,
                              lambda_cols,
                              lambda_rows,
                              ibn);
        
        //sprintf(outputString,"fprintf('%s: Initial fitness = %f\\n');",log_prefix_global,I);
        //mexEvalString(outputString);
        
        
        //sprintf(outputString,"fprintf('calculate_Fitness %f\\n');",I);
        //mexEvalString(outputString);
        
        //sprintf(outputString,"fprintf('current_block\\n');");
        //mexEvalString(outputString);
        
        //while (current_block <= *max_element(labels_final_col.begin(), labels_final_col.end())){
        while (current_block <= max_number_blocks){
            
           // sprintf(outputString,"fprintf('%s: splitting block = %i\\n');",log_prefix_global,current_block);
            //mexEvalString(outputString);
            
            //Building the new matrix from one of the two partitions
            int number_of_blocks=max_number_blocks;
            
            vector<int> indices_col;
            vector<int> indices_row;
            
            //sprintf(outputString,"fprintf('inside current_block->1 for\\n');");
            //mexEvalString(outputString);
            
            for (int i=0; i<labels_final_col.size(); i++){
                if (labels_final_col[i]==current_block){
                    indices_col.push_back(i);
                }
            }
            
            //sprintf(outputString,"fprintf('inside current_block->1 for\\n');");
            //mexEvalString(outputString);
            
            for (int l=0; l<labels_final_row.size(); l++){
                if (labels_final_row[l]==current_block){
                    indices_row.push_back(l);
                }
            }
            
            //sprintf(outputString,"fprintf('indices_row.size() = %lu && indices_col.size() = %lu\\n');",indices_row.size(),indices_col.size());
            //mexEvalString(outputString);
            
            if(indices_row.size()>1 && indices_col.size()>1){
                
                //keep a copy of old labels
                vector<int> old_labels_final_col = labels_final_col;
                vector<int> old_labels_final_row = labels_final_row;
                
                vector<vector<int> > sub_matrix; //new matrix
                vector<int> sub_matrix_k_cols(indices_col.size(),0);
                vector<int> sub_matrix_k_rows(indices_row.size(),0);
                for (int j=0; j<indices_row.size(); j++){
                    vector<int> aux;
                    for (int k=0; k<indices_col.size(); k++){
                        aux.push_back(adjacency_matrix[indices_row[j]][indices_col[k]]);
                    }
                    sub_matrix.push_back (aux);
                }
                //calculate degress
                for (int i_rows = 0; i_rows < indices_row.size(); i_rows++) {
                    for (int j_cols = 0; j_cols < indices_col.size(); j_cols ++) {
                        sub_matrix_k_rows[i_rows]+=sub_matrix[i_rows][j_cols];
                        sub_matrix_k_cols[j_cols]+=sub_matrix[i_rows][j_cols];
                    }
                }
                
               // sprintf(outputString,"fprintf('%s: Splitting block %i in two and number_of_blocks = %i\\n');",log_prefix_global,current_block,number_of_blocks);
               // mexEvalString(outputString);
               // sprintf(outputString,"fprintf('%s: sub_matrix size (%lu,%lu)\\n');",log_prefix_global,indices_row.size(),indices_col.size());
               // mexEvalString(outputString);
                
                //applying the first step function to the new matrix
                int n_row=int(sub_matrix.size());
                int n_col=int(sub_matrix[0].size());
                int malpha2=int((n_col+n_row)*alpha_scale_factor);
                
                //sprintf(outputString,
                //        "fprintf('inside current_block->bipartition = %i, n_col = %i, n_row = %i, alpha_scale_factor = %f\\n');",
                //        malpha2,n_col,n_row,alpha_scale_factor);
                //mexEvalString(outputString);
                
                int newBlockId = number_of_blocks + 1;
                vector<int> labels_cols(indices_col.size(),0);
                vector<int> labels_rows(indices_row.size(),0);
                
                bipartition(sub_matrix,labels_cols,labels_rows,newBlockId);
                max_number_blocks = newBlockId + 1 + 1;
                
                //sprintf(outputString,"fprintf('inside current_block->bipartition\\n');");
                //mexEvalString(outputString);
                
                first_step(
                           sub_matrix,
                           sub_matrix_k_cols,
                           sub_matrix_k_rows,
                           labels_cols,
                           labels_rows,
                           malpha2,
                           labels_cols,
                           labels_rows,
                           newBlockId,
                           newBlockId + 1,
                           max_number_blocks,
                           ibn);
                
                //sprintf(outputString,"fprintf('inside current_block->first_step\\n');");
                //mexEvalString(outputString);
                
                update_partitions_labels(
                                         labels_final_col,
                                         labels_final_row,
                                         labels_cols,
                                         labels_rows,
                                         current_block);
                
                //sprintf(outputString,"fprintf('inside current_block->update_partitions_labels\\n');");
                //mexEvalString(outputString);
                
                lambda_i(adjacency_matrix,
                         k_cols,
                         k_rows,
                         labels_final_col,
                         labels_final_row,
                         lambda_cols,
                         lambda_rows,
                         max_number_blocks,
                         ibn);
                
                //sprintf(outputString,"fprintf('inside current_block->lambda_i\\n');");
                //mexEvalString(outputString);
                
                In = calculate_Fitness(
                                       adjacency_matrix,
                                       k_cols,
                                       k_rows,
                                       lambda_cols,
                                       lambda_rows,ibn);
                
             //   sprintf(outputString,"fprintf('%s: New fitness = %f, max fitness %f\\n');",log_prefix_global,(double)In,(double)I);
               // mexEvalString(outputString);
                
                if (In>I){
                    I=In;
                  //  sprintf(outputString,"fprintf('%s: Bipartition accepted, new fitness improved = %f\\n');",log_prefix_global,(double)I);
                   // mexEvalString(outputString);
                    
                }else{
                    labels_final_col = old_labels_final_col;
                    labels_final_row = old_labels_final_row;
                    
                  //  sprintf(outputString,"fprintf('%s: Rejecting bi-partition\\n');",log_prefix_global);
                   // mexEvalString(outputString);
                }
            }
            
            current_block++;
            
            //recompute the maximum number of blocks
            max_number_blocks_cols=*max_element(labels_final_col.begin(), labels_final_col.end());
            max_number_blocks_rows=*max_element(labels_final_row.begin(), labels_final_row.end());
            max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows) + 1;
            
        }
        
        if (I > If){
            If = I;
            col_labels_out = labels_final_col;
            row_labels_out = labels_final_row;
            
        }
        
    }
    
    int min_block;
    int min_block_col=*min_element(col_labels_out.begin(), col_labels_out.end());
    int min_block_row=*min_element(row_labels_out.begin(), row_labels_out.end());
    if (min_block_col < min_block_row){
        min_block = min_block_col;
    }else{
        min_block = min_block_row;
    }
    
    for (int k=0; k<col_labels_out.size();k++){
        col_labels_out[k]=(col_labels_out[k]-min_block);
    }
    for(int m=0; m<row_labels_out.size();m++){
        row_labels_out[m]=(row_labels_out[m]-min_block);
    }
    
    out.push_back(row_labels_out);
    out.push_back(col_labels_out);
    
    return out;
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

//data_out_EO extremal_optimization(const mxArray *is, double alpha_parameter, int repetitions, bool ibn,char * log_prefix){
//    data_out_EO var;
//    vector<vector<int> > M;
//    ifstream is(Filename);
//    load_matrix(&is, &M);
//    //log_prefix_global = log_prefix;
//    //load_matrix_sparseMatlab(is, M);
//
//    //network degress by rows and columns
//    int N_cols=M[0].size();
//    int N_rows=M.size();
//    vector<int> k_rows(N_rows,0);
//    vector<int> k_cols(N_cols,0);
//    for (int i = 0; i < N_rows; ++i) {
//        for (int j = 0; j < N_cols; ++j) {
//            k_rows[i]+=M[i][j];
//            k_cols[j]+=M[i][j];
//        }
//    }
//
//    vector<vector<int> > partitions = recursive_step(
//                                                     M,
//                                                     k_cols,
//                                                     k_rows,
//                                                     alpha_parameter,
//                                                     repetitions, ibn);
//
//    int max_number_blocks_cols=*max_element(partitions[1].begin(), partitions[1].end());
//    int max_number_blocks_rows=*max_element(partitions[0].begin(), partitions[0].end());
//    int max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows) + 1;
//
//    vector<vector<double> > lambdas = call_lambda_i(
//                                                    M,
//                                                    k_cols,
//                                                    k_rows,
//                                                    partitions[1],
//                                                    partitions[0],
//                                                    max_number_blocks,
//                                                    ibn);
//
//    double Q=calculate_Fitness(
//                               M,
//                               k_cols,
//                               k_rows,
//                               lambdas[0],
//                               lambdas[1],
//                               ibn);
//
//    var.QI_eo= Q;
//    var.partitions_result=partitions;
//    return var;
//}

PYBIND11_MODULE(extremal_bi,m) {
    //py::module m("example", "Generating primes in c++ with python bindings using pybind11");
    m.def("recursive_step", &recursive_step, "Extremal optimization algorithm to detect communities");
    m.def("call_lambda_i", &call_lambda_i, "A function that calculates the fitness contribution of each node");
    m.def("calculate_Fitness", &calculate_Fitness, "A function that calculate the Ieo and Qeo of the whole network");
    m.def("bipartition", &bipartition, "Function to perform the bipartition of a given array");
    m.def("shuffle_partitions", &shuffle_partition, "Function to perform the shiffle on the bipartition of a given array");
    m.def("first_step", &first_step, "Function to move the node with lowest fitness until optimal Q/I is reached");
    // return m.ptr();
}
// compilation on linux g++ -O3 -Wall -shared -std=c++11 -fPIC `python -m pybind11 --includes` EO_functions.cpp -o extremal_bi.so

//compilation on mac: g++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python -m pybind11 --includes` EO_functions.cpp -o extremal_bi.so
