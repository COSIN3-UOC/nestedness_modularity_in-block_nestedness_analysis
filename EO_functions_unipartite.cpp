/*
 * IBN_bi.cpp
 *
 *  Created on: 12 feb. 2018
 *      Author: mariapalazzi
 */
////mex -largeArrayDims IBN_bi.cpp
        
//#include <//mex.h>
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
#include <cfloat>
#include <cstdio>
#include <numeric>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/functional.h>
#include <pybind11/stl_bind.h>
//#include "EO_aux_functions_opt1.hpp"

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) < (Y)) ? (Y) : (X))

using namespace std;
//void log_file( const string &text )
//{
//    if (false) {
//        ofstream log_file("log_file.log", ios_base::out | ios_base::app );
//        log_file << text << endl;
//    }
//}

//char outputString[500];
//char * log_prefix_global;

vector<int> myargsort(vector<double> x){
    int la = x.size();
    vector<double> c = x;
    vector<int> y (la);
    int i=0,ii,ti;
    for(i=0;i<la;i++)
        y[i]=i;
    i=0;
    double td;
    while (i<(la-1)){
        ii=i+1;
        while (ii<(la)){
            if (c[ii]<c[i]){
                td = c[i];
                c[i] = c[ii];
                c[ii] = td;
                ti = y[i];
                y[i] = y[ii];
                y[ii] = ti;
            }
            ii++;
        }
        i++;
    }
    return y;
}

vector<int> indexes_sorted_vector(vector<double> & x){
    
    std::vector<int> vs(x.size());
    std::iota(vs.begin(), vs.end(), 0);    
    auto comparator = [&x](int a, int b){ return x[a] < x[b]; };
    std::sort(vs.begin(), vs.end(), comparator);    
    //for (int i=0; i<5; i++){
    //    cout<< vs[i] <<endl;
    //}    
    
    return vs;
}

void shuffle_partitions(vector<int> & labels){

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
                
                shuffle_partitions(label_col);
            }
        }else{
            label_col[0] = nextBlockId;
            label_col[1] = nextBlockId + 1;
        }
        
//     int N_rows=input_matrix.size();
//     int N_cols=input_matrix[0].size();
//     vector<int> label_row2 (N_rows,0);
//     vector<int> label_col2 (N_cols,0);
//     /* Splitting in two and assigning the labels*/
//     for (unsigned int i=0; i<(N_rows);i++){
//         if (i<(N_rows/2))
//             label_row2[i]=1;
//         else
//             label_row2[i]=2;
//     }
//     for (unsigned int i=0; i<(N_cols);i++){
//         if (i<N_cols/2)
//             label_col2[i]=1;
//         else
//             label_col2[i]=2;
//     }
//     /* Shuffling the labels*/
//     random_device rd; //non fixed seed
//     mt19937 gen(rd());
//     shuffle ( label_col2.begin(), label_col2.end(), gen);
//     shuffle ( label_row2.begin(), label_row2.end(), gen);
// 
//     label_row=label_row2;
//     label_col=label_col2;
}

void lambdas_inblock(
        vector<vector<int> > & input_matrix,
        vector<int> & k_cols,
		vector<int> label_cols,
		vector<double>& lambda_cols,
        int max_number_blocks){
    
    int N_cols=input_matrix[0].size();
    int N_rows=input_matrix.size();
    lambda_cols= vector<double>(N_cols,0);
//    vector<int>::iterator it;
    //vector<int> k_cols(N_cols,0);
    //vector<int> k_rows(N_rows,0);

    // getting current number of blocks
    //int max_number_blocks_cols=*max_element(label_cols.begin(), label_cols.end());
    //int max_number_blocks_rows=*max_element(label_rows.begin(), label_rows.end());
    //int max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows);
    
    vector<vector<int> > list_of_blocks_cols(max_number_blocks+1,vector<int>());

    //sprintf(outputString,"fprintf('inside current_block->first_step->lambdas_inblock_rows->start (%lu,%lu), max_blocks %i\\n');",
    //         label_rows.size(),
    //         label_cols.size(),
    //        max_number_blocks);
    //mexEvalString(outputString);
    
    //creating a vector that separate the nodes according to the block the belong to columns
    //     for (int i=0; i<= max_number_blocks; i++){
    //         vector<int> appending_cols_nodes;
    //         for (int j = 0; j < label_cols.size(); ++j){
    //             if (label_cols[j]==i){
    //                 appending_cols_nodes.push_back(j);
    //             }
    //         }
    //         list_of_blocks_cols.push_back(appending_cols_nodes);
    //     }
    for (int j = 0; j < label_cols.size(); ++j){
        list_of_blocks_cols[label_cols[j]].push_back(j);
    }

    //sprintf(outputString,
    //        "fprintf('inside current_block->first_step->lambdas_inblock_rows->list_of_blocks_cols size %lu\\n');",
    //        label_cols.size());
    //mexEvalString(outputString);
    
    //rows
    //     for (unsigned int i=0; i<= max_number_blocks; i++){
    //             		vector<int> appending_rows_nodes;
    //             		for (unsigned int j = 0; j < label_rows.size(); ++j){
    //             			if (label_rows[j]==i){
    //             				appending_rows_nodes.push_back(j);
    //             			}
    //             		}
    //             		list_of_blocks_rows.push_back(appending_rows_nodes);
    //             }

    //sprintf(outputString,
    //        "fprintf('inside current_block->first_step->lambdas_inblock_rows->list_of_blocks_rows size %lu\\n');",
    //        label_rows.size());
    //mexEvalString(outputString);
    
    //getting the k degrees <--- this is now an input parameter, this does not change in all the iterations
    // for (unsigned int i = 0; i < label_rows.size(); ++i) {
    //     for (unsigned int j = 0; j < label_cols.size(); ++j) {
    //         k_rows[i]+=input_matrix[i][j];
    //         k_cols[j]+=input_matrix[i][j];
    //     }
    // }

     //computing the column nodes contribution (pair overlap)
    for (int l=0; l<list_of_blocks_cols.size(); l++){
        int size_blocks_cols=(double)list_of_blocks_cols[l].size();
        
        for (int i=0; i<size_blocks_cols; i++){
            int ii = list_of_blocks_cols[l][i];
            
            for (int j=0; j<size_blocks_cols; j++){
                
                int jj =list_of_blocks_cols[l][j];
                double PO_col_i=0;
                double normalization = ((double)(k_cols[jj]))*((double)(size_blocks_cols-1.0));//double normalization = (2/(N_rows+N_cols))*(1/((double)(k_cols[jj])*(double)(size_blocks_cols-1.0)));

                if ( (k_cols[ii]>=k_cols[jj])&(k_cols[jj]>0) & (normalization!=0.0)&(ii!=jj)){
                		for (int k=0; k<list_of_blocks_cols[l].size(); k++){
						int kk = list_of_blocks_cols[l][k];
						if ((input_matrix[kk][ii]*input_matrix[kk][jj])==1){
								PO_col_i++;
						}
					 }
                		if (k_cols[ii]==k_cols[jj]){
                    double null_model = ((double)(k_cols[ii]*k_cols[jj]))/((double)(N_rows));
                    PO_col_i=(PO_col_i-null_model)/(2.0*normalization);
                    lambda_cols[ii]+=PO_col_i;
                		} else{
                			double null_model = ((double)(k_cols[ii]*k_cols[jj]))/((double)(N_rows));
                			PO_col_i=(PO_col_i-null_model)/normalization;
                			lambda_cols[ii]+=PO_col_i;
                		}
                }
            }
        }
    }

    //sprintf(outputString,"fprintf('inside current_block->first_step->lambdas_inblock_cols \\n');");
    //mexEvalString(outputString);

    //sprintf(outputString,"fprintf('inside current_block->first_step->lambdas_inblock_rows \\n');");
    //mexEvalString(outputString);
    
}

void lambdas_modularity(
        vector<vector<int> > & input_matrix,
        vector<int> & k_cols,
		vector<int> label_cols, 
		vector<double>& lambda_cols, 
        int max_number_blocks){
    
        int N_cols=input_matrix[0].size();
        vector<int>::iterator it;
        vector<int> blocks_cols(label_cols);
        lambda_cols=vector<double>(N_cols,0);
        vector<double> kappa_cols(N_cols,0);
        //vector<int> k_rows(N_rows,0);
        //vector<int> k_cols(N_cols,0);

        //obtaining current number of blocks
        //cols
        sort(blocks_cols.begin(), blocks_cols.end());
        it=unique(blocks_cols.begin(), blocks_cols.end());
        blocks_cols.resize(distance(blocks_cols.begin(),it) );
        //vector<double> links_blocks_cols(*max_element(blocks_cols.begin(), blocks_cols.end()),0);
        vector<double> links_blocks_cols(max_number_blocks,0);
        double total_links=0.;

        //getting the total links
         for (unsigned int i = 0; i < label_cols.size(); ++i) {
            for (unsigned int j = 0; j < label_cols.size(); ++j) {
                total_links+=input_matrix[i][j];
            }
         }
//        total_links = accumulate(k_cols.begin(), k_cols.end(), 0.0);
        
       // getting the kappa's
       for (unsigned int i = 0; i < label_cols.size(); ++i) {
           for (unsigned int j = 0; j < label_cols.size(); ++j) {
               if  ((input_matrix[i][j]==1) & (label_cols[i]==label_cols[j])){
                   kappa_cols[j]+=1./(double)k_cols[j];
               }
           }
       }
       // getting the a_r(i)'s

       for (unsigned int i = 0; i < label_cols.size(); ++i) {
    	   links_blocks_cols[label_cols[i]-1]+=k_cols[i]/total_links;
       }
       //getting the lambda fitness
       for (unsigned int i = 0; i < label_cols.size(); ++i) {
           lambda_cols[i]=kappa_cols[i] - links_blocks_cols[label_cols[i]-1];
       }
}

double calculate_inblock_nestedness(
        vector<vector<int> > & input_matrix,
		vector<double> lambda_cols){

	double I;
	int N_cols=input_matrix[0].size();
	double i_col=0;

	//for (unsigned int i=0; i<N_cols; i++){
	//	i_col+=lambda_cols[i];
	//}
    i_col = accumulate(lambda_cols.begin(), lambda_cols.end(), 0.0);

	//for (unsigned int j=0; j<N_rows; j++){
	//	i_row+=lambda_rows[j];
	//}
	I=(2.0/((double)N_cols))*(i_col);
        //log_file("\tleaving calculate_inblock_nestedness function");
        //log_file("\t-- I: "+ to_string(I));
        return I;
}

double calculate_modularity(
        vector<vector<int> > & input_matrix,
        vector<int> & k_cols,
		vector<double> lambda_cols){

	double Q;
	double q_cols=0;
	int N_cols=input_matrix[0].size();
	//vector<int> k_rows(N_rows,0);
	//vector<int> k_cols(N_cols,0);
	double total_links=0;
	for (unsigned int i = 0; i < N_cols ; ++i) {
		for ( int j = 0; j < N_cols; ++j) {
			total_links+=input_matrix[i][j];
		}
	}
//    total_links = accumulate(k_rows.begin(), k_rows.end(), 0.0);

	for (int i=0; i<N_cols; i++){
		q_cols+=k_cols[i]*lambda_cols[i];
	}
	Q=(q_cols)/(total_links);
	//    cout<<Q<<endl;
        //log_file("\tleaving calculate_modularity function");
	return Q;

}

void lambda_i(
        vector<vector<int> > & input_matrix,
        vector<int> & k_cols,
		vector<int> label_cols, 
		vector<double> &lambda_cols, 
        int max_number_blocks,
        bool ibn){

    if (ibn==true){
        /* Function that calculate the fitness contributions of nodes (lambda_i) for the in-block nestedness*/
        lambdas_inblock(
            input_matrix,
            k_cols,
            label_cols,
            lambda_cols,
            max_number_blocks);
    }else {
        /* Function that calculate the fitness contributions of nodes (lambda_i) for the modularity*/
        lambdas_modularity(
            input_matrix,
            k_cols,
            label_cols,
            lambda_cols,
            max_number_blocks);
    }
    //log_file("\tleaving lambda_i function");
}
vector<double>  call_lambda_i(
        vector<vector<int> > & input_matrix,
        vector<int> & k_cols,
		vector<int> label_cols,
        int max_number_blocks,
		bool ibn){
	vector<double> out;
	vector<double> lambda_cols(label_cols.size(),0);
	lambda_i(
            input_matrix, 
            k_cols,
            label_cols,
            lambda_cols,
            max_number_blocks,
            ibn);
	out=lambda_cols;
	return out;
}

double calculate_Fitness(
        vector<vector<int> > & input_matrix,
        vector<int> & k_cols,
		vector<double> lambda_cols,
        bool ibn){
    //log_file("\tentering calculate_Fitness function");
    /*  Calculate the EO metric of the whole network */
	double metric;
    if (ibn==true) {
    	// in block nestedness
    	metric=calculate_inblock_nestedness(input_matrix, lambda_cols);
    }else{ // modularity 
    	metric=calculate_modularity(
                input_matrix, 
                k_cols,
                lambda_cols);
    }   

//     double (* calculate_genericFitness)(vector<vector<int> > &,vector<int> &,vector<int> &,vector<double>,vector<double>);
//     
//     metric=calculate_genericFitness(
//                 input_matrix, 
//                 k_cols,
//                 k_rows,                       
//                 lambda_cols,
//                 lambda_rows);
            
    //log_file("\t++ Metric: "+ to_string(metric));
    //log_file("\tleaving calculate_Fitness function");
    return metric;
}


/*Obtain the node that is going to be moved from its block, based on
 probability distribution. tau-EO like in Duch et al 2005.*/
int low_fitness_node(
        vector<double> lambda_cols) {

    int low_node;
    int N=lambda_cols.size();
    double tau = 1.+1./log(N); //tau-exponent
    vector<double> probabilities(N,0);
    double p_norm=0;
    vector<int> lambda_sorted(N);

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
    lambda_sorted = indexes_sorted_vector(lambda_cols);
    low_node=lambda_sorted[distribution(gen)];

    //log_file("\tleaving low_fitness_node function");
    return low_node;
}

void first_step(
                vector<vector<int> > &input_matrix,
                vector<int> & k_cols,
                vector<int> label_cols,
                int Malpha,
                vector<int> &label_cols_new,
                int blockId1,
                int blockId2,
                int max_number_blocks,
                bool ibn){
    /*Given a matrix arbitrarily partitioned in two, move the nodes with "lowest" fitness
     until optimal Ieo is reached. The stopping criteria for the optimal Ieo is given by
     Malpha.*/
    //log_file("****entering first_step function");
    vector<double> lambda_cols(label_cols.size(),0);
    vector<double> lambda_cols_copy(label_cols.size(),0);

    vector<int> label_temporal_cols;
    double Io=-1.0;
    double In=-1.0;
    int alpha=0;
    int max_alpha=0;
    int low_node=0;

    //moving the nodes of a current partition from one block to the other
    //    int c=0;
    //if(Malpha > 50){
    //    mexErrMsgTxt("Malpha error inside recursive_step\n");
    //}
    
    //compute the first lambda_i
    label_temporal_cols = label_cols;
    
    //sprintf(outputString,"fprintf('inside current_block->first_step->alpha %i, Malpha %i, max_blocks %i\\n');",
    //        alpha,Malpha,max_number_blocks);
    //mexEvalString(outputString);
    
    lambda_i(
             input_matrix,
             k_cols,
             label_temporal_cols,
             lambda_cols,
             max_number_blocks,
             ibn);
    
    while (alpha<Malpha){
        
        label_temporal_cols = label_cols;

        lambda_cols_copy = lambda_cols;

        //sprintf(outputString,"fprintf('inside current_block->first_step->alpha %i, Malpha %i \\n');",alpha,Malpha);
        //mexEvalString(outputString);
        
        //sprintf(outputString,"fprintf('inside current_block->first_step->lambda_i \\n');");
        //mexEvalString(outputString);
        
        low_node=low_fitness_node(lambda_cols);
        
        //sprintf(outputString,"fprintf('inside current_block->first_step->low_fitness_node \\n');");
        //mexEvalString(outputString);
        
        //TODO: do this update of lambda_i's for IBN in a nicer way
            // changing the node label
            // test instead of the if, replace the operation by using booleans
            //consider low_node ranges from 0 to (N+M)-1
		if (label_temporal_cols[low_node]==blockId1){
				label_temporal_cols[low_node]=blockId2;
			}else {
				label_temporal_cols[low_node]=blockId1;
			}

		lambda_i(
				 input_matrix,
				 k_cols,
				 label_temporal_cols,
				 lambda_cols,
				 max_number_blocks,
				 ibn);

        
        //sprintf(outputString,"fprintf('inside current_block->first_step->lambda_i_2 \\n');");
        //mexEvalString(outputString);
        
        In = calculate_Fitness(
                               input_matrix,
                               k_cols,
                               lambda_cols,
                               ibn);
        
        //sprintf(outputString,"fprintf('inside current_block->first_step->calculate_fitness \\n');");
        //mexEvalString(outputString);
        
        int suma_col_1=0,suma_col_2=0;
        for(int i=0; i<label_temporal_cols.size();i++){
            if (label_temporal_cols[i]==blockId1)
                suma_col_1++;
            else
                suma_col_2++;
        }

        //sprintf(outputString,"fprintf('inside current_block->first_step->suma_row_1 (%i,%i,%i,%i) \\n');",
        //        suma_col_1,
        //        suma_col_2,
        //        suma_row_1,
        //        suma_row_2);
        //mexEvalString(outputString);
        
        // if an increase is found, save the changes
        //if ((In>Io) & (suma_col_1>1) & (suma_col_2>1)& (suma_row_1>1) & (suma_row_2>1)){
        if (In>Io){
            Io=In;
            label_cols = label_temporal_cols;

            alpha=0;
        }else{
            lambda_cols = lambda_cols_copy; //this and the load of the copy of label cols has the same
            //objective, just thath the variable names operate differnt, the labels_temporal_cols logic
            //shold be improved
            alpha++;
        }
        if(max_alpha<alpha){
            max_alpha = alpha;
//            sprintf(outputString,"fprintf('%s: max_alpha %i of Malpha %i \\n');",log_prefix_global,max_alpha,Malpha);
//            mexEvalString(outputString);
        }
    }
    label_cols_new = label_cols;
    //log_file("**** ** New value of In: " + to_string(In));
}

void update_partitions_labels(
        vector<int> &total_in_label_cols, 
		vector<int> label_partition_cols,
		int current_block){
        //log_file("\tentering update_partitions_labels function");
	int j=0;
	//updating the label vector with the new partitions
	for (unsigned int i=0; i<total_in_label_cols.size(); i++){
		if (total_in_label_cols[i]==current_block){
			total_in_label_cols[i]=label_partition_cols[j];
			j++;
		}
	}
        //log_file("\tleaving update_partitions_labels function");
}

//typedef struct partition_component {
//    vector<int> row_ids;
//    vector<int> col_ids;
//    int pos_partition_rows;
//    int pos_partition_cols;
//} partition_component;

vector<int>  recursive_step(
        vector<vector<int> > & adjacency_matrix, 
        vector<int> & k_cols,
        double alpha_scale_factor, 
        int repetitions,
		bool ibn){
    //log_file("\tentering recursive_step function");
    vector<vector<int> >  out;
    int N_cols=adjacency_matrix[0].size();
//    int Malpha=int((N_cols+N_rows)*alpha_scale_factor);
    
    vector<int> col_labels_out(N_cols,1);
    vector<int > labels_temp_col(N_cols,1);
    vector<double > lambda_cols(N_cols,0);
            
    //printf("Begin recursive_step function and calling other functions %i repetitions\n",repetitions);
    ////sprintf(outputString,"fprintf('Begin recursive_step function and calling other functions %i repetitions\\n');",repetitions);
    ////mexEvalString(outputString);
    
    //if(Malpha > 50){
    //    //mexErrMsgTxt("Malpha error recursive_step\n");
    //}        
    //std::deque<partition_component> partitions_to_explore();
    //vector<vector<int>>current_parition;
    
    double If=-DBL_MAX;
    double I;
    int current_block=1;    
    for (int i=0; i<repetitions; i++){
        
        ////sprintf(outputString,"fprintf('Begin recursive_step function and calling other functions %i repetitions\\n');",i);
        ////mexEvalString(outputString);
        //sprintf(outputString,"fprintf('%s: Initializing 1st block\\n');", log_prefix_global);
//        //mexEvalString(outputString);
        
        double In=-DBL_MAX;

        vector<int> labels_final_col(N_cols,1);
        
        int max_number_blocks_cols=*max_element(labels_final_col.begin(), labels_final_col.end());
        int max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_cols) + 1;
                        
        lambda_i(adjacency_matrix,
                k_cols,
                labels_final_col,
                lambda_cols,
                max_number_blocks,
                ibn);

        I = calculate_Fitness(
                adjacency_matrix,
                k_cols,
                lambda_cols,
                ibn);
//        col_labels_out = labels_final_col;
//        row_labels_out = labels_final_row;
                        
        //while (current_block <= *max_element(labels_final_col.begin(), labels_final_col.end())){
        while (current_block <= max_number_blocks){
            
            //Building the new matrix from one of the two partitions
            //int number_of_blocks=*max_element(labels_final_col.begin(), labels_final_col.end());
            int number_of_blocks=max_number_blocks;
                        
            vector<int> indices_col;

            for (int i=0; i<labels_final_col.size(); i++){
                if (labels_final_col[i]==current_block){
                    indices_col.push_back(i);
                }
            }

            if( indices_col.size()>1){
                
                //keep a copy of old labels
                vector<int> old_labels_final_col = labels_final_col;
                vector<vector<int> > sub_matrix; //new matrix
                vector<int> sub_matrix_k_cols(indices_col.size(),0);
                for (int j=0; j<indices_col.size(); j++){
                    vector<int> aux;
                    for (int k=0; k<indices_col.size(); k++){
                        aux.push_back(adjacency_matrix[indices_col[j]][indices_col[k]]);
                    }
                    sub_matrix.push_back (aux);
                }
                //calculate degress
                for (int i_cols = 0; i_cols < indices_col.size(); i_cols++) {
                    for (int j_cols = 0; j_cols < indices_col.size(); j_cols ++) {
                         sub_matrix_k_cols[j_cols]+=sub_matrix[i_cols][j_cols];
                   }
               }
                //applying the first step function to the new matrix

                int n_col=int(sub_matrix[0].size()); 
                int malpha2=int((n_col)*alpha_scale_factor);


                int newBlockId = number_of_blocks + 1;
                vector<int> labels_cols(indices_col.size(),0);
                
                bipartition(sub_matrix,labels_cols,newBlockId);
                max_number_blocks = newBlockId + 1 + 1;
                
                ////sprintf(outputString,"fprintf('inside current_block->bipartition\\n');");
                ////mexEvalString(outputString);
                
                //if(malpha2 > 50){
                //    //mexErrMsgTxt("Malpha error current_block->recursive_step\n");
                //}                        
                first_step(
                        sub_matrix, 
                        sub_matrix_k_cols,
                        labels_cols, 
                        malpha2,
                        labels_cols, 
                        newBlockId,
                        newBlockId + 1,
                        max_number_blocks,
                        ibn);
                                                
                ////sprintf(outputString,"fprintf('inside current_block->first_step\\n');");
                ////mexEvalString(outputString);

                update_partitions_labels(
                        labels_final_col,
                        labels_cols,
                        current_block);                
                                                
                ////sprintf(outputString,"fprintf('inside current_block->update_partitions_labels\\n');");
                ////mexEvalString(outputString);

                lambda_i(adjacency_matrix, 
                        k_cols,
                        labels_final_col,
                        lambda_cols,
                        max_number_blocks,
                        ibn);                
                
                ////sprintf(outputString,"fprintf('inside current_block->lambda_i\\n');");
                ////mexEvalString(outputString);
                
                In = calculate_Fitness(
                        adjacency_matrix,
                        k_cols,
                        lambda_cols,ibn);

                //sprintf(outputString,"fprintf('%s: New fitness = %f, max fitness %f\\n');",log_prefix_global,(double)In,(double)I);
                //mexEvalString(outputString);
             
                ////log_file("New value of In " + to_string(In) );
                if (In>I){
                    I=In;
                    //labels_temp_col = vector<int>(labels_final_col);
                    //labels_temp_row = vector<int>(labels_final_row);
                    //sprintf(outputString,"fprintf('%s: Bipartition accepted, new fitness improved = %f\\n');",log_prefix_global,(double)I);
                    //mexEvalString(outputString);
                                        
                }else{
                    labels_final_col = old_labels_final_col;
                    
//                    //sprintf(outputString,"fprintf('%s: Rejecting bi-partition\\n');",log_prefix_global);
//                    //mexEvalString(outputString);
                }                
            }
            
            current_block++;   
            
            //recompute the maximum number of blocks
            max_number_blocks_cols=*max_element(labels_final_col.begin(), labels_final_col.end());
            max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_cols) + 1;
            
        }
        
        if (I > If){
            If = I;
            col_labels_out = labels_final_col;
            //printf("new Io = %f\n",Io);
                        
        }        
    }
    int min_block;
    int min_block_col=*min_element(col_labels_out.begin(), col_labels_out.end());
    min_block = min_block_col;

    min_block=(min_block-1);
    for (int k=0; k<col_labels_out.size();k++){
       col_labels_out[k]=(col_labels_out[k]-min_block);
    }
    
    // int max_number_blocks_cols=*max_element(col_labels_out.begin(), col_labels_out.end());
    // int max_number_blocks_rows=*max_element(row_labels_out.begin(), row_labels_out.end());  
    // int max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows) + 1;    
    // lambda_i(adjacency_matrix, 
    //         k_cols,
    //         k_rows,                               
    //         col_labels_out,
    //         row_labels_out,
    //         lambda_cols,
    //         lambda_rows,
    //         max_number_blocks,
    //         ibn);                
    // double In = calculate_Fitness(
    //         adjacency_matrix,
    //         k_cols,
    //         k_rows,                                   
    //         lambda_cols,
    //         lambda_rows,ibn);                    
    // //sprintf(outputString,"fprintf('%s: final Q value = %f, max_number_blocks = %i\\n');",
    //         log_prefix_global,(double)In,max_number_blocks);
    // //mexEvalString(outputString);
    
    //log_file("leaving recursive_step function");
    return col_labels_out;
}

//void load_matrix_sparseMatlab(const mxArray* is,vector< vector<int> > & matrix){
//
//    double *Gpr = mxGetPr(is);
//    size_t *Gir = mxGetIr(is);
//    size_t *Gjc = mxGetJc(is);
//
//    size_t M = mxGetM(is);
//    size_t N = mxGetN(is);
//
//    //sprintf(outputString,"fprintf('%s: Reading matrix of size (%i,%i)\\n');",log_prefix_global,(int)M,(int)N);
//    //mexEvalString(outputString);
//
//    size_t sIndEdges;
//    size_t eIndEdges;
//
//    // clear first
//    matrix.clear();
//    matrix.resize(N,vector<int>(M,0));
//
//    int w;
//    for(int v=0; v<N; v++){ //<--- should this increment on 2
//        ////sprintf(outputString,"fprintf('Row %i\\n');",(int)v);
//        ////mexEvalString(outputString);
//
//        //matrix.push_back(vector<int>(M,0));
//
//        sIndEdges = Gjc[v];
//        eIndEdges = Gjc[v+1];
//        for( int nInd = sIndEdges ; nInd<eIndEdges ; nInd++ ){
//            w = Gir[nInd];
//            matrix[v][int(w)] = 1;
//            //printf("row %i, col %i\n",v,w);
//        }
//    }
//    //sprintf(outputString,"fprintf('%s: end read matrix\\n');",log_prefix_global);
//    //mexEvalString(outputString);
//}

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
//    log_file("File read: pass");
}

//data_out fitness_of_partition(const mxArray *is, vector<int> partition_rows, vector<int> partition_cols, int ibn,char * log_prefix){
//
//    data_out var;
//    vector<vector<int> > M;
//    //ifstream is(Filename);
//    //load_matrix(&is, &M);
//    log_prefix_global = log_prefix;
//    load_matrix_sparseMatlab(is, M);
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
//    vector<vector<int> > partitions;
//    partitions.push_back(partition_rows);
//    partitions.push_back(partition_cols);
//
//    int max_number_blocks_cols=*max_element(partitions[1].begin(), partitions[1].end());
//    int max_number_blocks_rows=*max_element(partitions[0].begin(), partitions[0].end());
//    int max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows) + 1;
//
//    vector<vector<double> > lambdas = call_lambda_i(
//            M,
//            k_cols,
//            k_rows,
//            partitions[1],
//            partitions[0],
//            max_number_blocks,
//            ibn);
//
//    double Q=calculate_Fitness(
//            M,
//            k_cols,
//            k_rows,
//            lambdas[0],
//            lambdas[1],
//            ibn);
//
//    var.QI_eo= Q;
//    var.partitions_result=partitions;
//    return var;
//}

//data_out extremal_optimization(const mxArray *is, double alpha_parameter, int repetitions, bool ibn,char * log_prefix){
//    data_out var;
//    vector<vector<int> > M;
//    //ifstream is(Filename);
//    //load_matrix(&is, &M);
//    log_prefix_global = log_prefix;
//    load_matrix_sparseMatlab(is, M);
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
//            M,
//            k_cols,
//            k_rows,
//            alpha_parameter,
//            repetitions, ibn);
//
//    int max_number_blocks_cols=*max_element(partitions[1].begin(), partitions[1].end());
//    int max_number_blocks_rows=*max_element(partitions[0].begin(), partitions[0].end());
//    int max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_rows) + 1;
//
//    vector<vector<double> > lambdas = call_lambda_i(
//            M,
//            k_cols,
//            k_rows,
//            partitions[1],
//            partitions[0],
//            max_number_blocks,
//            ibn);
//
//    double Q=calculate_Fitness(
//            M,
//            k_cols,
//            k_rows,
//            lambdas[0],
//            lambdas[1],
//            ibn);
//
//    var.QI_eo= Q;
//    var.partitions_result=partitions;
//    return var;
//}
//int main(){
////	data_out results;
////	results=extremal_optimization("ibnb2.dat",.015,60,true);
////        cout << results.QI_eo << endl ;
////        for (int i=0; i<results.partitions_result[0].size(); i++){
////               cout<< i <<' '<< results.partitions_result[0][i]<<' '<< results.partitions_result[1][i]<<endl;
////               }
//	vector<vector<int> > M;
//	string Filenme = "ibn.dat";
//	ifstream is(Filenme);
//	load_matrix(&is, &M);
//	bool ibn=true;
//	  //network degress by rows and columns
//	    int N_cols=M[0].size();
//	    vector<int> k_cols(N_cols,0);
//	    for (int i = 0; i < N_cols; ++i) {
//	        for (int j = 0; j < N_cols; ++j) {
//	            k_cols[j]+=M[i][j];
//	        }
//	    }
//	vector<int> partitions;
//	partitions=recursive_step(M, k_cols,1 ,5, ibn);
//	for (unsigned int i=0; i<partitions.size(); i++){
//			   cout<< i <<' '<< partitions[i]<<endl;
//			   }
//		   //    cout<<Q<<endl;
//
//	int max_number_blocks_cols=*max_element(partitions.begin(), partitions.end());
//	int max_number_blocks = MAX(max_number_blocks_cols,max_number_blocks_cols) + 1;
//	vector<double> lambdas;
//	lambdas=call_lambda_i(M,  k_cols, partitions,max_number_blocks,ibn);
//	double Qeo=calculate_Fitness(M,  k_cols, lambdas,ibn);
//	   cout << Qeo <<" 0"<< endl;
//	return 0;
//}


PYBIND11_MODULE(extremal_uni,m) {
    //py::module m("example", "Generating primes in c++ with python bindings using pybind11");
    m.def("recursive_step", &recursive_step, "Extremal optimization algorithm to detect communities");
    m.def("call_lambda_i", &call_lambda_i, "A function that calculates the fitness contribution of each node");
    m.def("calculate_Fitness", &calculate_Fitness, "A function that calculate the Ieo and Qeo of the whole network");
    m.def("bipartition", &bipartition, "Function to perform the bipartition of a given array");
    m.def("shuffle_partitions", &shuffle_partitions, "Function to perform the shiffle on the bipartition of a given array");
    m.def("first_step", &first_step, "Function to move the node with lowest fitness until optimal Q/I is reached");
    // return m.ptr();
}
// compilation on linux g++ -O3 -Wall -shared -std=c++11 -fPIC `python -m pybind11 --includes` IBN_u_opt.cpp -o extremal.so

//compilation on mac: g++ -O3 -Wall -shared -std=c++11 -undefined dynamic_lookup `python -m pybind11 --includes` IBN_u_opt.cpp -o extremal.so
