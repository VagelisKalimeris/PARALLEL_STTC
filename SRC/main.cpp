/******************************************************************************
*******************************************************************************
*                                                                             *
* PROJECT NAME: STTC Analyses                                                 *
*                                                                             *
* FILE NAME: main.cpp                                                         *
*                                                                             *
*******************************************************************************
******************************************************************************/


#include <iostream>
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cstdio>
#include <math.h>

#include "../INCLUDE/common.hpp"
#include "../INCLUDE/p_p_null_dist.hpp"
#include "../INCLUDE/cond_null_dist.hpp"
#include "../INCLUDE/pairs_STTC.hpp"
#include "../INCLUDE/triplets_STTC.hpp"
#include "../INCLUDE/print.hpp"

using namespace std;

/******************************************************************************
* FUNCTION NAME: main                                                         *
*                                                                             *
* ARGUMENTS: The total numbers Circular Shifts and the size of the tile Dt.   *
*                                                                             *
* PURPOSE: This main function is for testing purposes only.                   *
*                                                                             *
* RETURNS: 0 on program completion.                                           *
*                                                                             *
* I/O: Opens/Reads a file containing the spike trains.                        *
*                                                                             *
******************************************************************************/
int main(int argc, char const *argv[])
{
// Command line arguments check
    if (argc != 4) {
        cout<<"Error! Wrong parameter count!"<<endl;
        return 0;
    }
// Command Line Arguments. First give random sample size, then tile size. 
    const int circ_shifts_num = atoi(argv[1]), Dt = atoi(argv[2]);
    
// Open File
    ifstream data, astros;
    data.open((string("DATASETS/") + argv[3]).c_str(), ifstream::in);
    if (!data.is_open()) {
        cout<<"Error opening input dataset file!"<<endl;
        return 0;
    }
    astros.open((string("ASTROCYTES/") + argv[3]).c_str(), ifstream::in);
    if (!astros.is_open()) {
        cout<<"\nProblem opening input astrocytes file!"<<endl
            <<"\nContinuing without astrocyte removal!"<<endl;
    }
    string line;
    
// Get astrocytes
    vector<int> astrocytes;
    while (getline(astros, line)) {
        astrocytes.push_back(atoi(line.c_str()) - 1);
    }
    const int astro_size = astrocytes.size();
    
// Get total number of neurons from file
    getline(data, line);
    const int ttl_neurons = line.length();
    const int neurons = line.length() - astro_size;
    data.seekg(0, data.beg);
    
// Our main data structure and astrocyte list
    vector<int> spike_trains[ttl_neurons];
    
// Store each neuron's firing (1's) to the data structure
    int total_time_samples = 0;
    while (getline(data, line)) {
        int astros_count = 0;
        int astrocyte = 0;
        if (astro_size) {
            astrocyte = astrocytes[0];
        }
        int push_count = 0;
        for (int neur = 0; neur < ttl_neurons; ++neur) {
            int pos;
            if (astro_size && neur == astrocyte) {
                pos = neurons + astros_count++;
                astrocyte = astrocytes[astros_count];
            }
            else {
                pos = push_count++;
            }
            if (line[neur] == '1') {
                spike_trains[pos].push_back(total_time_samples);
            }
        }
        total_time_samples++;
    }
    
// Make the mapping from virtual to real neuron's number
    int map[neurons];
    int astro = 0;
    int astrocyte = 0;
    if (astro_size) {
        astrocyte = astrocytes[0];
    }
    for (int neur = 0; neur < neurons; ++neur) {
        while (astro_size && (neur + astro) == astrocyte) {
            astrocyte = astrocytes[(++astro) % astro_size];
        }
        map[neur] = neur + astro + 1;
    }

// Close input files
    data.close();
    astros.close();
    
// Print message that computation is starting
    cout<<"\nComputing dataset "<<argv[3]<<" with Dt = "<<argv[2]
        <<" and control group = "<<argv[1]<<"."<<endl;

// Print the data structure and total number of firings in experiment
    char str[33];
    sprintf(str, "%d", circ_shifts_num);
    const string shifts_s = string(str);
    sprintf(str, "%d", Dt);
    const string Dt_s = string(str);
    ofstream info;
    info.open(("RESULTS/" + string(argv[3]) + "_" + shifts_s + "-shifts_" + 
                                    Dt_s + "-dt_neurons_info.txt").c_str());
    if (!info.is_open()) {
        cout<<"Error opening results neurons info file!"<<endl;
        return 0;
    }
    print_all_spikes(spike_trains, ttl_neurons, astrocytes, info, 
                                            string(argv[3]), shifts_s, Dt_s);
    
// Start random sequence
    srand(time(NULL));
    
// Time lines' arrays
    int tl_sizes[neurons];
    int* tl_array[neurons];
    int tl_size_max = 0;
    vector<int> spike_train;
    for (int neur = 0; neur < neurons; ++neur) {
        spike_train = spike_trains[neur];
        int tl_size = spike_train.size();
        tl_sizes[neur] = tl_size;
        tl_array[neur] = (int *)malloc(tl_size * sizeof(int));
        for (int ts = 0; ts < tl_size; ++ts) {
            tl_array[neur][ts] = spike_train[ts];
        }
        if (tl_size_max < tl_size) {
            tl_size_max = tl_size;
        }
    }
    
// All T for pairs
    double T_Aplus[neurons];
    double T_Bminus[neurons];
    #pragma omp parallel for
    for(int neur = 0; neur < neurons; ++neur) {
        int* tl = tl_array[neur];
        int tl_size = tl_sizes[neur];
        T_Aplus[neur] = T_A_plus(tl, tl_size, total_time_samples, Dt);
        T_Bminus[neur] = T_B_minus(tl, tl_size, total_time_samples, Dt);
    }
    
// Significant limit
//    int sgnfcnt_limit[neurons][neurons];
    
// Reduced spiketrain T for triplets
//    double T_Aplus_tripl[neurons][neurons];
    
// Calculate per pair STTC
    ofstream pairs;
    pairs.open(("RESULTS/" + string(argv[3]) + "_" + shifts_s + "-shifts_" + 
                                            Dt_s + "-dt_pairs.csv").c_str());
    if (!pairs.is_open()) {
        cout<<"Error opening results pairs file!"<<endl;
        return 0;
    }
    pairs<<"NeuronA,NeuronB,STTC,CtrlGrpMean,CtrlGrpStDev,CtrlGrpMedian,Percentile\n";
    ofstream pairs_cg;
    pairs_cg.open(("RESULTS/" + string(argv[3]) + "_" + shifts_s + "-shifts_" + 
                                            Dt_s + "-dt_pairs_cg.csv").c_str());
    if (!pairs_cg.is_open()) {
        cout<<"Error opening results pairs_cg file!"<<endl;
        return 0;
    }
    pairs_cg<<"NeuronA,NeuronB,STTC,CtrlGrpMean,CtrlGrpStDev,CtrlGrpMedian,Percentile";
    for (int i = 0; i < circ_shifts_num; ++i) {
        pairs_cg<<",STTC_"<<i+1;
    }
    pairs_cg<<"\n";
    for (int a = 0; a < neurons; a++) { // Neuron A
        int* tl_A = tl_array[a];
        int tl_A_size = tl_sizes[a];
        if (tl_A_size == 0) {continue;}
        int a_real = map[a];
        #pragma omp parallel
        {
            double tAp = T_Aplus[a];
        // Shifted spike trains will be copied here
            int* to_shift = (int *)malloc(tl_size_max * sizeof(int));
        // STTC values of shifted spike trains
            double shifted_res_arr[circ_shifts_num];
            #pragma omp for
            for (int b = 0; b < neurons; b++) { // Neuron B
                if (a == b) {continue;} // Skip same neurons
                int* tl_B = tl_array[b];
                int tl_B_size = tl_sizes[b];
                if (tl_B_size == 0) {continue;}
/*                sgnfcnt_limit[a][b] = sign_trpl_limit(tl_A, tl_A_size, tl_B, 
                                                                tl_B_size, Dt);
                if (sgnfcnt_limit[a][b] > 5) {
                    T_Aplus_tripl[a][b] = T_A_plus_tripl(tl_A, tl_A_size, 
                                    tl_B, tl_B_size, total_time_samples, Dt);
                }*/
                double tBm = T_Bminus[b];
                double pair_sttc = STTC_A_B(tl_A, tl_A_size, tl_B, tl_B_size, 
                                                                Dt, tBm, tAp);
                if (pair_sttc == 2.0) {continue;}
                int denominator = circ_shifts_num;
                double mean = 0;
                for (int shift = 0; shift < circ_shifts_num; shift++) {
                    unsigned int random = random_gen(total_time_samples);
                    circular_shift(to_shift, tl_A, tl_A_size, random, 
                                                        total_time_samples);
                    double tAp_s = T_A_plus(to_shift, tl_A_size, total_time_samples, 
                                                                        Dt);
                    shifted_res_arr[shift] = STTC_A_B(to_shift, tl_A_size, 
                                                tl_B, tl_B_size, Dt, tBm, tAp_s);
                    if (shifted_res_arr[shift] == 2.0) {
                        --denominator;
                    }
                    else {
                        mean += shifted_res_arr[shift];
                    }
                }
            // Pair is valid if control group has at least 80% valid STTC values
                if (double(denominator) < (0.8 * circ_shifts_num)) {continue;}
            // Mean is equal to sum of valid divided by size of valid null STTC values
                mean /= denominator;
            // Sorting of null STTC values
                sort(shifted_res_arr, (shifted_res_arr + circ_shifts_num));
                int b_real = map[b];
                
                double st_dev = 0.0;
                string null_STTC;
                char buffer[32];
                for (int i = 0; i < denominator; i++) {
                    double value = shifted_res_arr[i];
                    st_dev += pow(value - mean, 2.0);
                    
                    sprintf(buffer, "%f", value);
                    null_STTC += ',' + string(buffer);
                }
                sprintf(buffer, "%f", 2.0);
                for (int i = denominator; i < circ_shifts_num; i++) {
                    null_STTC += ',' + string(buffer);
                }
                st_dev = sqrt(st_dev / denominator);
                
                double median;
                if (denominator % 2) {
                    median = shifted_res_arr[denominator / 2];
                }
                else {
                    median = (shifted_res_arr[denominator / 2 - 1] + 
                                shifted_res_arr[denominator / 2]) / 2.0;
                }
                
                int pos = 0; 
                while (pos < denominator && 
                                    shifted_res_arr[pos] <= pair_sttc) {
                    ++pos;
                }
                double percentile = pos / double(denominator);
                
                #pragma omp critical
                pairs << a_real << ',' << b_real << ',' << pair_sttc 
                                << ',' << mean << ',' << st_dev 
                                << ',' << median << ',' << percentile << '\n';
                #pragma omp critical
                pairs_cg << a_real << ',' << b_real << ',' << pair_sttc 
                                << ',' << mean << ',' << st_dev 
                                << ',' << median << ',' << percentile 
                                << null_STTC << '\n';
            }
            free(to_shift);
        }
    }
    pairs.close();
    pairs_cg.close();
    
    /*
// Calculate conditional STTC
    ofstream triplets;
    triplets.open(("RESULTS/" + string(argv[3]) + "_" + shifts_s + 
                            "-shifts_" + Dt_s + "-dt_triplets.csv").c_str());
    if (!triplets.is_open()) {
        cout<<"Error opening results triplets file!"<<endl;
        return 0;
    }
    triplets<<"NeuronA,NeuronB,NeuronC,STTC,CtrlGrpMean,CtrlGrpMedian,Percentile\n";
    for (int a = 0; a < neurons; a++) { // Neuron A
        int* tl_A = tl_array[a];
        int tl_A_size = tl_sizes[a];
        if (tl_A_size == 0) {continue;}
        int a_real = map[a];
        #pragma omp parallel
        {
        // Shifted spike trains will be copied here
            int* to_shift = (int *)malloc(tl_size_max * sizeof(int));
        // STTC values of shifted spike trains
            double shifted_res_arr[circ_shifts_num];
            #pragma omp for
            for (int c = 0; c < neurons; c++) { // Neuron C
                if (a == c) {continue;} // Skip same neurons
                int* tl_C = tl_array[c];
                int tl_C_size = tl_sizes[c];
                if (tl_C_size == 0) {continue;}
                int tl_redA_size = sgnfcnt_limit[a][c];
                if (tl_redA_size <= 5) {
                    continue; // Reduced A spike train has < 5 spikes
                }
                double tApt = T_Aplus_tripl[a][c];
                int c_real = map[c];
                for (int b = 0; b < neurons; b++) { // Neuron B
                    if (b == a || b == c) {continue;} // Skip same neurons
                    int* tl_B = tl_array[b];
                    int tl_B_size = tl_sizes[b];
                    if (tl_B_size == 0) {continue;}
                    double tBm = T_Bminus[b];
                    double trip_sttc = STTC_AB_C(tl_A, tl_A_size, 
                                        tl_B, tl_B_size, tl_C, tl_C_size, 
                                        Dt, tBm, tApt, tl_redA_size);
                    if (trip_sttc == 2.0) {continue;}
                    int denominator = circ_shifts_num;
                    double mean = 0;
                    for (int shift = 0; shift < circ_shifts_num; shift++) {
                        unsigned int random = random_gen(total_time_samples);
                        circular_shift(to_shift, tl_C, tl_C_size, random, 
                                                        total_time_samples);
                        int tl_redA_size_s = sign_trpl_limit(tl_A, tl_A_size, 
                                                    to_shift, tl_C_size, Dt);
                        if (tl_redA_size_s > 5) {
                            double tApt_s = T_A_plus_tripl(tl_A, tl_A_size, to_shift, 
                                            tl_C_size, total_time_samples, Dt);
                            shifted_res_arr[shift] = STTC_AB_C(tl_A, tl_A_size, 
                                        tl_B, tl_B_size, to_shift, tl_C_size, 
                                        Dt, tBm, tApt_s, tl_redA_size_s);
                            if (shifted_res_arr[shift] == 2.0) {
                                --denominator;
                            }
                            else {
                                mean += shifted_res_arr[shift];
                            }
                        }
                        else {
                            --denominator;
                        }
                    }
                    if (double(denominator) < (0.8 * circ_shifts_num)) {continue;}
                    mean /= denominator;
                    sort(shifted_res_arr, (shifted_res_arr + circ_shifts_num));
                    double median;
                    if (denominator % 2) {
                        median = shifted_res_arr[denominator / 2];
                    }
                    else {
                        median = (shifted_res_arr[denominator / 2 - 1] + 
                                    shifted_res_arr[denominator / 2]) / 2.0;
                    }
                    int pos = 0; 
                    while (pos < denominator && 
                                    shifted_res_arr[pos] <= trip_sttc) {
                        ++pos;
                    }
                    int b_real = map[b];
                    double percentile = pos / double(denominator);
                    #pragma omp critical
                    triplets << a_real << ',' << b_real << ',' << c_real 
                                << ',' << trip_sttc << ',' << mean 
                                << ',' << median << ',' << percentile << '\n';
                }
            }
            free(to_shift);
        }
    }
    triplets.close();
    */
    
// Free memory
    for (int neur = 0; neur < neurons; ++neur) {
        free(tl_array[neur]);
    }
    
// Close output files
    info.close();
    cout<<"\nComputation has ended successfully. Output files are located "
        <<"inside RESULTS directory."<<endl;
    return 0;
}