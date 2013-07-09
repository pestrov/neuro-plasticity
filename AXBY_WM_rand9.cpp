/*************************************************************************************/
/*Step 0: History of Program                                                         */
/*************************************************************************************/

// 4 stimulus recurrent random sparse spiking network. A, X, B, Y. 
// Updated with network Pfister & Gerstner STDP rule (network only)
// 16 input pool structured.  Mark B. 10/14/08
// inputs rewritten (started) 11/23/08 Mark B. to be malleable
// They are no longer hard-coded to # pools with fixed input cells
// pool size is plastic along with inputs being probalisitic to any input 
// 9/10/09:  loads of errors fixed Jul-Aug 2010. - Voltage LTPi "if" statements corrected for order of precedence,
// bool to only potentiate existing synapses.
// 8/11/10: limit dwSum after multiplying by dw.
// 10/20/10: 4 trials of no plasticity every 80 trials to check specificity.
// 1/1/11 - cueFR adjusted for 3 epochs, normalized by length (cue1, cue1+cue2, delay, 2nd cue + overlaying delay).  DL_Out only during 2nd stimulus, uses correct 2nd cueFR
// 3/15/11: - Fixed Input triplet STDP normalization to pW02, not pW0, (kept for reading matrices), Corrected itime (was integer) in LTPi window code.  
// 7/22/11 - split tuning of AL Ecells/Icells, initialized eBinRate to zero for bin 0, heterogeneity is now (rand1() - 0.5) with a printout for Matlab to verify output, added binned meanFR for all AL cells viewable in Matlab with script, other changes/additions...
// 6/18/13 - First complete draft of updated version done by William.  Created functions and header files, and improved readability.

/*************************************************************************************/
/* Step 1: Include Statements                                                        */
/*************************************************************************************/

// Right now, we're using namespace std.  I might take this out in the future because it's probably better
// programming practice to declare our std functions explicitly each time.
using namespace std;
// Here, we include our basic preprocessing directives.  
#include <iostream> 
#include <fstream> 
#include <iomanip>
#include <string>
#include <time.h>
#include <stdio.h>
#include <vector>
#include <cstdlib>
#include <algorithm>
#include <cmath>
#include <sstream>
//#include <mathimf.h> // Intel math library for the cluster only

// Now, we include our custom header files.  These contain the constants used
// throughout our program.  Organization can definitely be improved.
// When I was writing a Java version of this program, I used getter/setter 
// methods for these constants, using a Scanner object to draw the values
// from .DAT files.  That being said, after research and some thought, I don't
// really see why these methods are particularly useful or necessary in this
// situation; getter/setter methods are often frowned upon in the programming
// community, anyway.  Thus, the header files simply contain const declarations
// of the values.  In order to alter any given value, one should go into the
// given header file and change the value there. 

#include "BoolInputs.h" // Boolean inputs 
#include "DL_net.h"  // Decision layer, and AL-to-DL weight parameters
#include "ExcitatoryAndInhibitory.h" // Excitatory and inhibitory constants
#include "Time.h" // Time-related constants
#include "SynapticReversal.h" // Synaptic reversal constants
#include "MBdir_AL_DL.h"  // Network size, gain, homeostasis
#include "ExponentialDecays.h" // Exponential decay constants
#include "FacilitationAndDepression.h" // Facilitation and depression constants
#include "HomeostasisAndLTPi.h" // Homeostasis and LTPi constants
#include "InputCellsAndProb.h" // Input cells and input probability arrays 
#include "InputSpecific.h" // Input specific constants
#include "LeakyIntegrateAndFire.h" // Leaky integrate and fire constants
#include "HeterogeneityVectors.h" // Heterogeneity vectors
#include "RewardPredictionHistory.h" // Reward Prediction Error code - history component
#include "StructInputs.h" // Structural inputs
#include "SynapticCond.h" // Synaptic conductance constants
#include "PGSTDPAmplitude.h" // PGSTD Amplitude constants
#include "Weights.h" // Weight constants
#include "BinsAndStim.h" // Bins and stim constants
#include "TrialValues.h" // Trial-specific constants (constants that have to be reset for each new trial)

// MersenneTwister.h is used to create our random number generators
// (rand1, rand2, and rand3) that are used throughout the program.
// MersenneTwister will provide more randomized results than the default 
// random number generator provided by C++.
#include "MersenneTwister.h"
MTRand rand1;  // This instance is used when working with connections.
MTRand rand2;  // This instance is used when working with noise.
MTRand rand3; //  This instance is used when working with Poisson inputs.
// The seeds below will be used in our rand objects.  Note that these are 
// initially assigned default values; they will be altered below according 
// to the user's command line input.
int seed1 = 1;
int seed2 = seed1; 
int seed3 = seed1; 
// Currently, seed3 remains 1 throughout the program.

/*************************************************************************************/
/* Step 2: Declaration of Functions                                                  */
/*************************************************************************************/

// Here, we declare the functions that will be used throughout the program.

//// SEED, NUMIN, AND DIGITS FUNCTIONS ////
// Here, we find the seeds (at least the first two!) for our MersenneTwister objects.
void findSeeds(int& seed1, int& seed2, int& numIn) {
    if (numIn >= 10000) {
        seed2 = numIn / 10000;
        numIn = numIn - 10000 * seed2;
    }
    if (numIn >= 100) {
        seed1 = numIn / 100;
        numIn = numIn - 100 * seed1;
    }
    // Let's print numIn so we can make sure that we received the command line
    // parameter correctly.
    cout << "numIn " << numIn << endl;
}

// Here, we derive the digits we'll use to choose the number of input cells and the input
// probability from our numIn value.
void findDigits(int& digit1, int& digit2, int& numIn) {
    if (numIn > 0) {
        digit1 = (numIn - 1) / 5;
        digit2 = numIn - 5 * digit1;
    }
    digit1 += 1;
    // We print the two digits.
    cout << "digit1: " << digit1 << endl;
    cout << "digit2: " << digit2 << endl;
}

// Here, we use the digits we just found to find the number of input cells and the input
// probability.
void findinputs(int& digit1, int& digit2, int& inputCells, double& inputProb, int* inputCellArray, double* inputProbArray) {
    if (digit1 > 0) {
        inputCells = inputCellArray[digit1 - 1];
    }
    if (digit2 > 0) {
        inputProb = inputProbArray[digit2 - 1];
    }
    // We print the inputs for our cells and probability values.
    cout << "inputCells: " << inputCells << endl;
    cout << "inputProb: " << inputProb << endl;
}

//// REWARD HISTORY FUNCTIONS ////

// Here, we fill out the preliminary reward history.
void fillHistFactor(vector<double>& histFactor, double& sumHistFactor, double histDecay, int nHist) {
    histFactor[0] = (1 - histDecay) / (1 - exp(log(histDecay) * nHist)); // Exponential decay trace
    sumHistFactor = histFactor[0];
    for (int i = 1; i < nHist; i++) {
        histFactor[i] = histFactor[i-1] * histDecay; 
        sumHistFactor += histFactor[i];
    }
    // Our sumHistFactor should be 1.  We'll check that this is true in our print statement below.
    cout << "sumHistFactor = 1? " << sumHistFactor << endl;
}

//// TIME FUNCTIONS ////

// Here, we fill out our time vector.  We also initialize cue length, delay length, and stimulus delivery times.
void fillTime(vector<double>& t, vector<double>& tOn, vector<double>& tOff, double& cue1Length, double& cue2Length, double& delay1Length, double& delay2Length, int tMax, int nT, double dt) {
    for (int i = 0; i <= nT; i++) {
        t[i] = dt * double(i);
    }
    // Stimulus delivery times, ex. A&B if no persistence, A-delay-B if persistence
    tOn[0] = 0.250;
    tOff[0] = 0.75;
    tOn[1] = 1.25;
    tOff[1] = 1.75;
    // Cue and delay lengths
    cue1Length = tOff[0] - tOn[0];
    cue2Length = tOff[1] - tOn[1];
    delay1Length = tOn[1] - tOff[0];
    delay2Length = tMax - tOff[1];
}

//// GOAL FUNCTIONS ////

// Here, we fill our goal rate vectors, with noise for heterogeneity.
void fillGoalVectors(vector<double>& r_Goal, vector<double>& r_Goal_I_to_E, vector<double>& r_DGoal, int NE, int NI, int dNE, MTRand& rand1, int rgE, double rgIE0, int rgI, int rgD) {
    for (int j = 0; j < NE; j++) {
        r_Goal[j] = rgE + 5 * (rand1() - 0.5);
        r_Goal_I_to_E[j] = rgIE0 + r_Goal[j];
    }
    for (int j = NE; j < NE + NI; j++) {
        r_Goal[j] = rgI + 2.5 * (rand1() - 0.5);
    }
    // Goal rate vector for Associative input->Decision layer
    for (int j = 0; j < dNE; j++) {
        r_DGoal[j] = rgD;
    }
}

//// FILLING OUT PRELIMINARY WEIGHT MATRICES ////

// Here, we fill our NE-NE, NI-NE, NI-NI, and NE-NI loops, thus creating our weight matrix W.
void fillLoops(int NE, int neurons, double connProb, double wEE_xFactor, double W0, double ranWStrengths, double wIE_ConnProb, double wEI_ConnProb, double wIE0, double wII0, double wEI0, int dNeurons, MTRand& rand1, vector< vector<double> >& W, vector< vector<double> >& dWstruct) {
    //NE-NE loop
    for (int i = 0; i < NE; i++) {
        for (int j = 0; j < NE; j++) {
            if (rand1() < connProb) { // WEE_conn_prob
                W[i][j] = wEE_xFactor * W0 * (1 + ranWStrengths * (rand1() - 0.5));
            }
        }
    W[i][i] = 0; // No self-synapses
    }
    // NI->NE
    for (int i = NE; i < neurons; i++) {
        for (int j = 0; j < NE; j++) {
            if (rand1() < wIE_ConnProb) {
                W[i][j] = wIE0 * (1 + ranWStrengths * (rand1() - 0.5));
            }
        }
    }
    // NI-NI
    for (int i = NE; i < neurons; i++) {
        for (int j = NE; j < neurons; j++) {
                W[i][j] = wII0* (1 + ranWStrengths * (rand1() - 0.5));
            }
    }
    //NE->NI (Recurrent I)
    for (int i = 0; i < NE; i++) {
        for (int j = NE; j < neurons; j++) {
            if (rand1() < wEI_ConnProb) {
                W[i][j] = wEI0 * (1 + ranWStrengths * (rand1() - 0.5));
            }
        }
    }
}

// Here, we randomize eW.
void randomize_eW(int nIn, int NE, int NI, MTRand& rand1, double inputProb, vector< vector<double> >& eW) {
    //A
    for (int i = 0; i < int(nIn / 4); i++) {
        for (int j = 0; j < NE + NI; j++) {
            if (rand1() < inputProb) {
                eW[i][j] = 1; // 0-1, 2-3, 4-5, 6-7
            }
        }
    }
    //X   
    for (int i = int(nIn / 4); i < int(nIn / 2); i++) {
        for (int j = 0; j < NE + NI; j++) {
            if (rand1() < inputProb) {
                eW[i][j] = 1; // 0-1, 2-3, 4-5, 6-7
            }
        }
    }
    //B
    for (int i = int(nIn / 2); i < int(3 * nIn / 4); i++) {
        for (int j = 0; j < NE + NI; j++) {
            if (rand1() < inputProb) {
                eW[i][j] = 1; // 0-1, 2-3, 4-5, 6-7
            }
        }
    }       
    //Y
    for (int i = int(3 * nIn / 4); i < nIn; i++) {
        for (int j = 0; j < NE + NI; j++) {
            if (rand1() < inputProb) {
                eW[i][j] = 1; // 0-1, 2-3, 4-5, 6-7
            }
        }
    }
}

// Here, we account for the heterogeneity of pW.
void heterogeneity_pW(int nIn, int NE, double pW0, double pW0I, double ranWStrengths, MTRand& rand1, int neurons, int dNeurons, vector< vector<double> >& pW, vector< vector<double> >& eW, bool feedForwardI) {
    for (int i = 0; i < nIn; i++) {
        for (int j = 0; j < NE; j++) {
            pW[i][j] =(eW[i][j]) * pW0 * (1 + ranWStrengths * (rand1() - 0.5));
        }
        if (feedForwardI) {
            for (int j = NE; j < NE + NI; j++) {
                pW[i][j] =(eW[i][j]) * pW0I * (1 + ranWStrengths * (rand1() - 0.5));
            }
        }
    }
}

// Here, we initialize external Poisson inputs for pW (if we are using Poisson inputs).
void poisson_pW(int nIn, int neurons, int dNE, int dNEPool, double pool1Weight, double pool2Weight, vector< vector<double> >& pW) {
    // Testing so bias external Poisson inputs
    for (int i = 0; i < nIn; i++) {
        for (int j = neurons; j < neurons + dNE; j++) {
            if (j < neurons + dNEPool) {
                pW[i][j] = pool1Weight;
            } else {
                pW[i][j] = pool2Weight;
            }
        }
    }
}

// Here, we account for our pools in our weight matrix W.
void pools_W(int NE, int neurons, int dNeurons, int dNE, double dLayerInput_W0, int dNPools, int dNEPool, double WEE, double dW0, int dNIPool, double WEI, double WIE, double WII, vector< vector<double> >& W) {
   for (int i = 0; i < NE; i++) {
        for (int j = neurons; j < neurons + dNE; j++) {
             W[i][j] = dLayerInput_W0;
        }
    }  
    for (int pool = 1; pool <= dNPools; pool++) {
        for (int i = neurons + dNEPool * (pool-1); i < neurons + dNEPool * pool; i++) {
            for (int j = neurons + dNEPool * (pool-1); j < neurons + dNEPool * pool; j++) {
                W[i][j] = wEE * dW0; // 1.9 2.45 2.7 (FRss ~= 10Hz) or Poisson 2.1-2.2 // 2.75
            }
        }
    }
    // NEpools  -> Inhibitory pool - WEI
    for (int i = neurons; i < neurons + dNEPool; i++) { //column
        for (int j = neurons + dNE + dNIPool; j < neurons + dNeurons; j++) { // row
            W[i][j] = wEI * dW0; //2.0 2.40 // 5.0
        }
    }
    // NEpools  -> Inhibitory pool - WEI
    for (int i = neurons + dNEPool; i < neurons + dNE; i++) { //column
        for (int j = neurons + dNE; j < neurons + dNE + dNIPool; j++) { // row
            W[i][j] = wEI * dW0; //2.0 2.40 // 5.0
        }
    }
    // Inhibitory pool -> NEpools - WIE
    for (int i = neurons + dNE; i < neurons + dNE + dNIPool; i++) {
        for (int j = neurons; j < neurons + dNEPool; j++) {
            W[i][j] = wIE * dW0;  //1.0 1.45 2.5
        }
    } 
    // Inhibitory pool -> NEpools - WIE
    for (int i = neurons + dNE + dNIPool; i < neurons + dNeurons; i++) {
        for (int j = neurons + dNEPool; j < neurons + dNE; j++) {
            W[i][j] = wIE * dW0;  //1.0 1.45 2.5
        }
    }
    // Inhibitory pool -> Inhibitory pool - WII
    for (int i = neurons + dNE; i < neurons + dNeurons; i++) {
        for (int j = neurons + dNE; j < neurons + dNeurons; j++) {
            W[i][j] = wII * dW0; //1.0
        }
    }        
}

// Here, we modify pW with our excitatory and inhibitory weight mods.
void modify_pW(int nIn, int NE, int NI, double pW02, double pWOI2, vector< vector<double> >& pW) {
    for (int i = 0; i < nIn; i++) {
        for (int j = 0; j < NE; j++) {
            pW[i][j] *= pW02; // Excitatory mod
        }
        for (int j = NE; j < NE + NI; j++) {
            pW[i][j] *= pW0I2; // Inhibitory mod
        }
    }
}

// Here, we normalize W and set up our C1 and C3 bool matrices.
void normalize_And_C1C3Bools(int neurons, int dNeurons, int nIn, vector< vector<double> >& W1, vector< vector<double> >& W, vector< vector<double> >& pW, bool C1[][totalNeurons], bool C3[][totalNeurons]) {
    // Network connection bool and normalize matrices
    for (int i = 0; i < neurons + dNeurons; i++) {
        for (int j = 0; j < neurons + dNeurons; j++) {
            W1[i][j] = W[i][j];   // for STDP normalization
            if (W[i][j] > 0) {
                C1[i][j] = 1; 
            }
            // Separate bool now for PG, STDP rules to only apply when input cells projects->neuron
            if (i < nIn) {
                if (pW[i][j] > 0) {
                    C3[i][j] = 1; // C3 is not being currently for anything!
                }
            }
        }
    }
}

// Here, we initialize our heterogeneity vectors.
void initialize_Heterogeneity_Vectors(int neurons, int NE, int dNE, int dNeurons, MTRand& rand1, vector<double>& E, double E_s, vector<double>& tau_M, double tau_M_S, vector<double>& g_L, double AL_g_L_s, vector<double>& vth_E, double vth_E_sA, vector<double>& vReset_E, double vReset_E_s, vector<double>& tRef_E, double tRef_E_s, double vth_I_sA, double vReset_I_s, double tRef_I_s, double g_L_s, double vth_E_sD, double vth_I_sD) {
    // Static value * (heterogeneous offset value, which is mean static value subtracting heterogeneous offset!) * (random number from the uniform set=> [-0.5...0.5] )
    for (int i = 0; i < neurons; i++) {
        E[i] = E_s + (2.5e-3) * (rand1() - 0.5);  //-70e-3
        if (i < NE) {
            tau_M[i] = tau_M_S + (0.75e-3) * (rand1() - 0.5);  // 10e-3
            g_L[i] = AL_g_L_s + (1.0e-6) * (rand1() - 0.5);  // 23e-6
            vth_E[i] = vth_E_sA + (2.0e-3) * (rand1() - 0.5);  // 50e-3
            vReset_E[i] = vReset_E_s + (2.0e-3) * (rand1() - 0.5);  // 60e-3
            tRef_E[i] = tRef_E_s + (0.25e-3) * (rand1() - 0.5);   // .002
        } else {
            tau_M[i] = tau_M_S + (0.75e-3) * (rand1() - 0.5);  // 10e-3
            g_L[i] = AL_g_L_s + (1.0e-6) * (rand1() - 0.5);  // 23e-6
            vth_I[i] = vth_I_sA + (2.0e-3) * (rand1() - 0.5);  // 50e-3
            vReset_I[i] = vReset_I_s + (2.0e-3) * (rand1() - 0.5);  // 60e-3
            tRef_I[i] = tRef_I_s + (0.25e-3) * (rand1() - 0.5);  // .001
        }
    }
    // Decision cell parameters (in red) are from Wang 02 Neuron
    // Decision pools: can't have heterogenity because of bias
    // Xiao-Jing's parameters
    for (int i = neurons; i < neurons + dNE; i++) {  // Ecells
        E[i] = E_s; // -70e-3;  // -70mV
        tau_M[i] = 20e-3; //tau_M_S; // 20e-3; // 20mS 
        g_L[i] = g_L_s; //25e-6; // 25nS
        vth_E[i] = vth_E_sD; //vth_E_s; //-48e-3;  // -50mV
        vReset_E[i] = -55e-3; //vReset_E_s; // -60e-3; // -55mV
        tRef_E[i] = tRef_E_s; //.002; // 2ms
    }
    for (int i = neurons + dNE; i < neurons + dNeurons; i++) { // Icells
        E[i] = E_s; //-70e-3;  // -70mV
        tau_M[i] = tau_M_S; //10e-3 ; // 10ms
        g_L[i] = 30e-6; //g_L_s; //20e-6;  // 20nS
        vth_I[i] = vth_I_sD; //-50e-3; // -50mV
        vReset_I[i] = -55e-3; //vReset_I_s; // -55e-3; // -55mV
        tRef_I[i] = tRef_I_s; // .001; // 1ms
    }        
}

// Here, we initialize stim shuffle, the vector used to choose our stim value.
void initialize_StimShuffle(int inputs, vector< int>& stimShuffle) {
    for (int i = 0; i < inputs; i++) {
        stimShuffle[i] = i + 1;
    }
}

// Here, we print out the values of W.
void print_W(int neurons, int dNeurons, vector< vector<double> >& W) {
    // Prints initial W
    ofstream W02out;
    W02out.open("W02.dat");
    for (int i = 0; i < neurons + dNeurons; i++) {
        for (int j = 0; j < neurons + dNeurons; j++) {
            W02out << W[i][j] << " ";
        }
        W02out << endl;
    }
    W02out.close();
}

// Here, we print out the values of pW.
void print_pW(int nIn, int NE, int NI, int dNeurons, vector< vector<double> >& pW) {
    // Write out Input weight matrix
    ofstream pW02_out;  // output initial input weight matrix
    pW02_out.open("pW02.dat");
    for (int i = 0; i < nIn; i++) {
        for (int j = 0; j < NE + NI + dNeurons; j++) {
            pW02_out << pW[i][j] << " ";
        }
        pW02_out << endl;
    }
    pW02_out.close();   
}

// Here, we initialize various starting values to 0 for each trial.
// This adds some time to each trial, but it also makes our main method much neater.
// If we want to make the program as efficient time-wise as possible, we can just 
// declare those values at the beginning of each trial.
void trial_Values_to_0(vector <vector<double> >& spikeTimes, vector< int>& numberSpikes, vector<double>& gSyn_I_Vec, vector<double>& gSyn_E_Vec, vector<double>& gpSyn_Vec, double& g_Tot, vector<double>& vInf, vector< vector<double> >& V, vector<double>& lastSpike, vector<double>& syn, vector<double>& syn_I, vector<double>& syn_E, double& tau_Eff, double& pSynMax, double& synMax_I, double& synMax_E, double& synMax_NMDA, int& someIndex, vector<double>& gSyn_NMDA, vector<double>& g_NMDA, vector<double>& Vth, vector<double>& syn_NMDA, vector<double>& g_Ref, vector<double>& meanFiringRate, vector<double>& W_Shift, vector<double>& W_Shift_I, vector<double>& DW_Shift, vector<double>& Fac, vector<double>& synE_Fac, vector<double>& synNMDA_Fac, vector<double>& Dep, vector<double>& synE_Dep, vector<double>& synNMDA_Dep, int neurons, int dNeurons, int max_CellSpikes, int nT) {
    // Zero out integers and doubles
    someIndex = 0;
    g_Tot = 0;
    tau_Eff = 0;
    pSynMax = 0;
    synMax_I = 0;
    synMax_E = 0;
    synMax_NMDA = 0;
    // Zero out 1D vectors
    for (int i = 0; i < neurons + dNeurons; i++) {
        numberSpikes[i] = 0;
        gSyn_I_Vec[i] = 0;
        gSyn_E_Vec[i] = 0;
        gpSyn_Vec[i] = 0;
        vInf[i] = 0;
        lastSpike[i] = 0;
        syn[i] = 0;
        syn_I[i] = 0;
        syn_E[i] = 0;
        gSyn_NMDA[i] = 0;
        g_NMDA[i] = 0;
        Vth[i] = 0;
        syn_NMDA[i] = 0;
        g_Ref[i] = 0;
        meanFiringRate[i] = 0;
        W_Shift[i] = 0;
        Fac[i] = 0;
        synE_Fac[i] = 0;
        synNMDA_Fac[i] = 0;
        Dep[i] = 0;
        synE_Dep[i] = 0;
        synNMDA_Dep[i] = 0;
    }
    for (int i = 0; i < NE; i++) {
        W_Shift_I[i] = 0;
    }
    for (int i = 0; i < dNeurons; i++) {
        DW_Shift[i] = 0;
    }
    // Zero out 2D vectors
    for (int i = 0; i < neurons + dNeurons; i++) {
        for (int j = 0; j < max_CellSpikes; j++) {
            spikeTimes[i][j] = 0;
         }
    }
    for (int i = 0; i < neurons + dNeurons; i++) {
        for (int j = 0; j < nT+1; j++) {
            V[i][j] = 0;
        }
    }
}

// Here, we initialize our lastSpike matrix.
void initialize_LastSpike(int neurons, int dNeurons, int NE, int dNE, vector<double>& tRef_E, vector<double>& tRef_I, vector<double>& lastSpike) {
    // Initialize lastSpike 
    for (int l = 0; l < neurons + dNeurons; l++) {
        lastSpike[l] = 2 * tRef_E[l];
        if ((l >= NE && l <= neurons) || l >= neurons + dNE) {
            lastSpike[l] = 2 * tRef_I[l];
        }
    }
}

// Here, we initialize our initial voltage.
void initialize_InitialVoltage(int neurons, int dNeurons, vector< vector<double> >& V, vector<double>& E) {
    // Initial Voltage is equal to the leak
    for (int j = 0; j < neurons + dNeurons; j++) {
        V[j][0] = E[j];
    }
}

// Here, we initialize our facilitation/depression factors.
void initialize_FacAndDep(bool facilitation, bool depression, int NE, vector<double>& Fac, vector<double>& Dep, double F0, double D0) {
    if (facilitation) {
        for (int j = 0; j < NE; j++) {
            Fac[j] = F0;
        }
    }
    if (depression) {
        for (int j = 0; j < NE; j++) {
            Dep[j] = D0;
        }
    }
}

// Here, we set our stim value for each trial.
void stimSetup(int& stim, int trial, vector< int>& stimShuffle) {
    if (trial == 1 || trial % 4 == 1) {
        random_shuffle(stimShuffle.begin(), stimShuffle.end());
    }
    cout << "Trial block set" << endl;
    for (int i = 0; i < inputs; i++) {
        cout << stimShuffle[i] << endl;
    }
    if (trial % 4 == 1) {
        stim = stimShuffle[0];
    }
    if (trial % 4 == 2) {
        stim = stimShuffle[1];
    }
    if (trial % 4 == 3) {
        stim = stimShuffle[2];
    }
    if (trial % 4 == 0) {
        stim = stimShuffle[3];
    }
    cout << "stim: " << stim << endl;
}

// Here, we account for AMPADecay in pSyn.
void pSyn_AMPADecay(int nIn, vector<double>& pSyn, double AMPADecay) {
    for (int cell = 0; cell < nIn; cell++) {
        pSyn[cell] *= AMPADecay;
    }
}

// Here, we initialize cellStart and cellStop (and cellStart2 and cellStop2).
void initialize_CellStartAndStop(double ti, double tOn0, double tOff0, double tOn1, double tOff1, bool& InputOn, bool& InputOn2, int stim, int& cellStart, int& cellStop, int& cellStart2, int& cellStop2, int inputCells, bool urgency, double& g_Urgency, double gMax_Urgency) {
    if (ti > tOn0 && ti <= tOff0) {
        InputOn = 1;
        if (stim == 1 || stim == 2) { //A
            cellStart = 0;
            cellStop = inputCells - 1;				  
        } else { // C
            cellStart = inputCells;
            cellStop = 2 * inputCells - 1;
        }
    }
    if (ti > tOn1 && ti <= tOff1) {
        InputOn2 = 1; 
        if (stim == 1 || stim == 3) { // B
            cellStart2 = 2 * inputCells;
            cellStop2 = 3 * inputCells - 1;
        } else { // D
            cellStart2 = 3 * inputCells;
            cellStop2 = 4 * inputCells - 1;
        }				
        if (urgency) { 
            g_Urgency = gMax_Urgency * (ti - tOn1) / (tOff1 - tOn1);
            if (g_Urgency < 0) { 
                g_Urgency = 0;
            }
        }
    }
}

// Here, we perform the calculations necessary for a given input.
void use_InputOn(int cellStart, int cellStop, double dt, MTRand& rand3, double r0, vector< int>& nPSpikes, int max_InputSpikes, double ti, vector< vector<double> >& pSpikes, double pSynMax, double s0, vector<double>& pSyn) {
    for (int cell = cellStart; cell <= cellStop; cell++) { // first stimulus of pair
        double dt1 = dt;
        double randt = rand3() / r0;
        while (randt < dt1) {
            nPSpikes[cell]++;
            if (nPSpikes[cell] >= max_InputSpikes) {
                cout << "cell " << cell << " input spikes " << nPSpikes[cell] << " r0 " << r0 << " t " << ti << endl;
            }
            pSpikes[cell][nPSpikes[cell]] = ti;     
            pSynMax = s0 * (1 - pSyn[cell]); 
            pSyn[cell] += pSynMax;
            dt1 -= randt;
            randt = rand3() / r0;
        }
	}
}

// Here, we perform our synapse conductance update.
void synapse_Cond_Update(int neurons, int dNeurons, int j, bool C1[][totalNeurons], int NE, int dNE, vector<double>& gSyn_I_Vec, vector< vector<double> >& W, vector<double>& syn_I, vector<double>& syn_E, bool facilitation, bool depression, vector<double>& gSyn_E_Vec, vector<double>& gSyn_NMDA, vector<double>& syn_NMDA, vector<double>& synE_Fac, vector<double>& synNMDA_Fac, vector<double>& synE_Dep, vector<double>& synNMDA_Dep) {
    for (int k = 0; k < neurons + dNeurons; k++) {  // maybe run this loop through the dNE's
        if (j != k && C1[k][j] == 1 && W[k][j] != 0) { // Inhibitory update
            if ((k >= NE && k < neurons) || k >= neurons + dNE) {  // Inhibitory cells
                gSyn_I_Vec[j] += W[k][j] * syn_I[k];  //gGABA update
            } 
            // Associative layer ==> Decision layer Static/facilitating/depressing synapses 
            if (k < NE && j >= neurons && j < neurons + dNE) {
                if (!facilitation && !depression) {
                    gSyn_E_Vec[j] += W[k][j] * syn_E[k];
                    gSyn_NMDA[j] += W[k][j] * syn_NMDA[k];
                }
                if (facilitation) {
                    gSyn_E_Vec[j] += W[k][j] * synE_Fac[k];
                    gSyn_NMDA[j] += W[k][j] * synNMDA_Fac[k];
                }
                if (depression) { // HAXED TO MAKE AL-TO-DL Facilitating Facilitating for this Simulation only!!!!!
                    gSyn_E_Vec[j] += W[k][j] * synE_Fac[k];
                    gSyn_NMDA[j] += W[k][j] * synNMDA_Fac[k];
                }  
            } else { // All other synapses      
                if (k < NE && j < NE) { // Facilitating recurrent excitatory synapses.
                    gSyn_E_Vec[j] += W[k][j] * synE_Dep[k];
                    gSyn_NMDA[j] += W[k][j] * synNMDA_Dep[k];
                } else {
                    gSyn_E_Vec[j] += W[k][j] * syn_E[k];
                    gSyn_NMDA[j] += W[k][j] * syn_NMDA[k];
                }
            }
        }
	} 
}

// Here, we perform our AL/DL conductance update.
void ALDL_Cond_Update(int j, vector<double>& gSyn_I_Vec, double gSyn_I, int NE, vector<double>& gSyn_E_Vec, double GAL_E_AMPA, vector<double>& gSyn_NMDA, double GAL_E_NMDA, int neurons, double GAL_I_AMPA, double GAL_I_NMDA, int dNE, double GDL_AMPA, double GDL_NMDA, double nSigma_I, double nSigma_E, MTRand& rand2, double dt, double sigma_I, double sigma_E, double ti, double tOn1) {
    gSyn_I_Vec[j] *= gSyn_I;  // Update GABA conductances
    // Update the Associative layer AMPA/NMDA conductances
	if (j < NE) { // Excitatory
        gSyn_E_Vec[j] *= GAL_E_AMPA;
        gSyn_NMDA[j] *= GAL_E_NMDA;
	}
	if (j >= NE && j < neurons) { // Inhibitory
        gSyn_E_Vec[j] *= GAL_I_AMPA;
        gSyn_NMDA[j] *= GAL_I_NMDA;
	}
	if (j >= neurons && j < neurons + dNE) {  // Update the Decision layer AMPA/NMDA conductances
        gSyn_E_Vec[j] *= GDL_AMPA; // Excitatory DL conductances
        gSyn_NMDA[j] *= GDL_NMDA;
	} else if (j >= neurons + dNE) { 
            gSyn_E_Vec[j] *= GDL_I_AMPA;  // Inhibitory DL conductances
            gSyn_NMDA[j] *= GDL_I_NMDA;
        }	  
	if (j < neurons) {
        gSyn_I_Vec[j] += nSigma_I * rand2() * sqrt(dt);  // inhibitory conductance noise term
        gSyn_E_Vec[j] += nSigma_E * rand2() * sqrt(dt);  // excitatory conductance noise term
	} else {
        gSyn_I_Vec[j] += sigma_I * rand2() * sqrt(dt);  // inhibitory conductance noise term
        gSyn_E_Vec[j] += sigma_E * rand2() * sqrt(dt);  // excitatory conductance noise term
	}
	if (j >= neurons && j < neurons + dNE && ti < tOn1) {  // Set to zero before presentation of 2nd stimulus. Only decide at 2nd stim
        gSyn_E_Vec[j] *= 0; // Excitatory DL conductances
        gSyn_NMDA[j] *= 0;
	}    
}

// Here, we update gpSyn_Vec.
void update_gpSyn_Vec(int j, int nIn, vector<double>& gpSyn_Vec, vector< vector<double> >& pW, double gpSyn, vector<double>& pSyn, bool urgency, int neurons, int dNE, double g_Urgency) {
    gpSyn_Vec[j] = 0;
	for (int stim = 0; stim < nIn; stim++) {
        if (pW[stim][j] != 0) {
            gpSyn_Vec[j] += gpSyn * pW[stim][j] * pSyn[stim]; // A, X, B, Y Poisson Input currents
        }
    }
	if (urgency) { // Implement ramping urgency signal here
        if  (j >= neurons && j < neurons + dNE) { // to E-cells in DL only
            gpSyn_Vec[j] += g_Urgency;
        }
	} 
}

// Here, we update our refractory decay.
void update_Refract_Cond(int j, int NE, int neurons, int dNE, vector<double>& g_Ref, double dt, double tRef_Ij, double tRef_Ej) {
    if ((j >= NE && j < neurons) || j >= neurons + dNE) {  // Inhibitory cells
        g_Ref[j] *= exp(-dt / tRef_Ij); // Inhibitory refractory decay
	} else {
        g_Ref[j] *= exp(-dt / tRef_Ej); // Excitatory refractory decay.
	}
}

// The gaussian function performs the Gaussian operations required to compute vInf.
double gaussian() {
    static int iset = 0;
    static double gset;
    double fac, rsq, v1, v2;
    if (iset == 0) {
        do {
            v1 = 2.0 * rand2() - 1;
            v2 = 2.0 * rand2() - 1;
            rsq = v1 * v1 + v2 * v2;
        } while (rsq >= 1 || rsq == 0);
        fac = sqrt(-2.0 * log(rsq) / rsq);
        gset = v1 * fac;
        iset = 1;
        return v2 * fac;
    } else {
        iset = 0;
        return gset;
    }
}

// Here, we update vInf.
void update_vInf(vector<double>& vInf, int j, int neurons, double g_Lj, double Ej, double g_Refj, double vReset, double gSyn_NMDAj, double g_NMDAj, double E_NMDA, double gSyn_E_Vecj, double E_AMPA, double gSyn_I_Vecj, double E_GABA, double gpSyn_Vecj, double nSigma, double sigma, double dt, double g_Tot) {
    if (j < neurons) {  // Associative layer uses nSigma
        vInf[j] = (g_Lj * Ej + g_Refj * vReset + gSyn_NMDAj * g_NMDAj * E_NMDA + gSyn_E_Vecj * E_AMPA + gSyn_I_Vecj * E_GABA + gpSyn_Vecj * E_AMPA + nSigma * gaussian() * sqrt(dt)) / g_Tot; // Vss
	} else {
        vInf[j] = (g_Lj * Ej + g_Refj * vReset + gSyn_NMDAj * g_NMDAj * E_NMDA + gSyn_E_Vecj * E_AMPA + gSyn_I_Vecj * E_GABA + gpSyn_Vecj * E_AMPA + sigma * gaussian() * sqrt(dt)) / g_Tot; // Vss
    }
}

// Here, we update our dynamic thresholds.
void update_Dynamic_Thresholds(int j, int NE, int neurons, int dNE, double& vth0, double vth_Ij, double vth_Ej, double& tRef, double tRef_Ij, double tRef_Ej, vector<double>& Vth, double vth_Max, double ti, double lastspikej) {
    if ((j >= NE && j < neurons) || j >= neurons + dNE) {  // Inhibitory cells
        vth0 = vth_Ij;
        tRef = tRef_Ij;
        Vth[j] = vth0 + (vth_Max - vth0) * exp(-(ti - lastspikej) / tRef); // Impossible Vth to counter LIF reset
	} else {
        vth0 = vth_Ej;
        tRef = tRef_Ej;
        Vth[j] = vth0 + (vth_Max - vth0) * exp(-(ti - lastspikej) / tRef); // Impossible Vth to counter LIF reset
	}
}

// Here, we perform our inhibitory refractory/reset.
void inhib_Excit_Refract_Reset(int j, int i, int NE, int neurons, int dNE, vector< vector<double> >& V, double vSpike, double vReset_Ij, double vReset_Ej, double Vthj, vector<double>& g_Ref, double delta_GRef, vector<double>& lastSpike, double ti) {
    if ((j >= NE && j < neurons) || j >= neurons + dNE) { // Inhibitory refractory/reset
        if (V[j][i-1] == vSpike) {  
            V[j][i] = vReset_Ij;
        }
        if (V[j][i] >= Vthj) {
            V[j][i] = vSpike;
            g_Ref[j] += delta_GRef;
            lastSpike[j] = ti;
        }
	} else { // Excitatory refractory/reset
        if (V[j][i-1] == vSpike) {
            V[j][i] = vReset_Ej;
        }
        if (V[j][i] >= Vthj) {
            V[j][i] = vSpike;
            g_Ref[j] += delta_GRef;
            lastSpike[j] = ti;
        }
	}
}
// Here, we record our spike times.
void record_spikeTimes(int j, double Vji, double vSpike, vector<int>& numberSpikes, int someIndex, vector< vector<double> >& spikeTimes, double ti) {
    if (Vji == vSpike) {
        numberSpikes[j]++;
        someIndex = numberSpikes[j];
        spikeTimes[j][someIndex] = ti;
	}   
}

// Here, we update our synaptic conductance (GABA, AMPA, NMDA).
void update_Syn_Cond(int j, vector<double>& syn_I, vector<double>& syn_E, double GABADecay, double AMPADecay, int neurons, vector<double>& syn_NMDA, double AL_NMDADecay, double DL_NMDADecay, vector<double>& synE_Fac, vector<double>& synNMDA_Fac, vector<double>& synE_Dep, vector<double>& synNMDA_Dep) {
	syn_I[j] *= GABADecay;
	syn_E[j] *= AMPADecay;
	if (j < neurons) {
        syn_NMDA[j] *= AL_NMDADecay;
	} else {
        syn_NMDA[j] *= DL_NMDADecay;
	}
	synE_Fac[j] *= AMPADecay;
	synNMDA_Fac[j] *= AL_NMDADecay;
	synE_Dep[j] *= AMPADecay;
	synNMDA_Dep[j] *= AL_NMDADecay;
}

// Here, we update our facilitation and depression decays.
void update_Fac_Dep_Decays(bool facilitation, bool depression, int j, int NE, vector<double>& Fac, vector<double>& Dep, double dt, double tau_Fac, double tau_Dep) {
    if (facilitation) {
        if (j < NE) {
            Fac[j] = 1+ (Fac[j] - 1) * exp(-dt / tau_Fac);
        }
    }
	if (depression) {
        if (j < NE) {
            Dep[j] = 1 - (1 - Dep[j]) * exp(-dt / tau_Dep);
        }
    }
}

// Here, we update our inhibitory and excitatory synapses.
void update_Inhib_Excit_Synapses(int j, double Vji, double vSpike, double synMax_I, double synMax_E, double synMax_NMDA, double s0, vector<double>& syn_I, vector<double>& syn_E, vector<double>& syn_NMDA, int NE, int neurons, int dNE, bool facilitation, bool depression, vector<double>& synE_Fac, vector<double>& synNMDA_Fac, vector<double>& Fac, vector<double>& synE_Dep, vector<double>& synNMDA_Dep, vector<double>& Dep, double alpha, double fMax, double dFrac) {
    if (Vji == vSpike) {
        synMax_I = s0 * (1-syn_I[j]);
        synMax_E = s0 * (1-syn_E[j]);
        synMax_NMDA = s0 * (1-syn_NMDA[j]);
        if ((j >= NE && j < neurons) || j >= neurons + dNE) { // Inhibitory cells
            syn_I[j] += synMax_I;
        } else {  // All excitatory neurons
            syn_E[j] += synMax_E;	  
            syn_NMDA[j] += synMax_NMDA;
            if (j < NE) {  // Only excitatory neurons in Associative layer Facilitate/Depress
                if (facilitation) {
                    synE_Fac[j] += s0 * (1 - synE_Fac[j]) * Fac[j];
                    synNMDA_Fac[j] += s0 * (1 - synNMDA_Fac[j]) * Fac[j];
                }
                if (depression) {
                    synE_Dep[j] += s0 * (1 - synE_Dep[j]) * Dep[j];
                    synNMDA_Dep[j] += s0 * (1 - synNMDA_Dep[j]) * Dep[j];
                }
            }	
        } 
        // Update facilitation & depression factors after spike
        if (facilitation) {
            if (j < NE) {
                Fac[j] += alpha * (fMax - Fac[j]);
            }
        }
        if (depression) {
            if (j < NE) {
                Dep[j] *= (1 - dFrac);
            }
        }
    }
}

// Here, we update our cue firing rate in cueFR.
void write_Cue_Firing_Rate(int neurons, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes, double tOn0, double tOff0, double tOn1, double tOff1, vector< int>& cueFR, double cue1Length, double cue2Length) {
    for (int j = 0; j < neurons; j++) {
        for (int k = 0; k < numberSpikes[j]; k++) {
            if ((spikeTimes[j][k] > tOn0 && spikeTimes[j][k] < tOff0) || (spikeTimes[j][k] > tOn1 && spikeTimes[j][k] < tOff1)) { // all spikes from start of first stim to end of 2nd stim.
                cueFR[j]++;
            }
        }
        cueFR[j] /= (cue1Length+cue2Length);
    }
}

// Here, we measure selectivity in 1st stim.
void measure_Selectivity_1st_Stim(int neurons, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes, double tOn0, double tOff0, double tOn1, double tOff1, vector< int>& cueFR_1stStim, double cue1Length) {
    for (int j = 0; j < neurons; j++) {
        for (int k = 0; k < numberSpikes[j]; k++) {
            if (spikeTimes[j][k] > tOn0 && spikeTimes[j][k] < tOff0) { // 2nd stim spikes only.
                cueFR_1stStim[j]++;
            }
        }
    cueFR_1stStim[j] /= cue1Length;
    }
}

// Here, we measure selectivity in 2nd stim (persistence).
void measure_Selectivity_2nd_StimPerst(int neurons, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes, double tOn1, double tOff1, vector< int>& cueFR_2ndStim_Perst, double cue2Length) {
    for (int j = 0; j < neurons; j++) {
        for (int k = 0; k < numberSpikes[j]; k++) {
            if (spikeTimes[j][k] > tOn1 && spikeTimes[j][k] < tOff1) { // 2nd stim spikes only.
                cueFR_2ndStim_Perst[j]++;
            }
        }
        cueFR_2ndStim_Perst[j] /= cue2Length;
    }
}

// Here, we measure selecitvity in persistence.
void measure_Selectivity_Perst(int NE, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes, double tOff0, double tOn1, vector< int>& cueFR_Perst, double delay1Length) {
    for (int j = 0; j < NE; j++) {
        for (int k = 0; k < numberSpikes[j]; k++) {
            if (spikeTimes[j][k] > tOff0 && spikeTimes[j][k] < tOn1) {// 2nd stim spikes only.
                cueFR_Perst[j]++;
            }
        }
        cueFR_Perst[j] /= delay1Length;
    }
}

// Here, we find the mean excitatory persistent delay activity.
void find_Mean_ENeuron_Delay_Rate(double& mean_ENeuron_Delay_Rate, int NE, vector< int>& cueFR_Perst) {                   
    for (int j = 0; j < NE; j++) {
        mean_ENeuron_Delay_Rate += cueFR_Perst[j];
    }
    mean_ENeuron_Delay_Rate /= NE;
    cout << "Mean excitatory persistent delay activity: " << mean_ENeuron_Delay_Rate << endl;
}
        
// Here, we measure the selectivity of 2nd perst.
void measure_Selectivity_2nd_Perst(int NE, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes, double tOff1, double tMax, vector< int>& cueFR_2nd_Perst, double delay2Length) {
    for (int j = 0; j < NE; j++) {
        for (int k = 0; k < numberSpikes[j]; k++) {
            if (spikeTimes[j][k] > tOff1 && spikeTimes[j][k] < tMax) { // 2nd delay spikes only.
                cueFR_2nd_Perst[j]++;
            }
        }
        cueFR_2nd_Perst[j] /= delay2Length;
    }
}

// Here, we find the mean excitatory 2nd persistent delay activity.
void find_Mean_ENeuron_2nd_Delay_Rate(double& mean_ENeuron_2nd_Delay_Rate, int NE, vector< int>& cueFR_2nd_Perst) {
    for (int j = 0; j < NE; j++) {
        mean_ENeuron_2nd_Delay_Rate += cueFR_2nd_Perst[j];
    }
    mean_ENeuron_2nd_Delay_Rate /= NE;
    cout << "Mean excitatory 2nd persistent delay activity: " << mean_ENeuron_2nd_Delay_Rate << endl;
}

// Here, we find the bin rate.
void find_Bin_Rate(int neurons, int nBins, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes, double bindT, vector< vector<double> >& AL_ECell_BinRate/*double AL_ECell_BinRate[][19]*/) {    
    for (int cell = 0; cell < neurons; cell++) {
        for (int ispike = 0; ispike < numberSpikes[cell]; ispike++) {
            int bin = int(spikeTimes[cell][ispike] / bindT);
            AL_ECell_BinRate[cell][bin]++;
        }
        for (int bin = 0; bin < nBins; bin++) {
            AL_ECell_BinRate[cell][bin] /= bindT; // Normalizes meanFR by bin size
        }
    }
}

// Here, we find mean_WR_ds.
void output_TrialRate_ALNeurons_DLWeights(int NE, int neurons, int dNEPool, vector<double>& WR_ds, vector<double>& Wds, double& mean_WR_ds, vector< vector<double> >& W, vector< int>& cueFR) {
    for (int j = 0; j < NE; j++) {
        Wds[j] = W[j][neurons + 1] - W[j][neurons + dNEPool+1];
        WR_ds[j] = cueFR[j] * (W[j][neurons + 1] - W[j][neurons + dNEPool + 1]);
    }    
    for (int j = 0; j < NE; j++) {
        mean_WR_ds += WR_ds[j];
    }
    mean_WR_ds /= NE;
    cout << "Mean cueFR* (Wdig-Wswitch) = " << mean_WR_ds << endl;
}

// Here, we find the mean mean firing rate.
void find_Mean_Mean_Firing_Rate(int neurons, int dNeurons, vector< int>& numberSpikes, double tMax, vector<double>& meanFiringRate, vector<double>& mean_Mean_Firing_Rate, double& mean_rE, double& mean_rI, int NE, int NI) {
    for (int j = 0; j < neurons + dNeurons; j++) {
        meanFiringRate[j] = double(numberSpikes[j]) / tMax;
        mean_Mean_Firing_Rate[j] += meanFiringRate[j];
    }
    for (int j = 0; j < NE; j++) {
        mean_rE += meanFiringRate[j];
    }
    mean_rE /= NE;
    for (int j = NE; j < NE + NI; j++) {
        mean_rI += meanFiringRate[j];
    }
    mean_rI /= NI;
}

// Here, we divide the mean mean firing rate by inputs.
void divide_Mean_Mean_Firing_Rate_By_Inputs(int neurons, int dNeurons, bool oneTrialHomeo, vector<double>& mean_Mean_Firing_Rate, int inputs) {
    for (int j = 0; j < neurons + dNeurons; j++) {
        if (!oneTrialHomeo) {
            mean_Mean_Firing_Rate[j] /= double(inputs); // meanFR = meanFR/#inputs
        }
    }
}

// Here, we find the excitatory input & recurrent excitatory homeostatic shift / update.
void input_And_Recurrent_Excitatory_Homeostatic(int NE, vector<double>& W_Shift, double e_E, vector<double>& r_Goal, vector<double>& mean_Mean_Firing_Rate, bool inputEHomeostasis, int nIn, vector< vector<double> >& pW, bool EE_Homeostasis, vector< vector<double> >& W, double& sumWShift) {
    for (int j = 0; j < NE; j++) {
        W_Shift[j] = e_E * (r_Goal[j] - mean_Mean_Firing_Rate[j]); // Excitatory Input & Recurrent Excitatory homeostatic shift
        mean_Mean_Firing_Rate[j] = 0;
        if (inputEHomeostasis) {
            for (int i = 0; i < nIn; i++) {
                pW[i][j] *= (1 + W_Shift[j]);  // Input excitatory homeostatic update
            }
        }
        if (EE_Homeostasis) {
            for (int i = 0; i < NE; i++) {
                W[i][j] *= (1 + W_Shift[j]); // Recurrent Excitatory homeostatic update
            }
        }
        sumWShift += W_Shift[j];  // measure the total homeostasis shift/change to Input->E & WEE
    }
    cout << "Mean E Wshift " << sumWShift / double(NE) << endl;
}

// Here, we find the inhibitory input homeostatic shift.
void input_Inhibitory_Homeostatic(int NE, int NI, vector<double>& W_Shift, double e_Input_I, vector<double>& r_Goal, vector<double>& mean_Mean_Firing_Rate, int nIn, vector< vector<double> >& pW, double& sumWIShift) {                                 
    for (int j = NE; j < NE + NI; j++) {
        W_Shift[j] = e_Input_I * (r_Goal[j] - mean_Mean_Firing_Rate[j]);  // Inhibitory Input homeostatic shift
        mean_Mean_Firing_Rate[j] = 0;  // initialize meanFR to zero.
        for (int i = 0; i < nIn; i++) {
            pW[i][j] *= (1 + W_Shift[j]);  // Inhibitory Input homeostatic shift
        }
        sumWIShift += W_Shift[j];
    }
	cout << "Mean I Wshift " << sumWIShift / double(NI) << endl;
}

// Here, we find the decision layer homeostatic shift.
void input_DL_Homeostatic(int neurons, int dNE, vector<double>& DW_Shift, double eD, double rgD, vector<double>& mean_Mean_Firing_Rate, int NE, vector< vector<double> >& W, double& meanDWShift) {
    for (int j = neurons; j < neurons + dNE; j++) {
        DW_Shift[j] = eD * (rgD - mean_Mean_Firing_Rate[j]);
        for (int i = 0; i < NE; i++) {
             W[i][j] *= (1 + DW_Shift[j]);
        } 
        meanDWShift += DW_Shift[j];
        mean_Mean_Firing_Rate[j] = 0;
    }
    cout << "Mean DWshift " << meanDWShift / double(dNE) << endl;
}

// Here, we implement LTPi homeostasis.
void LTPI_homeostasis(int NE, vector<double>& W_Shift_I, double e_I, vector<double>& r_Goal_I_to_E, vector<double>& meanFiringRate, int NI, vector< vector<double> >& W, double& sum_WShift_I) {
    for (int j = 0; j < NE; j++) {
        W_Shift_I[j] = e_I * (r_Goal_I_to_E[j] - meanFiringRate[j]);
        if (W_Shift_I[j] < 0) {
            W_Shift_I[j] = 0;
        }
        for (int i = 0; i < NI; i++) {
            W[NE+i][j] *= (1 - W_Shift_I[j]);
        }
        sum_WShift_I += W_Shift_I[j];
    }
    cout << "Mean W_Shift_I " << sum_WShift_I/double(NE) << endl;
}

// Here, we implement the PG rule for EESTDP.
void PG_Rule_EESTDP(int NE, bool C1[][totalNeurons], vector< int>& counter, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes, double tau3m, int endFlag, double deltaT, vector< vector<double> >& dwSum, double tau2m, double A2m, double A3m, double tau2p, double tau3p, double A2p, double A3p, double dW, vector< vector<double> >& W, double W0, double& sum_PGdw, int& countIJ) {
    for (int i = 0; i < NE; i++) {
        for (int j = 0; j < NE; j++) {
            if (i != j) { // no self-self plasticity Wii, Wjj
                if (C1[i][j] == 1) {     // Only bother doing this if there is a connection.
                    counter[j] = 1;
                    for (int ispike = 0; ispike < numberSpikes[i]; ispike++) {
                        double ispiketime = spikeTimes[i][ispike];
                        double ifactm = 0;
                        for (int iterm = 1; iterm < ispike; iterm++) {
                            ifactm += exp(-(ispiketime - spikeTimes[i][ispike-iterm]) / tau3m);
                        }
                        int jspike = counter[j];
                        endFlag= 0;
                        while (spikeTimes[j][jspike] < ispiketime && endFlag == 0) { // post-pre
                            deltaT = spikeTimes[j][jspike] - ispiketime;
                            dwSum[i][j] -= (exp(deltaT / tau2m)) * (A2m + A3m * ifactm);
                            if (jspike < numberSpikes[j]) {
                                jspike += 1;
                            } else {
                                endFlag = 1;
                            }
                        }
                        while ((spikeTimes[j][jspike] < ispiketime + 5*tau2p ) && (endFlag== 0)) { // pre-post, replaced 100ms with 5*tau2p for scaling, exp(-5)
                            deltaT = spikeTimes[j][jspike] - ispiketime;
                            double jfactp = 0;
                            for (int jterm = 1; jterm < jspike; jterm++) {
                                jfactp += exp(-(spikeTimes[j][jspike] - spikeTimes[j][jspike-jterm]) / tau3p);
                            }
                            dwSum[i][j] += exp(-deltaT / tau2p) * (A2p + A3p * jfactp);
                            if (jspike < numberSpikes[j]) { 
                                jspike += 1;
                            } else {
                                endFlag = 1;
                            }
                        }
                    } 
                    dwSum[i][j] *= dW;
                    if (dwSum[i][j] > 0.5) {
                        dwSum[i][j] = 0.5;
                    } else {
                        if (dwSum[i][j] < -0.5) { 
                            dwSum[i][j] = -0.5;
                        }
                    }	      
                    W[i][j] += dwSum[i][j] * W0;
                    sum_PGdw += dwSum[i][j]; // total change from triplet STDP across synapses
                    countIJ++; // number of synapses updated
	                if (W[i][j] < 0) {
                        W[i][j] = 0;
                    }
                    if (W[i][j] > 20 * W0) {    
                        W[i][j] = 20 * W0;  // abitrary ceiling on Wij (final), currently 5 or  20 x W0,
                    }             
                } 
            } 
        } 
    }
	cout << "Mean PG wEE dwSum " << sum_PGdw / double(countIJ) << endl;
}

// Here, we implement the PG rule for input ESTDP.
void PG_Rule_Input_ESTDP(int NE, int nIn, vector< vector<double> >& pW, vector< int>& counter, vector< vector<double> >& spikeTimes, vector< int>& numberSpikes, vector< int>& nPSpikes, vector< vector<double> >& pSpikes, double tau3m, int endFlag, double deltaT, vector< vector<double> >& dwEinsum, double tau2m, double A2m, double A3m, double tau2p, double tau3p, double A2p, double A3p, double pW02, double sum_InputPGdw, int countIJ) {
    for (int j = 0; j < NE; j++) {
        for (int i = 0; i < nIn; i++) {
            if (pW[i][j] != 0) {     // Cell specific input bool since random prob per cell.
                counter[j] = 1;
                for (int ispike = 0; ispike < nPSpikes[i]; ispike++) {
                    double ispiketime = pSpikes[i][ispike]; 
                    double ifactm = 0;
                    for (int iterm = 1; iterm < ispike; iterm++) {
                        ifactm += exp(-(ispiketime - pSpikes[i][ispike-iterm]) / tau3m);
                    }
                    int jspike = counter[j];
                    endFlag= 0;
                    while (spikeTimes[j][jspike] < ispiketime && endFlag == 0) { // post-pre
                        deltaT = spikeTimes[j][jspike] - ispiketime;
                        dwEinsum[i][j] -= exp(deltaT / tau2m) * (A2m + A3m * ifactm);
                        if (jspike < numberSpikes[j]) {
                            jspike += 1;
                        } else {
                            endFlag = 1;
                        }
                    }
                    while (spikeTimes[j][jspike] < ispiketime + 5 * tau2p && endFlag== 0) { // pre-post, replaced 100ms with 5*tau2p for scaling, exp(-5)
                        deltaT = spikeTimes[j][jspike] - ispiketime;
                        double jfactp = 0;
                        for (int jterm = 1; jterm < jspike; jterm++) {
                            jfactp += exp(-(spikeTimes[j][jspike] - spikeTimes[j][jspike-jterm]) / tau3p);
                        }
                        dwEinsum[i][j] += exp(-deltaT / tau2p) * (A2p + A3p * jfactp);
                        if (jspike < numberSpikes[j]) { 
                            jspike += 1;
                        } else {
                            endFlag = 1;
                        }
                    }
                } 
                dwEinsum[i][j] *= dW;	
                if (dwEinsum[i][j] > .5) {
                    dwEinsum[i][j] = 0.5;
                } else {
                    if (dwEinsum[i][j] < -0.5) {
                        dwEinsum[i][j] = -0.5;
                    }
                }
                pW[i][j] += dwEinsum[i][j] * pW02;
                if (pW[i][j] < 0) {
                    pW[i][j] = 0;
                }
                if (pW[i][j] > 20 * pW02) {
                    pW[i][j] = 20 * pW02;
                }
                sum_InputPGdw += dwEinsum[i][j];
                countIJ++;
            } 
        } 
	} 
    cout << "Mean PG Input-E dwSum " << sum_InputPGdw / double(countIJ) << endl;
 }

// Here, we implement the LTPi rule.
void LTPi_Rule(int NE, int neurons, bool C1[][totalNeurons], vector< int>& numberSpikes, vector< vector<double> >& spikeTimes, double dt, bool VLTPi, vector< vector<double> >& V, double VLTPi_Threshold, double LTPi_Window, vector< vector<double> >& dwSum_I, double idW, vector< vector<double> >& W, double wIE0, double& sum_I_dwSum, int countIJ) {
    for (int i = NE; i < neurons; i++) { // NI
        for (int j = 0; j < NE; j++) { // NE
            if (C1[i][j] == 1) {
                int jspike = 0;
                for (int k = 0; k < numberSpikes[i]; k++) { // # ispikes per NI
                    double itime = spikeTimes[i][k];  // inhibitory spike times, in seconds
                    int itimebin = spikeTimes[i][k] / dt; // inhibitory spike time bins, for accessing correct voltage.
                    double Vfrac = 1;
                    if (VLTPi) {
                        if (V[j][itimebin] < VLTPi_Threshold) { // post-synaptic excitatory cell voltage at time of pre-synaptic inhibitory spike.
                            Vfrac = 0;  // Do Not potentiate if post-synaptic excitatory cell is below a voltage threshold, usually -65mV (i.e. 5mV above Vleak)
                        }
                    }
                    while (jspike < numberSpikes[j] && spikeTimes[j][jspike] < itime-LTPi_Window) {
                        jspike++;
                    }
                    if ((spikeTimes[j][jspike] > itime + LTPi_Window) || (jspike == numberSpikes[j] && spikeTimes[j][jspike] < itime-LTPi_Window)) {
                        dwSum_I[i-NE][j] += idW * Vfrac;
                    }        
                }
                if (dwSum_I[i-NE][j] < 0) {
                    dwSum_I[i-NE][j] = 0; // min 0
                }
                if (dwSum_I[i-NE][j] > 0.5) {
                    dwSum_I[i-NE][j] = 0.5; // max .5
                }    
                W[i][j] += dwSum_I[i-NE][j] * wIE0; 
                sum_I_dwSum += dwSum_I[i-NE][j];
                countIJ++;
                if (W[i][j] < 0) {
                    W[i][j] = 0; // Wij cannot be negative
                }    
                if (W[i][j] > 20 * wIE0) {
                    W[i][j] = 20 * wIE0;  // arbitrary ceiling on Wij (final), currently 5 or  20 x W0,
                }
            }
        }
    }
    cout << "Mean dwSum_I " << sum_I_dwSum / double(countIJ) << endl;
}

// Here, we make sure to floor / ceiling any values of pW and W outside certain bounds.
void floor_And_Ceiling(int nIn, int neurons, vector< vector<double> >& pW, vector< vector<double> >& W) {
    for (int i = 0; i < nIn; i++) {
        for (int j = 0; j < neurons; j++) {
            if (pW[i][j] > 20 * pW0) { 
                pW[i][j] = 20 * pW0;
            }
			if (pW[i][j] < 0) {
				pW[i][j]  = 0;
			}
        }
   }
    for (int i = 0; i < neurons; i++) {
        for (int j = 0; j < neurons; j++) {
            if (W[i][j] > 20 * W0) {
                W[i][j] = 20 * W0;
            }
            if (W[i][j] < 0) {
                W[i][j]  = 0;
            }
        }
    }
}

// Here, we perform the first part of our structural plasticity implementation.
void structPlast_Operation_1(int NE, bool C1[][totalNeurons], vector< int>& counter, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes, int endFlag, double tauStruct, double deltaT, vector< vector<double> >& dWstruct, double structRateFact, double& mean_dWStruct, int& n_dW_Struct, double& max_dWStruct, double& min_dWStruct) {
    for (int i = 0; i < NE; i++) {
        for (int j = 0; j < NE; j++) {
            if ( C1[i][j] == 1) {     // Only bother doing this if there is a connection.
                counter[j] = 1;
                for (int ispike = 0; ispike<numberSpikes[i]; ispike++) {
                    double ispiketime = spikeTimes[i][ispike];
                    int jspike = counter[j];
                    endFlag = 0;
                    while (spikeTimes[j][jspike] < ispiketime && endFlag == 0) { // post-pre
                        if (jspike < numberSpikes[j]) {
                            jspike += 1;
                        } else {
                            endFlag = 1;
                        }
                    }
                    while (spikeTimes[j][jspike] < ispiketime + 5 * tauStruct && endFlag == 0) { // pre-post
                        deltaT = spikeTimes[j][jspike] - ispiketime;
                        dWstruct[i][j] += exp(-deltaT / tauStruct);
                        if (jspike < numberSpikes[j]) {
                            jspike += 1;
                        } else {
                            endFlag = 1;
                        }
                    }
                } 
                dWstruct[i][j] -= structRateFact;
                mean_dWStruct += dWstruct[i][j];
                n_dW_Struct++;
                if (dWstruct[i][j] > max_dWStruct) {
                    max_dWStruct = dWstruct[i][j];
                }
                if (dWstruct[i][j] < min_dWStruct) { 
                    min_dWStruct = dWstruct[i][j];
                }
            } 
        } 
    } 
    cout << "Mean dWstruct " << mean_dWStruct / n_dW_Struct << " N " << n_dW_Struct << " max " << max_dWStruct << " min " << min_dWStruct << endl;
}

// Here, we perform the second part of our structural plasticity implementation.
void structPlast_Operation_2(int NE, bool C1[][totalNeurons], vector< vector<double> >& dWstruct, double removeStrength, vector< int>& nPreChange, vector< int>& nPostChange, int& nTotChanges, vector< vector<double> >& W, MTRand& rand1, double W0, double ranWStrengths, double structMemoryFact, double structMax, double& mean_dWStruct, int& n_dW_Struct, double& max_dWStruct, double& min_dWStruct) {
    for (int i = 0; i < NE; i++) {
        for (int j = 0; j < NE; j++) {
            if (C1[i][j] == 1) {     // Only bother doing this if there is a connection.
                if (dWstruct[i][j] < removeStrength) {
                    nPreChange[i]++;
                    nPostChange[j]++;
                    nTotChanges++;
                    W[i][j] = 0;
                    C1[i][j] = 0;
                    dWstruct[i][j] = 0;
                    int newi = i;
                    int newj = j;
                    int preChange = 0;
                    int postChange = 0;
                    if (preShift && !postShift) {
                        preChange = 1;
                    }
                    if (!preShift && postShift) { 
                        postChange = 1;
                    }
                    if (preShift && postShift) {
                        if (rand1() < 0.5) { 
                            preChange = 1;
						} else { 
							postChange = 1;
						}
					}
                    if (preChange) {
                        bool oldconn = 1;
                        while (oldconn) {
                            newi = int(rand1() * NE);
                            oldconn = C1[newi][j];
                        }
                        W[newi][j] = W0 * (1 + ranWStrengths * (rand1() - 0.5));
                        C1[newi][j] = 1;
                        dWstruct[newi][j] = 0;
                    }
                    if (postChange) {
                        bool oldconn = 1;
                        while (oldconn) {
                            newj = int(rand1() * NE);
                            oldconn = C1[i][newj];
                        }
                        W[i][newj] = W0 * (1 + ranWStrengths * (rand1() - 0.5));
                        C1[i][newj] = 1;
                        dWstruct[i][newj] = 0;
                    }
                }  
                dWstruct[i][j] *= structMemoryFact;
                if (dWstruct[i][j] > structMax) {
                    dWstruct[i][j] = structMax;
                }
            } //end if i->j connection
        } // next j
    } // next i
	cout << "Mean dWstruct " << mean_dWStruct / n_dW_Struct << " N " << n_dW_Struct << " max " << max_dWStruct << " min " << min_dWStruct << endl; 
} 

// Here, we fill excitatory bin rate.
void fill_EBinRate(int dNPools, int neurons, int dNEPool, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes, double bindT, int nBins, vector< vector<double> >& eBinRate) {
	for (int pool= 0; pool < dNPools; pool++) {
		for (int cell = neurons + pool * dNEPool; cell < neurons + (pool+1) * dNEPool; cell++) {
			for (int ispike = 0; ispike < numberSpikes[cell]; ispike++) {
				int bin = int(spikeTimes[cell][ispike] / bindT);
				eBinRate[pool][bin]++;
			}
		}
		for (int bin = 0; bin<nBins; bin++) {
			eBinRate[pool][bin] /= (dNEPool * bindT);
		}
    }
    eBinRate[0][0] = 0; eBinRate[1][0] = 0;  //Initial bin is zero.
}

// Here, we print out excitatory bin rate.
void print_Excitatory_BinRate(int nBins, double bindT, int dNPools, vector< vector<double> >& eBinRate) {
    ofstream poolout1;
    ofstream poolout2;
    poolout1.open("pool1.dat");
    poolout2.open("pool2.dat"); 
    for (int bin = 1; bin < nBins - 1; bin++) {
        cout << "Time bin: " << bindT * bin << "   ";
        for (int pool = 0; pool < dNPools; pool++) {
            cout << eBinRate[pool][bin] << "   ";
        }
        cout << endl;
        poolout1 << bindT * bin << " " << eBinRate[0][bin] << endl;
        poolout2 << bindT * bin << " " << eBinRate[1][bin] << endl;
    }
    poolout1.close();
    poolout2.close();
}

// Here, we fill out dPool_Spikes and calculate Pool_ssFr.
void fill_poolssFR_And_Dpoolspikes(vector<double>& pool_ssFR, vector< int>& dPool_Spikes, double tOn1, double tOff1, int dNPools, vector< vector<double> >& eBinRate) {
    for (int bin = tOn1 * 10; bin < tOff1 * 10; bin++) { // last second
        for (int pool = 0; pool < dNPools; pool++) {
            dPool_Spikes[pool] += eBinRate[pool][bin];
        }
    }
    for (int bin = (tOn1 * 10); bin < tOff1 * 10; bin++) { // last second
        if (dPool_Spikes[0] - dPool_Spikes[1] > 0) {
            pool_ssFR[0] += 1;
        }
        if (dPool_Spikes[1] - dPool_Spikes[0] > 0) {
            pool_ssFR[1] += 1;
        }
    }
    cout << "pool_ssFR " << pool_ssFR[0] <<  "  " << pool_ssFR[1] << endl;
}

// Here, we find our mean cue rate.
void mean_Cue_Rate(int NE, double& mean_Network_ENeuron_Rate, vector< int>& cueFR_2ndStim_Perst) {
    for (int j = 0; j < NE; j++) {
        // mean_Network_ENeuron_Rate += cueFR[j] * cueFR[j];
        mean_Network_ENeuron_Rate += cueFR_2ndStim_Perst[j] * cueFR_2ndStim_Perst[j];
    }
    mean_Network_ENeuron_Rate /= NE;
    cout << "mean squared network rate " << mean_Network_ENeuron_Rate << endl;
}

// Here, we find out what our population is (for a 2 second trial).
void find_Population(double eBinRate0testBin, double eBinRate1testBin, int& population) {
    if (eBinRate0testBin > eBinRate1testBin + 20.0) {
        population = 1;
    }
    if (eBinRate1testBin > eBinRate0testBin + 20.0) {
        population = 2;
    }
}

// Here, we find our populations' bin reaction time.
void find_Reaction_Time(int nBins, vector< vector<double> >& eBinRate) {
    ofstream RTout;
    RTout.open("RT.dat", ios::app);
    for (int bin = 0; bin < nBins; bin++) { 
        if (eBinRate[0][bin] > eBinRate[1][bin] + 20.0 && eBinRate[0][bin] >= 10) {
            cout << "Population 1 bin Reaction Time: " << bin << endl;
            RTout << bin << endl;
            break;
        }
        if (eBinRate[1][bin] > eBinRate[0][bin] + 20.0 && eBinRate[1][bin] >= 10) {
            cout << "Population 2 bin Reaction Time: " << bin << endl;
            RTout << bin << endl;
            break;
        }
    }
    RTout.close();
}

// Here, we calculate our reward.
void find_Reward(int stim, int population, int& reward) {
    if (((stim == 1 || stim == 4) && population == 1) || ((stim == 2 || stim == 3) && population == 2)) {
        reward = 1;
    } else {
        reward = -1;
    }
}

// Here, we add reward prediction error factor into dopamine signal.
void add_Reward_Dopamine(int trialBatch, int nHist, double& expectReward, vector<double>& histFactor, vector< vector< int> >& decisionHist, int stim) {
    if (trialBatch >= nHist) {
        for (int i = 1; i <= nHist; i++) {
            expectReward += histFactor[i - 1] * decisionHist[trialBatch - i][stim - 1];
        }
    } else {
        for (int i = 1; i <= trialBatch; i++) {
            expectReward += histFactor[i - 1] * decisionHist[trialBatch - i][stim - 1];
        }
    }
    if (expectReward < 0) {
        expectReward = 0; // No punishment
    }
 }

// Here, we display our reward expectation error.
void display_Reward_Expectation_Error(double& rewardExpectationError, int reward, double expectReward, int stim) {
    rewardExpectationError = reward - expectReward;
    cout << "rewardExpectationError " << rewardExpectationError << " reward " << reward << " stim " << stim << endl;
}
 
// Here, we implement our Poisson DA reward rule.
void poisson_DA_Reward_Rule(int NE, vector< int>& cueFR_2ndStim_Perst, double mean_Network_ENeuron_Rate, double DA_Epsilon, double rewardExpectationError, int population, int neurons, int dNEPool, vector< vector<double> >& W, double dLayerInput_W0, int dNE) {
    for (int i = 0; i < NE; i++) {
        double deltaWDA = ((cueFR_2ndStim_Perst[i]*cueFR_2ndStim_Perst[i]) - mean_Network_ENeuron_Rate) * DA_Epsilon;
        if (deltaWDA > 0.5) {
            deltaWDA = 0.5; // max 50% change per trial to synapse
        }
        if (deltaWDA < -0.5) {
            deltaWDA = -0.5;
        }
        deltaWDA *= rewardExpectationError;
        if (population == 1) {
            for (int j = neurons; j < neurons + dNEPool; j++) {
                W[i][j] += deltaWDA * dLayerInput_W0;
            }
        }
        if (population == 2) {
            for (int j = neurons + dNEPool; j < neurons + dNE; j++) {
                W[i][j] += deltaWDA * dLayerInput_W0;
            }
        }
	}
}
 
// Here, we compute our max and min change in synapse.
void compute_Max_Min_Change_Synapse(int NE, int neurons, int dNE, vector< vector<double> >& W, double dLayerInput_W0, double&  mean_W_AD1, double& mean_W_AD2) {
    for (int i = 0; i < NE; i++) {
        for (int j = neurons; j < neurons + dNE; j++) {
            if (W[i][j] < 0) {
                W[i][j] = 0;
            }
            if (W[i][j] > 2.5 * dLayerInput_W0) {
                W[i][j] = 2.5 * dLayerInput_W0;
            }
            if (j < neurons + dNE / 2) { 
                mean_W_AD1 += W[i][j];
            } else {
                mean_W_AD2 += W[i][j];
            }
        }
	}
	cout << "Mean W_AD 1 then 2 " << mean_W_AD1 / double(NE * dNE / 2) << " " << mean_W_AD2 / double(NE * dNE / 2) << endl;
}

// Here, we print a selection of final information at the end of our plasticity if statement.
void print_Final_Plasticity_Info(double eBinRate0testBin, double eBinRate1testBin, int population, int reward, int stim) {
    cout << "population 1 last bin: " << eBinRate0testBin << endl;
	cout << "population 2 last bin: " << eBinRate1testBin << endl;
	if (population == 1) {
        cout << "Release " << endl;
    }
	if (population == 2) {
        cout << "Hold " << endl;
    }
	// Output reward
	if (reward == 1) {
        cout << "Correct response to stim: " << stim << endl;
	} else {
        cout << "Incorrect response to stim: " << stim << endl;
    }
	cout << "reward = " << reward << endl;
}

// Here, we print our WR_ds values.
void print_WR_ds(int NE, vector<double>& WR_ds) { 
    ofstream WDSrout;
    WDSrout.open("WDSr.dat", ios::app);
    for (int j = 0; j < NE; j++) {
        WDSrout << WR_ds[j] << " " << j << endl;
    }
    WDSrout.close();
}

// Here, we print our Wds values.
void print_Wds(int trial, int NE, vector<double>& Wds) {
    if (trial == 1 || trial % 4 == 1 || trial == nTrials) {
        ofstream Wdsout3;
        Wdsout3.open("Wds.dat", ios::app);
        for (int j = 0; j < NE; j++) {
            Wdsout3 << Wds[j] <<  " " << j << endl;
        }
        Wdsout3.close();
    }
}

// Here, we record which trials are rewarded.
void record_Trials_Rewarded(int reward) {
    ofstream Rewardout;
    Rewardout.open("Reward.dat", ios::app);
    Rewardout << reward << endl;
    Rewardout.close();
}

// Here, we print out our excitatory bin rate.
void eOut_Trial(int nBins, double bindT, int dNPools, vector< vector<double> >& eBinRate) {
    ofstream Ebintrialout;
    Ebintrialout.open("Ebintrial.dat", ios::app);
    for (int bin= 0; bin < nBins; bin++) {
        Ebintrialout << bindT * bin << " ";
        for (int pool = 0; pool < dNPools; pool++) {
            Ebintrialout << eBinRate[pool][bin] << " ";
        }
        Ebintrialout << endl;
    }
    Ebintrialout.close();
}

// Here, we print our stimuli order.
void print_Stimuli_Order(int inputs, vector< int>& stimShuffle) {
    ofstream Stimuliout;
    Stimuliout.open("Stimuli.dat", ios::app);
    for (int i = 0; i < inputs; i++) {
        Stimuliout << stimShuffle[i] << endl;
    }
    Stimuliout.close();
}

// The createString method will be used to create the trial headers used 
// throughout the program.  The method's parameters include a string (which will 
// be our name prefix) and two integers: num (which will be the number inserted 
// into the header) and digits (which will be how many digits the header will 
// contain).
string createString(string prefix, int num, int digits) {
	// First, we insert our number into the display string.
	string display;
    stringstream sstm;
    sstm << num;
    display = sstm.str();
	// Then, while we haven't yet reached the number of digits we want, we add 0s, one at a time.
	while (display.length() < digits) {
		display = "0" + display;
	}
	// Finally, we add our name prefix to the beginning of the string, print it, and return it.
	display = prefix + display;
	return display;
}

// Here, we output our spike raster.
void output_Spike_Raster(int trial, int neurons, int dNeurons, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes) {
    ofstream Rasterout;
    string Rasterstring = createString("Raster_", trial, 3);
    Rasterstring += ".dat";
    Rasterout.open(Rasterstring.c_str());
    // Raster output
    for (int i = 0; i < neurons + dNeurons; i++) {
        for (int j = 0; j < numberSpikes[i]; j++) {
            Rasterout << spikeTimes[i][j] << " " << i << endl;
        }
    }
    Rasterout.close();
}

// Here, we output our current trial raster.
void output_Current_Trial_Raster(int neurons, int dNeurons, vector< int>& numberSpikes, vector< vector<double> >& spikeTimes) {
    ofstream currentRaster;  // outputs current trial raster, use for PC simulations to monitor persistence, not cluster.
    currentRaster.open("currentRaster.dat");
    for (int i = 0; i < neurons + dNeurons; i++) {
        for (int j = 0; j < numberSpikes[i]; j++) {
            currentRaster << spikeTimes[i][j] << " " << i << endl;
        }
    }
    currentRaster.close();
}

// Here, we output our mean excitatory rate.
void write_CueFr_1(int NE, vector< int>& cueFR) {
    ofstream cueFRout2;
    cueFRout2.open("cueFR.dat", ios::app);
    double meanErate = 0;
    for (int j = 0; j < NE; j++) {
        cueFRout2 << cueFR[j] << " " << j << endl;
        meanErate += cueFR[j];
    }
    cueFRout2.close();
    cout << "mean E rate " << meanErate / double(NE);
}

// Here, we output our mean I rate.
void write_CueFr_2(int NE, int neurons, vector< int>& cueFR) {
    ofstream cueFRoutI;
    cueFRoutI.open("cueFRI.dat", ios::app);
    double meanIrate = 0;
    for (int j = NE; j < neurons; j++) {
        cueFRoutI << cueFR[j] << " " << j << endl;
        meanIrate += cueFR[j];
    }
    cueFRoutI.close();
    cout << " mean I rate " << meanIrate / double(NI);
}

// Here, we output our cueFR_1stStim values.
void write_CueFr_1stStim(int NE, vector< int>& cueFR_1stStim) {
    ofstream cueFRout5;
    cueFRout5.open("cueFR_1stStim.dat", ios::app);
    for (int j = 0; j < NE; j++) {
        cueFRout5 << cueFR_1stStim[j] << " " << j << endl;
    }
    cueFRout5.close();
}

// Here, we output our cueFR perst values.
void write_CueFr_Perst(int NE, vector< int>& cueFR_Perst) {
    ofstream cueFRout4;
    cueFRout4.open("cueFR_Perst.dat", ios::app);
    for (int j = 0; j < NE; j++) {
        cueFRout4 << cueFR_Perst[j] << " " << j << endl;
    }
    cueFRout4.close();
}

// Here, we output our cueFR 2nd perst values.
void write_CueFr_2nd_Perst(int NE, vector< int>& cueFR_2nd_Perst) {
    ofstream cueFRout6;
    cueFRout6.open("cueFR_2nd_Perst.dat", ios::app);    
    for (int j = 0; j < NE; j++) {
        cueFRout6 << cueFR_2nd_Perst[j] << " " << j << endl;
    }
    cueFRout6.close();
}

// Here, we output our cueFR 2nd stim perst values.
void write_CueFr_2ndStim_Perst(int NE, vector< int>& cueFR_2ndStim_Perst) {
    ofstream cueFRout3;
    cueFRout3.open("cueFR_2ndStim_Perst.dat", ios::app);
    for (int j = 0; j < NE; j++) {
        cueFRout3 << cueFR_2ndStim_Perst[j] << " " << j << endl;
    }
    cueFRout3.close();
}

// Here, we write our mean excitatory neuron delay rate and our mean excitatory neuron 2nd delay rate.
void write_Mean_PerstFR(double mean_ENeuron_Delay_Rate, double mean_ENeuron_2nd_Delay_Rate) {
     // First
     ofstream mean_perstFR;
     mean_perstFR.open("mean_perstFR.dat", ios::app);
     mean_perstFR << mean_ENeuron_Delay_Rate << endl;
     mean_perstFR.close();
     // Second
     ofstream mean_2nd_perstFR;
     mean_2nd_perstFR.open("mean_2nd_perstFR.dat", ios::app);
     mean_2nd_perstFR << mean_ENeuron_2nd_Delay_Rate << endl;
     mean_2nd_perstFR.close();
}

// Here, we write out our mean AL excitatory cell bin rate for a trial.
void AL_ECell_BinTrial(int nBins, double bindT, int neurons, vector< vector<double> >& AL_ECell_BinRate) {
    ofstream AL_ECell_BinTrialout;
    AL_ECell_BinTrialout.open("Mean_AL_ECell_BinRate.dat", ios::app);
    for (int bin = 0; bin < nBins; bin++) {
        AL_ECell_BinTrialout << bindT * bin << " ";
        for (int cell = 0; cell < neurons; cell++) {
            AL_ECell_BinTrialout << AL_ECell_BinRate[cell][bin] << " ";
        }
        AL_ECell_BinTrialout << endl;
    }
    AL_ECell_BinTrialout.close();
}

// Here, we write out our mean AL Excitatory cell bin rate for a trial.
void output_AL_ECell_BinRate(int trial, int nBins, double bindT, int neurons, vector< vector<double> >& AL_ECell_BinRate) {
    ofstream AL_ECell_BinRate_out;
    string AL_ECell_BinRate_string = createString("Mean_AL_ECell_BinRate_out_", trial, 3);
    AL_ECell_BinRate_string += ".dat";
    AL_ECell_BinRate_out.open(AL_ECell_BinRate_string.c_str());
    for (int bin= 0; bin < nBins; bin++) {
        AL_ECell_BinRate_out << bindT * bin << " ";
        for (int cell = 0; cell < neurons; cell++) {
            AL_ECell_BinRate_out <<  AL_ECell_BinRate[cell][bin] << " ";
        }
        AL_ECell_BinRate_out << endl;
    }
    AL_ECell_BinRate_out.close();
}

// Here, we output our final W and pW values.
void print_Final_W_And_pW(int neurons, int dNeurons, vector< vector<double> >& W, vector< vector<double> >& pW, int nIn) {
    ofstream W_out;
    W_out.open("W.dat");
    for (int i = 0; i < neurons + dNeurons; i++) {
        for (int j = 0; j < neurons + dNeurons; j++) {
            W_out << W[i][j] << " "; //<< " i " << i << " j " << j << endl;
        }
        W_out << endl;
    }
    W_out.close();
    ofstream pW_out;  // output initial input weight matrix
    pW_out.open("pW.dat");
    for (int i = 0; i < nIn; i++) {
        for (int j = 0; j < neurons + dNeurons; j++) {
            pW_out << pW[i][j] << " ";
        }
        pW_out << endl;
    }
    pW_out.close();
 }

// Here, we output our W values each 10 trials (as well as final trial).
void output_Regularly_W(int trial, int neurons, int dNeurons, vector< vector<double> >& W) {
    ofstream Wout;
    string Wstring = createString("Wout_", trial, 3);
    Wstring += ".dat";
    Wout.open(Wstring.c_str());  
    for (int i = 0; i < neurons + dNeurons; i++) {
        for (int j = 0; j < neurons + dNeurons; j++) {
            Wout << W[i][j] << " ";
        }
        Wout << endl;
    }
    Wout.close();
}

// Here, we output our pW values each 10 trials (as well as final trial).
void output_Regularly_pW(int trial, int nIn, int neurons, int dNeurons, vector< vector<double> >& pW) {
    ofstream pWout2;
    string pWstring = createString("pWout_", trial, 3);
    pWstring += ".dat";
    pWout2.open(pWstring.c_str());  
    for (int i = 0; i < nIn; i++) {
        for (int j = 0; j < neurons + dNeurons; j++) {
            pWout2 << W[i][j] << " ";
        }
        pWout2 << endl;
    }
    pWout2.close();
}

// Here, we output our voltage data.
void output_Voltage_Data(int nT, int neurons, int dNeurons, vector< vector<double> >& V) {   
    ofstream V_out;
    V_out.open("V.dat");	
    for (int i = 0; i <= nT; i++) {
        for (int j = neurons; j < neurons + dNeurons; j++) {
            V_out << V[j][i] << " ";
        }
        V_out << endl;
    }
    V_out.close();
}

// Calculate and print meanWEE.
void print_Mean_wEE(int NE, vector< vector<double> >& W) {
    double sumwEE = 0;
    int count_wEE = 0;
    for (int i = 0; i < NE; i++) {
        for (int j = 0; j < NE; j++) {
            if (W[i][j] > 0) {
                sumwEE += W[i][j];
                count_wEE++;
            }
        }
    }
    cout << " meanwEE " << sumwEE / double(count_wEE);
}

// Calculate and print meanWIE.
void print_Mean_wIE(int NE, int neurons, vector< vector<double> >& W) {
    double sumwIE = 0;
    int count_wIE = 0;
    for (int i = NE; i < neurons; i++) {
        for (int j = 0; j < NE; j++) {
            if (W[i][j] > 0) {
                sumwIE += W[i][j];
                count_wIE++;
            }
        }
    }
    cout << " meanwIE " << sumwIE / double(count_wIE);
}

// Calculate and print meanWInE.
void print_Mean_wInE(int nIn, int NE, vector< vector<double> >& pW) {
    double sumWInE = 0;
    int countWInE = 0;
    for (int i = 0; i < nIn; i++) {
        for (int j = 0; j < NE; j++) {
            if (pW[i][j] > 0) {
                sumWInE += pW[i][j];
                countWInE++;
            }
        }
    }
    cout << " meanWInE " << sumWInE / double(countWInE);
}

// Calculate and print meanWInI.
void print_Mean_wInI(int nIn, int neurons, vector< vector<double> >& pW) {
    double sumWInI = 0;
    int countWInI = 0;
    for (int i = 0; i < nIn; i++) {
        for (int j = NE; j < neurons; j++) {
            if (pW[i][j] > 0) {
                sumWInI += pW[i][j];
                countWInI++;
            }
        }
    }
    cout << " meanWInI " << sumWInI / double(countWInI) << endl;
}

// Here, we output our nTotalchange value as well our nPreChange and nPostChange values.
void output_N_Changes(int nTotChanges, int NE, vector< int>& nPreChange, vector< int>& nPostChange) { 
    cout << " nTotalchange " << nTotChanges << endl;
    ofstream NchangesOut;
    NchangesOut.open("Nchanges.dat");
    for (int i = 0; i < NE; i++) {
        NchangesOut << i << " " << nPreChange[i] << " " << nPostChange [i] << endl;
    }
}

/*************************************************************************************/
/* Step 3: Initializations Before Trials                                             */
/*************************************************************************************/

// Note that the main function accepts a parameter from the user; this will be
// used as the value numIn below.
int main(int argc, char **argv) {
    ////////// STEP 3A: NUMIN, SEEDS, AND DIGITS //////////
    // To begin, let's print the starting time for our program.
    time_t rawTime;
    time(&rawTime);
    printf("The current local time is: %s", ctime (&rawTime));
    // We declare numIn below.  This will be used to generate our seeds, as well as
    // to find the indices we'll use in our input cells and input probability arrays.
    int numIn = 0;
    // If the user gives a command line argument, that argument will be inserted into
    // numIn.
    if(argc > 1) {
        numIn = atoi(argv[1]);
    }
    // Now, let's find our seed values, using the findSeeds function defined above.
    findSeeds(seed1, seed2, numIn);
    // If numIn is greater than 25, we return an error because that means that our digits won't work
    // with our InputCells and InputProbability matrices.
    // Also, numIn cannot be 0, so we should also give an error if this is the case.
    if (numIn > 25 || numIn == 0) {
        cout << "Error: your numIn is 0 or is greater than 25 after seed reductions." << endl;
        exit(1);
    }
    // Assuming that our numIn is less than or equal to 25 (and greater than 0), we will find 
    // the digits for the indices of our input cells (digit1) and input probability (digit2) arrays.  
    int digit1 = 0;
    int digit2 = 0;
    // We will use the findDigits function defined above.
    findDigits(digit1, digit2, numIn);
    // Now, we'll use these digits to find the input for our cells and probability.
	int inputCells = 0;
    double inputProb = 0;
    // We use the findinputs function defined above.
    findinputs(digit1, digit2, inputCells, inputProb, inputCellArray, inputProbArray);
    // We find our nIn value by multiplying the number of inputs by the number of input cells
    const int nIn = inputs * inputCells;
    // We'll also compute our r0 value, which is the input rate (Hz) per input cell.
    double r0 = 480.0 / inputCells / (3 * inputProb);
    // Now, having computed our seed values, we set the seeds for our 
    // MersenneTwister random objects created above.
    rand1.seed(seed1);
    rand2.seed(seed2);
    rand3.seed(seed3);
    // We will use srand when using the random_shuffle function below
    // to select our stim value for each trial.  
    srand(seed1);       
    ////////// STEP 3B: REWARD PREDICTION, TIME, AND GOAL RATE //////////
    // REWARD PREDICTION HISTORY
    fillHistFactor(histFactor, sumHistFactor, histDecay, nHist);
    // TIME VECTORS (t, preliminary values for tOn and tOff)
    fillTime(t, tOn, tOff, cue1Length, cue2Length, delay1Length, delay2Length, tMax, nT, dt);
    // GOAL RATE VECTORS, WITH NOISE FOR HETEROGENEITY
    fillGoalVectors(r_Goal, r_Goal_I_to_E, r_DGoal, NE, NI, dNE, rand1, rgE, rgIE0, rgI, rgD);
    ////////// STEP 2C: WEIGHT-RELATED INITIALIZATIONS //////////
    // Initialize weight matrices.  We have to initialize these in the main method because
    // their size depends on nIn, a number computed from our numIn value.
    vector< vector<double> > eW(nIn, DDNeuronVec);
    vector< vector<double> > pW(nIn, DDNeuronVec);
    bool C3[nIn][neurons + dNeurons];  // pW after P(removal input cell)
    // Now, let's fill our loops (NE-NE, NI-NE, NI-NI, NE-NI (recurrent I))
    fillLoops(NE, neurons, connProb, wEE_xFactor, W0, ranWStrengths, wIE_ConnProb, wEI_ConnProb, wIE0, wII0, wEI0, dNeurons, rand1, W, dWstruct);
    // We'll also randomize our input weight matrix eW by input cell number
    randomize_eW(nIn, NE, NI, rand1, inputProb, eW);
    // Let's also introduce the heterogeneity of pW
    heterogeneity_pW(nIn, NE, pW0, pW0I, ranWStrengths, rand1, neurons, dNeurons, pW, eW, feedForwardI);
    // If we want to do the Poisson test, we implement that test in our pW matrix.
    if (DL_Poisson_Test) {
        poisson_pW(nIn, neurons, dNE, dNEPool, pool1Weight, pool2Weight, pW);
    }
    // Weight modification with pools
    pools_W(NE, neurons, dNeurons, dNE, dLayerInput_W0, dNPools, dNEPool, wEE, dW0, dNIPool, wEI, wIE, wII, W);
    // Modify input weights
    modify_pW(nIn, NE, NI, pW02, pW0I2, pW);
    // Network connection bool and normalize matrices
    normalize_And_C1C3Bools(neurons, dNeurons, nIn, W1, W, pW, C1, C3);
    // Initialize heterogeneity vectors
    initialize_Heterogeneity_Vectors(neurons, NE, dNE, dNeurons, rand1, E, E_s, tau_M, tau_M_S, g_L, AL_g_L_s, Vth_E, vth_E_sA, VReset_E, VReset_E_s, TRef_E, TRef_E_s, vth_I_sA, vReset_I_s, tRef_I_s, g_L_s, vth_E_sD, vth_I_sD);
    // Initialize stim shuffle vector
    initialize_StimShuffle(inputs, stimShuffle);
    
    /*************************************************************************************/
    /* Step 4: Printing Important Information for User                                   */
    /*************************************************************************************/ 
    
    // We want the user to know if he/she is using facilitation and/or
    // depressing synapses, so we'll report that information.
    if (facilitation) {
        cout << "Using Facilitating synapses" << endl;
    }
    if (depression) {
        cout << "Using Depressing synapses" << endl;
    }
    if (structPlast) {
        cout << "Structural Plasticity On" << endl;
    }
    // Let's print our sigma values.
    cout << "Sigma = " << sigma << " " << sigma_E << " " << sigma_I << endl;
    // Let's print the number of neurons in our network and the number of 
    // trials we'll be performing.
    cout << "Number of trials: " << nTrials << endl;
    cout << "Neurons: " << neurons << endl;
    // Let's print some decision layer values.
    cout << "Decision neurons: " << dNeurons << " dNEPool: " << dNEPool << " dNIPool: " << dNIPool << endl;
    cout << "r0: " << r0 << endl;
    cout << "cue1Length = " << cue1Length << " cue2Length = " << cue2Length << " 1st delay1Length = " << delay1Length << " 2nd delay1Length = " << delay2Length << endl;
    cout << "rgIE0 = " << rgIE0 + rgE << endl;
    cout << "pW0 " <<  pW0  << endl;
    if (DL_Poisson_Test) {
        cout << "Using Poisson inputs for testing" << endl; 
        cout << "Stimulus 1 bias = " << stimBiasPerc << endl;
    }
    cout << "Associative Input to Decision layer weight " << dLayerInput_W0 << endl;
    cout << "pop1 wEI" << endl;
    cout << "pop2 wEI " << endl;
    cout << "dLayerInput = " << dLayerInput << endl;
    cout << "wEE = " << wEE << endl;
    cout << "wEI = " << wEI << endl;
    cout << "wIE = " << wIE << endl;
    cout << "wII = " << wII << endl;
    cout << "Input weights modified: " << pW02 << " " << pW0I2 << endl; 
    // Print out W and pW to files
    print_W(neurons, dNeurons, W);
    print_pW(nIn, NE, NI, dNeurons, pW);
    
    /*************************************************************************************/
    /* Step 5: Trial Loop                                                                */
    /*************************************************************************************/
    
    for (int trial = 1; trial <= nTrials; trial++) { 
        cout << "Trial " << trial << endl; // print the current trial out
        // Initialize vector/matrices to zero for each individual trial. 
        trial_Values_to_0(spikeTimes, numberSpikes, gSyn_I_Vec, gSyn_E_Vec, gpSyn_Vec, g_Tot, vInf, V, lastSpike, syn, syn_I, syn_E, tau_Eff, pSynMax, synMax_I, synMax_E, synMax_NMDA, someIndex, gSyn_NMDA, g_NMDA, Vth, syn_NMDA, g_Ref, meanFiringRate, W_Shift, W_Shift_I, DW_Shift, Fac, synE_Fac, synNMDA_Fac, Dep, synE_Dep, synNMDA_Dep, neurons, dNeurons, max_CellSpikes, nT);
        vector< vector<double> > pSpikes(nIn, vec_Max_InputSpikes);
        vector<double> pSyn(nIn);
        // Initialize lastSpike vector
        initialize_LastSpike(neurons, dNeurons, NE, dNE, TRef_E, tRef_I, lastSpike);
        // Initializing initial voltage
        initialize_InitialVoltage(neurons, dNeurons, V, E);
        // Initialize facilitation / depression factors
        initialize_FacAndDep(facilitation, depression, NE, Fac, Dep, F0, D0);        
        // Using a fixed stimulus below
        stimSetup(stim, trial, stimShuffle);
        // Poisson inputs
        vector<int> nPSpikes(nIn);
        ////// Time Integration //////
        for (int i = 1; i <= nT; i++) {
            // Define pSyn during each input, A, B,// AB
            pSyn_AMPADecay(nIn, pSyn, AMPADecay);        	
            bool InputOn = 0;
            bool InputOn2 = 0; 
            int cellStart = 0;
            int cellStop = 0;
            int cellStart2 = 0;
            int cellStop2 = 0;
            // conductance that ramps up for DL cells when urgency == 1 
            double g_Urgency = 0;
            // Initialize cellStart / cellStop / cellStart2 / cellStop2
            initialize_CellStartAndStop(t[i], tOn[0], tOff[0], tOn[1], tOff[1], InputOn, InputOn2, stim, cellStart, cellStop, cellStart2, cellStop2, inputCells, urgency, g_Urgency, gMax_Urgency);
            // First stimulus of pair, then second one
            if (InputOn) {
                use_InputOn(cellStart, cellStop, dt, rand3, r0, nPSpikes, max_InputSpikes, t[i], pSpikes, pSynMax, s0, pSyn);
            }
            if (InputOn2) {
                use_InputOn(cellStart2, cellStop2, dt, rand3, r0, nPSpikes, max_InputSpikes, t[i], pSpikes, pSynMax, s0, pSyn);
            }
            ////// Network Code //////
            for (int j = 0; j < neurons + dNeurons; j++) {  // probably need to run over both layers 
                if (numberSpikes[j] > max_CellSpikes) {
                    cout << "Pre G " << j <<  " " << numberSpikes[j] <<  " " << t[i] << endl;
                }
                gSyn_I_Vec[j] = 0;
                gSyn_E_Vec[j] = 0;
                gSyn_NMDA[j] = 0;
                // Icells network -> NE-neurons or neurons + dNE:neurons + dNeurons
                // Ecells network -> 0-NE or neurons-dNE
                synapse_Cond_Update(neurons, dNeurons, j, C1, NE, dNE, gSyn_I_Vec, W, syn_I, syn_E, facilitation, depression, gSyn_E_Vec, gSyn_NMDA, syn_NMDA, synE_Fac, synNMDA_Fac, synE_Dep, synNMDA_Dep);
                // Separate conductances for Associate (AL) /Decision (DL) Layer
                ALDL_Cond_Update(j, gSyn_I_Vec, gSyn_I, NE, gSyn_E_Vec, GAL_E_AMPA, gSyn_NMDA, GAL_E_NMDA, neurons, GAL_I_AMPA, GAL_I_NMDA, dNE, GDL_AMPA, GDL_NMDA, nSigma_I, nSigma_E, rand2, dt, sigma_I, sigma_E, t[i], tOn[1]);
                // Update gpSyn_Vec
                update_gpSyn_Vec(j, nIn, gpSyn_Vec, pW, gpSyn, pSyn, urgency, neurons, dNE, g_Urgency);             
                // Refractory conductance decay for Inhibitory - Excitatory cells.
                update_Refract_Cond(j, NE, neurons, dNE, g_Ref, dt, tRef_I[j], TRef_E[j]);
                // Total input conductance
                g_Tot = g_L[j] + g_Ref[j] + gSyn_NMDA[j] * g_NMDA[j] + gSyn_E_Vec[j] + gSyn_I_Vec[j] + gpSyn_Vec[j];
                // Associative layer and Decision layer use different voltage noise terms {nSigma, and Sigma} respectively
                update_vInf(vInf, j, neurons, g_L[j], E[j], g_Ref[j], vReset, gSyn_NMDA[j], g_NMDA[j], E_NMDA, gSyn_E_Vec[j], E_AMPA, gSyn_I_Vec[j], E_GABA, gpSyn_Vec[j], nSigma, sigma, dt, g_Tot);           
                tau_Eff = (tau_M[j] * g_L[j]) / g_Tot; // Effective time constant changes with conductances (g_Tot)
                V[j][i] = vInf[j] + (V[j][i - 1] - vInf[j]) * exp(-dt / tau_Eff); // V(Neuron) 	
                g_NMDA[j] = 1 / (1 + mg_Ext * exp(-0.062 * (-V[j][i]) / 3.57e-3));  // Compte...Wang 02 g_NMDA
                // Dynamic Thresholds, resets, and refractory conductances
                update_Dynamic_Thresholds(j, NE, neurons, dNE, vth0, vth_I[j], Vth_E[j], tRef, tRef_I[j], TRef_E[j], Vth, vth_Max, t[i], lastSpike[j]);
                // Resets and refractory conductances
                inhib_Excit_Refract_Reset(j, i, NE, neurons, dNE, V, vSpike, vReset_I[j], VReset_E[j], Vth[j], g_Ref, delta_GRef, lastSpike, t[i]);
                if (numberSpikes[j] > max_CellSpikes) { 
                    cout << "Done Thresh " << j <<  " " << numberSpikes[j] << " " << t[i] << endl;
                }
                // Record spike times
                record_spikeTimes(j, V[j][i], vSpike, numberSpikes, someIndex, spikeTimes, t[i]);
                // Update synaptic conductances: GABA, AMPA, NMDA
                update_Syn_Cond(j, syn_I, syn_E, GABADecay, AMPADecay, neurons, syn_NMDA, AL_NMDADecay, DL_NMDADecay, synE_Fac, synNMDA_Fac, synE_Dep, synNMDA_Dep);
                // facilitation and depression decays
                update_Fac_Dep_Decays(facilitation, depression, j, NE, Fac, Dep, dt, tau_Fac, tau_Dep);
                // Update inhibitory and excitatory synapses
                update_Inhib_Excit_Synapses(j, V[j][i], vSpike, synMax_I, synMax_E, synMax_NMDA, s0, syn_I, syn_E, syn_NMDA, NE, neurons, dNE, facilitation, depression, synE_Fac, synNMDA_Fac, Fac, synE_Dep, synNMDA_Dep, Dep, alpha, fMax, dFrac);
            }  // end network code ****
        }  // end time integration
        vector<int> cueFR(neurons + dNeurons);
        // Write the cue firing rate of each neuron during stimulus presentation.
        write_Cue_Firing_Rate(neurons, numberSpikes, spikeTimes, tOn[0], tOff[0], tOn[1], tOff[1], cueFR, cue1Length, cue2Length);
        // Activity during the first stimulus to measure selectivity  from 1st stimulus.
        vector<int> cueFR_1stStim(neurons);
        measure_Selectivity_1st_Stim(neurons, numberSpikes, spikeTimes, tOn[0], tOff[0], tOn[1], tOff[1], cueFR_1stStim, cue1Length);
        // Activity during the second stimulus to measure selectivity with persistent activity from 1st stimulus. (When persistent selective activity has arisen).
        vector<int> cueFR_2ndStim_Perst(neurons);
        measure_Selectivity_2nd_StimPerst(neurons, numberSpikes, spikeTimes, tOn[1], tOff[1], cueFR_2ndStim_Perst, cue2Length);
        // Activity during the second stimulus to measure selectivity with persistent activity from 1st stimulus. (When persistent selective activity has arisen).
        vector<int> cueFR_Perst(neurons);
        measure_Selectivity_Perst(NE, numberSpikes, spikeTimes, tOff[0], tOn[1], cueFR_Perst, delay1Length);
        // Find mean excitatory persistent delay activity
        double mean_ENeuron_Delay_Rate = 0;
        find_Mean_ENeuron_Delay_Rate(mean_ENeuron_Delay_Rate, NE, cueFR_Perst);
        // Activity during the second delay period to measure selectivity with persistent activity from 1st+2nd stimulus. (When persistent selective activity has arisen).
        vector<int> cueFR_2nd_Perst(neurons);
        measure_Selectivity_2nd_Perst(NE, numberSpikes, spikeTimes, tOff[1], tMax, cueFR_2nd_Perst, delay2Length);
        // Finding mean excitatory 2nd persistent delay activity         
        double mean_ENeuron_2nd_Delay_Rate = 0;
        find_Mean_ENeuron_2nd_Delay_Rate(mean_ENeuron_2nd_Delay_Rate, NE, cueFR_2nd_Perst);
        // Find bins and binRate
        cout << "Number of bins = " << nBins << endl;
        vector<double> binSize(nBins);
        vector< vector<double> > AL_ECell_BinRate(neurons, binSize);
        find_Bin_Rate(neurons, nBins, numberSpikes, spikeTimes, bindT, AL_ECell_BinRate);       
        // Output the trial rate of AL neurons to DL weights to see if there is a shift
        vector<double> WR_ds(NE);
        vector<double> Wds(NE);
        double mean_WR_ds = 0;
        double mean_rE = 0;
        double mean_rI = 0;
        output_TrialRate_ALNeurons_DLWeights(NE, neurons, dNEPool, WR_ds, Wds, mean_WR_ds, W, cueFR);
        // Find mean firing rate and mean mean firing rate
        find_Mean_Mean_Firing_Rate(neurons, dNeurons, numberSpikes, tMax, meanFiringRate, mean_Mean_Firing_Rate, mean_rE, mean_rI, NE, NI);           
        // Excitatory mean firing rate averaged across Ecells and Icells
        double rgIE2 = (mean_rI * tMax * exp(-LTPi_Window * mean_rE) * idW) / e_I + mean_rE;  // currently uses existing mean rates
        // print out potential goal rates and mean network excitatory/inhibitory rates
        cout << "rgIE2: " << rgIE2 << " <rE> " << mean_rE << " <rI> " << mean_rI << endl; 
        if (plasticityOn) {  // Synaptic Plasticity Bool     
            if (oneTrialHomeo || trial % inputs == 0) {  // If oneTrialHomeo
                divide_Mean_Mean_Firing_Rate_By_Inputs(neurons, dNeurons, oneTrialHomeo, mean_Mean_Firing_Rate, inputs);     
                if (EE_Homeostasis || inputEHomeostasis) {
                    double sumWShift = 0;
                    input_And_Recurrent_Excitatory_Homeostatic(NE, W_Shift, e_E, r_Goal, mean_Mean_Firing_Rate, inputEHomeostasis, nIn, pW, EE_Homeostasis, W, sumWShift);
                }
                if (inputIHomeostasis) {
                    double sumWIShift = 0;
                    input_Inhibitory_Homeostatic(NE, NI, W_Shift, e_Input_I, r_Goal, mean_Mean_Firing_Rate, nIn, pW, sumWIShift); 
                }
                if (DA_Reward_Homeostasis) { // Decision layer input homeostasis
                    double meanDWShift = 0;
                    input_DL_Homeostatic(neurons, dNE, DW_Shift, eD, rgD, mean_Mean_Firing_Rate, NE, W, meanDWShift);
                }
            }
            // LTPi homeostasis
            if (LTPi) {               
                double sum_WShift_I = 0;
                LTPI_homeostasis(NE, W_Shift_I, e_I, r_Goal_I_to_E, meanFiringRate, NI, W, sum_WShift_I);
            }
            // Pfister & Gerster EESTDP rule
            if (PG_EESTDP) { 
                double sum_PGdw = 0;
                int countIJ = 0;
                vector< vector<double> > dwSum(NE, NEVec);				
                vector<int> counter(NE);
                double deltaT = 0;
                int endFlag = 0;
                PG_Rule_EESTDP(NE, C1, counter, numberSpikes, spikeTimes, tau3m, endFlag, deltaT, dwSum, tau2m, A2m, A3m, tau2p, tau3p, A2p, A3p, dW, W, W0, sum_PGdw, countIJ); 
            }
            // Pfister & Gerster Input to E cells STDP rule
            if (PG_Input_ESTDP) {         
                double sum_InputPGdw = 0;
                int countIJ = 0;
                vector < vector<double> > dwEinsum(nIn, NEVec);
                vector<int> counter(NE);
                int endFlag = 0;
                double deltaT = 0;
                PG_Rule_Input_ESTDP(NE, nIn, pW, counter, spikeTimes, numberSpikes, nPSpikes, pSpikes, tau3m, endFlag, deltaT, dwEinsum, tau2m, A2m, A3m, tau2p, tau3p, A2p, A3p, pW02, sum_InputPGdw, countIJ);
            } 
            // LTPi rule Maffei et al.
            if (LTPi) {
                vector<vector<double> > dwSum_I(NI, NEVec);  // initialize dwSum_I per trial.
                vector<vector<double> > LTPi_Shift(NI, NEVec);
                double sum_I_dwSum = 0;
                int countIJ = 0;
                LTPi_Rule(NE, neurons, C1, numberSpikes, spikeTimes, dt, VLTPi, V, VLTPi_Threshold, LTPi_Window, dwSum_I, idW, W, wIE0, sum_I_dwSum, countIJ);
            }
            //Floor and ceiling for pW and W
            floor_And_Ceiling(nIn, neurons, pW, W);
            // Section to remove and add new synapses
            if (structPlast) {
                double deltaT = 0;
                int endFlag = 0;
                double mean_dWStruct = 0;
                double max_dWStruct = -100;
                double min_dWStruct = 100;
                int n_dW_Struct = 0;
                vector<int> counter(NE);
                structPlast_Operation_1(NE, C1, counter, numberSpikes, spikeTimes, endFlag, tauStruct, deltaT, dWstruct, structRateFact, mean_dWStruct, n_dW_Struct, max_dWStruct, min_dWStruct);
                structPlast_Operation_2(NE, C1, dWstruct, removeStrength, nPreChange, nPostChange, nTotChanges, W, rand1, W0, ranWStrengths, structMemoryFact, structMax, mean_dWStruct, n_dW_Struct, max_dWStruct, min_dWStruct);

            } // end if structural plasticity is on
        }
        // Basic Hebbian DA reward learning rule
        if (DL_Out) { 		
            cout << "Number of bins = " << nBins << endl; // Edit this so it shows up as 20 instead of 19
            int testBin = int(tOff[1] / bindT) + 1;
            if (testBin > nBins - 1) {
                testBin = nBins - 1; // Time bin to test rate of "winning" pool
            }
            vector<double> binRate(nBins);
            vector< vector<double> > eBinRate(dNPools, binRate);
            fill_EBinRate(dNPools, neurons, dNEPool, numberSpikes, spikeTimes, bindT, nBins, eBinRate);
            // Print out excitatory bin rate
            print_Excitatory_BinRate(nBins, bindT, dNPools, eBinRate);
            // Fill pool_ssFR and dPool_Spikes
            vector<double> pool_ssFR(dNPools);
            vector<int> dPool_Spikes(dNPools);
            fill_poolssFR_And_Dpoolspikes(pool_ssFR, dPool_Spikes, tOn[1], tOff[1], dNPools, eBinRate);
            // Hebbian DA reward based plasticity rule
            double mean_Network_ENeuron_Rate = 0;
            int population = 0;
            int reward = 0;
            // Mean cue rate of neurons by each neurons total spikes, mean pop/N
            mean_Cue_Rate(NE, mean_Network_ENeuron_Rate, cueFR_2ndStim_Perst);
            // Find population for 2 second trial
            find_Population(eBinRate[0][testBin], eBinRate[1][testBin], population); 
            // Reaction Time - localized to bin (Reaction Time defined by first bin in which one population is firing at least 20Hz above the other, can also add threshold (e.g. also be above 20Hz or something)
            find_Reaction_Time(nBins, eBinRate);
            // if AB/CD are input and "Release" population wins => reward OR if AD/CB are input and "Hold" population wins => reward
            find_Reward(stim, population, reward);
            int trialBatch = (trial - 1) / inputs;
            decisionHist[trialBatch][stim - 1] = reward;
            if (plasticityOn) {	
                if (DA_Reward) { 
                    double expectReward = 0;
                    double rewardExpectationError = 0;
                    if (rewardExpectation) {
                        // ADD REWARD PREDICTION ERROR FACTOR INTO DOPAMINE SIGNAL
                        add_Reward_Dopamine(trialBatch, nHist, expectReward, histFactor, decisionHist, stim);
                    }
                    // Compute and display reward expectation error
                    display_Reward_Expectation_Error(rewardExpectationError, reward, expectReward, stim);              
                    // Basic Poisson DA reward rule.
                    poisson_DA_Reward_Rule(NE, cueFR_2ndStim_Perst, mean_Network_ENeuron_Rate, DA_Epsilon, rewardExpectationError, population, neurons, dNEPool, W, dLayerInput_W0, dNE);
                    // Max/min total change to synapse
                    double mean_W_AD1 = 0;
                    double mean_W_AD2 = 0;
                    compute_Max_Min_Change_Synapse(NE, neurons, dNE, W, dLayerInput_W0, mean_W_AD1, mean_W_AD2);     
                    print_Final_Plasticity_Info(eBinRate[0][testBin], eBinRate[1][testBin], population, reward, stim);
                } 
            }             
            // Data output code placement
            print_WR_ds(NE, WR_ds);
            print_Wds(trial, NE, Wds);
            //Record which trials were rewarded or not
            record_Trials_Rewarded(reward);
            // Attempt at EOut trial version
            eOut_Trial(nBins, bindT, dNPools, eBinRate);
        }
        // Outputs stimuli order over all trials
        // only write out when new stimuli block is delivered
        if (trial == 1 || trial % 4 == 1) { 
            print_Stimuli_Order(inputs, stimShuffle);
        }
        if (debugFlagOn && rasterOutput) {
            // Output Spike Raster in Izhikevich format i.e. 2 columns [Neuron, Spike_time]
            output_Spike_Raster(trial, neurons, dNeurons, numberSpikes, spikeTimes);
        }
        output_Current_Trial_Raster(neurons, dNeurons, numberSpikes, spikeTimes);
        // START OF CUE LENGTH NORMALIZED CUEFR OUTPUT DATA TYPES 
        // 1.) CUEFR OF ECELLS DURING STIM 1 & 2. 
        // 2.) CUEFR OF ICELLS DURING STIM 1 & 2. 
        // 3.) CUEFR OF ECELLS DURING STIM 1
        // 4.) CUEFR OF ECELLS DURING 1ST DELAY
        // 5.) CUEFR OF ECELLS DURING 1ST DELAY
        // 6.) CUEFR OF ECELLS DURING 2ND STIMULUS, WITH DELAY. IF PRESENT.
        // 7, 8.) Mean 1st and 2nd delay persistent firing rate across NE.
        // cueFR append - i.e. one file (tOn[0]-tOff[0] and tOn[1]-tOff[1])
        write_CueFr_1(NE, cueFR);
        // cueFR append - i.e. one file (tOn[0]-tOff[0] and tOn[1]-tOff[1])
        write_CueFr_2(NE, neurons, cueFR);
        // cueFR_1stStim - i.e. activity during 1st stimulus, (tOn[0] - tOff[0])
        write_CueFr_1stStim(NE, cueFR_1stStim);
        // 1st Delay persistent neural activity - i.e. tOff[0] - tOn[1]
        write_CueFr_Perst(NE, cueFR_Perst);
        // 2nd Delay persistent neural activity - i.e. tOff[1] - tMax
        write_CueFr_2nd_Perst(NE, cueFR_2nd_Perst);
        // cueFR_2ndStim_Perst - i.e. activity during 2nd stimulus, with persistent first stimulus activity (tOn[1] - tOff[1])
        write_CueFr_2ndStim_Perst(NE, cueFR_2ndStim_Perst);
        // output mean persistent delay activity and mean 2nd persistent delay activity per trial as -- mean_perstFR append - i.e. one file
        write_Mean_PerstFR(mean_ENeuron_Delay_Rate, mean_ENeuron_2nd_Delay_Rate);
        // END OF CUEFR OUTPUT DATA TYPES 
        // ASSOCIATIVE LAYER BINNED FIRING RATES (100 MS BINS) START
        //AL_Ecellout append trial version
        AL_ECell_BinTrial(nBins, bindT, neurons, AL_ECell_BinRate);
        // Single trial output
        output_AL_ECell_BinRate(trial, nBins, bindT, neurons, AL_ECell_BinRate);
        // END OF ASSOCIATE LAYER BINNED FIRING RATES
        // Output final input/network weight matrices: pW, W
        if (trial == nTrials) {
            print_Final_W_And_pW(neurons, dNeurons, W, pW, nIn);
        }
        // Output weight matrix, Wij, every 10 trials, 1st and last 2.
        if (trial % 40 == 0 || trial == nTrials) {	
            output_Regularly_W(trial, neurons, dNeurons, W);
            output_Regularly_pW(trial, nIn, neurons, dNeurons, pW);
        }
        // Output Voltage data to files
        if (trial == nTrials && voltageDat) {
            output_Voltage_Data(nT, neurons, dNeurons, V);
        }
        // Calculate and print meanWEE
        print_Mean_wEE(NE, W);
        // Calculate and print meanWIE
        print_Mean_wIE(NE, neurons, W);
        // Calculate and print meanWInE
        print_Mean_wInE(nIn, NE, pW);
        // Calculate and print meanWInI
        print_Mean_wInI(nIn, neurons, pW);
    }
    // Print nTotalchange and output nPreChange and nPostChange
    output_N_Changes(nTotChanges, NE, nPreChange, nPostChange);
    // Print final information
    cout << "Finished" << endl;
    time (&rawTime);
    printf ("The current local time is: %s", ctime (&rawTime) );
    return 0;
} 
