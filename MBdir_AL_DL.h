// Size of network and number of trials
const int nTrials = 200;
const int neurons = 600;
const int interneurons = int(neurons / 5);
const int NE = neurons - interneurons;
const int NI = interneurons;
const int nPools = 16; // number of pools
const int inputs = 4; // number of stimul
vector<int> inputsVec(inputs); // vector for inputs
const int trialStim = 2; // Number of inputs per trial 

// DL Network noise amplitudes
//const double sigma = 20.0e-6; // DL Voltage noise amplitude
const double sigma_E = 4e-3;// 1.1e-3; // DL_E conductance noise amplitude
const double sigma_I = 4e-3;// 1.1e-3; // DL_I conductance noise amplitude
// AL Network noise values
const double nSigma = 60.0e-6;  // Associative layer Voltage noise amplitude
const double nSigma_E = 1.2e-3;//1.2e-3; // AL_E conductance noise amplitude
const double nSigma_I = 1.2e-3;//1.2e-3; // AL_I conductance noise amplitude

//Goal rates
const int rgE = 8;  // excitatory
const int rgI = 10; // inhibitory
const int rgD = 20;   // Decision layer Ecells
const double rgIE0 = 8 - rgE;// (rgI*tmax*exp(-LTPi_window*rgE)*idW)/e;  // currently hand-tuned.  subtract rgE since it's added below.

const double g_L_s = 35e-6;//35e-6;//23e-6; // Leak conductance - increase -> harder to fire
const double AL_g_L_s = 35e-6;//35e-6;//23e-6; // Leak conductance - increase -> harder to fire
const double pW02 = 4;//8.0; //2max // how much to amplify the read-in pW by.
const double pW0I2 = 4;//10.0;//3.50--12.0 worked but, made too strong V deflections // how much to amplify the read-in pW by.

vector<double> DDNeuronVec(neurons + dNeurons);
vector<double> NEVec(NE);

//Goal rate vectors
vector<double> r_Goal(neurons);
vector<double> r_Goal_I_to_E(neurons);
vector<double> r_DGoal(dNeurons);

