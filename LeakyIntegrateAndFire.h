// Leaky Integrate-and-Fire constants
const double E_s = -70e-3; // Leak reversal
const double tau_M_S = 10e-3; // Membrane time constant
const double vSpike = 30e-3;
const double vReset = -80e-3;  
double tRef = 0.002;
// Refractory Conductance
const double delta_GRef = 2e-3; // increase in refractory conductance per spike
// Max spike constants
const int max_CellSpikes = 400 * int(tMax);
vector<double> vec_Max_CellSpikes(max_CellSpikes);
  
const int max_InputSpikes = 5000;
vector<double> vec_Max_InputSpikes(max_InputSpikes);

// Mean mean firing rate
vector<double> mean_Mean_Firing_Rate(neurons + dNeurons);