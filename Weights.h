// Weights
const double wEE = 1.0;
const double wEI = 50.0;//90.0;
const double wIE = 8.0;//10.0;
const double wII = 2.0;
//Network Weight Matrices and values
const double wEE_ConnProb = 0.1;
const double wIE_ConnProb = 0.25;
const double wEI_ConnProb = 0.25;//0.25;//0.1;
const double wII_ConnProb = 1.0;
const double wEE_xFactor = 2.50;
// We'll scale synaptic strength by network size.
const double wEE0_Default = 100;
const double wEE0 = 0.05;//100 / double(neurons);
const double wIE0 = 0.05;//(1.75 * wEE0_Default) / double(neurons);//0.2*WEE0;
const double wEI0 = (1.75 * wEE0_Default) / double(neurons);// 0.5*WEE0;
const double wII0 = 0.038;//(0.5 * wEE0_Default) / double(neurons);   //1.0*W0;
const double W0 = wEE0;  // base synaptic strength
const double connProb = wEE_ConnProb;

// Weight vectors
vector< vector<double> > W(neurons + dNeurons, DDNeuronVec);
vector< vector<double> > dWstruct(neurons + dNeurons, DDNeuronVec);
vector< vector<double> > W1(neurons + dNeurons, DDNeuronVec); // Normalization Matrix used by STDP code

// Initialize input weights
const double pW0 = 1.00;  // Excitatory input weights
const double pW0I = 1.00; // Inhibitory input weights

// Only used for Poisson test
const double pool1Weight = 1.10;  // bottom of .85, no lower
const double pool2Weight = 1.00;
const double stimBiasPerc =  (1 - 1.0 / (pool1Weight / pool2Weight)) * 100.0;

// Decision layer weights
const double dW0 = 1 * 5 / double(dNeurons);  // Decision layer within weights
const double dLayerInput_W0 = 0.088;//dLayerInput * W0;

// bool matrices
const int totalNeurons = neurons + dNeurons; // This is used when passing C1 as a parameter
bool C1[neurons + dNeurons][neurons + dNeurons]; //THE MATRIX, showing if the synapse is active or not