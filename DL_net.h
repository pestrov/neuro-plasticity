// Decision Layer parameters
const double sigma = 30e-6;
const double dLayerInput = 0.35; //Not used anywhere now
// Decision layer neurons and associated values
const int dNeurons = 250;  // Number of decision layer neurons 4:1, E:I  
const int dInterneurons = int(dNeurons / 5);
const int dNE = dNeurons - dInterneurons;
const int dNI = dInterneurons;
const int dNPools = 2; // number of pools
const int dNEPool = int(dNE / dNPools);
const int dNIPool = int(dNI / dNPools);