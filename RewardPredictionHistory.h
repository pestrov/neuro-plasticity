// Reward Prediction Error code - history component
const int nHist = 5; // No. of trial responses to produce "expectation" of reward.
const double histDecay = 0.5;
double sumHistFactor = 0;
vector< vector<int> > decisionHist(nTrials, inputsVec);  // Initialize Reward expectancy history.
vector<double> histFactor(nHist);  // History factor for REWARD PREDICTION ERROR FACTOR INTO DOPAMINE SIGNAL, exponential decay.
//DA Reward constant
const double DA_Epsilon = 0.0001;  // 0.001 is too strong