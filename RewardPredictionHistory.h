// Reward Prediction Error code - history component
const int historyLength = 5; // No. of trial responses to produce "expectation" of reward.
const double histDecay = 0.5;
double histFactorSum = 0;
vector< vector<int> > decisionHist(nTrials, inputsVec);  // Initialize Reward expectancy history.
vector<double> histFactor(historyLength);  // History factor for REWARD PREDICTION ERROR FACTOR INTO DOPAMINE SIGNAL, exponential decay.
//DA Reward constant
const double DA_Epsilon = 0.0001;  // 0.001 is too strong