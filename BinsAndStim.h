// Bins
const double bindT = 0.100; // bin size
const double ceilingBins = ceil(tMax / bindT);
const int nBins = int(ceilingBins); // number of bins
// Stim
int stim = 0;
// Vector that contains the permutation of stimulus presentation order: AB, AY, XB, XY
vector<int> stimShuffle(inputs);