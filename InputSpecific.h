//** Input specific parameters
const double gpSyn = 4e-6; // Input "AMPA" conductance
vector<double> tOn(trialStim);// allocate entries for ton
vector<double> tOff(trialStim);

const double gMax_Urgency = 5e-6;//20e-6;  // maximum extra excitatory conductance to DL E-cells

// Cue and delay lengths
double cue1Length = 0;
double cue2Length = 0;
double delay1Length = 0;
double delay2Length = 0;

vector<int> nPreChange(NE);
vector<int> nPostChange(NE);
int nTotChanges = 0;