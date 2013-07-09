// Vectors reset at the beginning of each trial
vector< vector<double> > spikeTimes(neurons + dNeurons, vec_Max_CellSpikes);
vector<int> numberSpikes(neurons + dNeurons);
vector<double> gSyn_I_Vec(neurons + dNeurons);
vector<double> gSyn_E_Vec(neurons + dNeurons);
vector<double> gpSyn_Vec(neurons + dNeurons);
vector<double> vInf(neurons + dNeurons);
vector< vector<double> > V(neurons + dNeurons, t);	
vector<double> lastSpike(neurons + dNeurons);
vector<double> syn(neurons + dNeurons);
vector<double> syn_I(neurons + dNeurons);
vector<double> syn_E(neurons + dNeurons);

vector<double> gSyn_NMDA(neurons + dNeurons);
vector<double> g_NMDA(neurons + dNeurons);
vector<double> Vth(neurons + dNeurons);
vector<double> syn_NMDA(neurons + dNeurons);
vector<double> g_Ref(neurons + dNeurons);	
vector<double> meanFiringRate(neurons + dNeurons);
vector<double> W_Shift(neurons + dNeurons);
vector<double> W_Shift_I(NE);
vector<double> DW_Shift(dNeurons);
    
vector<double> Fac(neurons + dNeurons);
vector<double> synE_Fac(neurons + dNeurons);
vector<double> synNMDA_Fac(neurons + dNeurons);
vector<double> Dep(neurons + dNeurons);
vector<double> synE_Dep(neurons + dNeurons);
vector<double> synNMDA_Dep(neurons + dNeurons);

double g_Tot;
double tau_Eff;
double pSynMax;
double synMax_I;
double synMax_E;
double synMax_NMDA;
    
int someIndex;
