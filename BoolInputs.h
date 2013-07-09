// This file contains the various boolean inputs used throughout the program.

// Homeostasis bools
const bool EE_Homeostasis = 1;
const bool inputEHomeostasis = 1;
const bool inputIHomeostasis = 0;
const bool oneTrialHomeo = 0;

// P&G STDP bools
const bool PG_EESTDP = 1;
const bool PG_Input_ESTDP = 1;
  
//LTPi bool
const bool LTPi = 1;
const bool VLTPi = 1;
  
// DA reward based plasticity
const bool DA_Reward = 1;
const bool DA_Reward_Homeostasis = 1; // not actively used
const bool urgency = 1; // Ramping "Urgency SIgnal" for decision-making
const bool rewardExpectation = 1; // Make DA based on Reward Expectaion Error

// Short-term plasticity bools
const bool facilitation = 1;
const bool depression = 1;
  
// Structural Change bool
const bool structPlast = 1;  

// Structural inputs bools
const bool preShift = 1;
const bool postShift = 1;
const bool ranWStrengths = 1; // heterogeneity

// Input bools
const bool recurrentI = 1;
const bool feedForwardI = 1;
  
// data bools
const bool input_W = 0; // read in pWij and Wij from successful simulation
const bool voltageDat = 0;  // write out a voltage file (large in size)
const bool debugFlagOn = 1;
const bool DL_Out = 1; 

// weight bools
const bool weightsIn = 1;

// If you are testing Decision layer with Poisson Inputs
const bool DL_Poisson_Test = 0;

// Just before the trial loop!
const bool rasterOutput = 1;
const bool plasticityOn = 1;