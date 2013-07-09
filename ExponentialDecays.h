// Exponential Decays: Synaptic
const double AMPADecay = exp(-dt / tau_AMPA);
const double GABADecay = exp(-dt / tau_GABA);
const double AL_NMDADecay = exp(-dt / AL_Tau_NMDA);
const double DL_NMDADecay = exp(-dt / DL_Tau_NMDA);