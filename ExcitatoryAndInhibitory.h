// Excitatory
const double vth_E_sA = -50e-3; // -45 Excitatory Vth
const double vth_E_sD = -48e-3; // -45 Excitatory Vth
const double VReset_E_s = -58e-3; // Excitatory reset
const double TRef_E_s = 0.00225; // Excitatory tau_ref
// Inhibitory
const double vth_I_sA = -48e-3; // Inhibitory Vth
const double vth_I_sD = -50e-3; // Inhibitory Vth
const double vReset_I_s = -58e-3; // Inhibitory Vreset
const double tRef_I_s = 0.00125; // Inhibitory tau_ref
// Persistence & NMDA specific
double vth0 = -45e-3; // same as Vth_e, may need to be changed.  We actually alter this in our program, so
                      // it can't be const.
const double vth_Max = 150e-3; // max Vth 