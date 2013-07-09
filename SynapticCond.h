// Synaptic Conductances (AL = Associative layer, DL=Decision Layer)
const double GAL_I_AMPA = 3.0e-6; //1.5e-6;
const double GAL_I_NMDA = 4.0e-6;// 
const double GAL_E_AMPA = 3.0e-6; //1.5e-6;
const double GAL_E_NMDA = 6.0e-6;// 

const double GDL_AMPA = 2.0e-6;
const double GDL_NMDA = 6.0e-6;
const double GDL_I_AMPA = 1.0e-6;
const double GDL_I_NMDA = 1.5e-6;
  
const double mg_Ext = 1;//1.5; // External [Mg2+] mM - Fig 5.16 Abbott & Dayan - range 1-2
const double gSyn_I = 100e-6; // GABA conductance
const double s0 = 0.5; // Synaptic saturation variable