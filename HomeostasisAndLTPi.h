//** Homeostasis and LTPi parameters
const double LTPi_Window = 0.020; // LTPi window size
const double VLTPi_Threshold = -65e-3; // Depolarization threshold for LTPi
const double idW = 0.001;//0.005; // change per non-vetoed spike
const double e_E = 0.001;// Homeostasis rate constant
const double e_I = 0.001;//(.005/10.0); // Homeostasis rate constant TOASK: the same as in excititary?

const double e_Input_I = 0.001;
const double eD = 0.001;