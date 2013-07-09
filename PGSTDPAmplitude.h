// PG STDP Amplitude and time constants
const double A2p = 5e-5;//0.0;// 1e-9; // Amplitude of pre-post pair LTP
const double A2m = 7e-3;//7.1e-3;//0.01;// 0.005; // Amplitude of post-pre pair LTD
const double A3p = 6.2e-3;//6.5e-3;//0.01;// 0.005; // Amplitude for triplet LTP
const double A3m = 2.3e-4;//0.0;// 5.0e-4; // Amplitude for triplet LTD
  
const double tau2p = 0.0168;// 0.030; // tau_plus:  Time constant pre-post pair LTP, seconds
const double tau2m = 0.0337;//0.050; // tau_minus:  Time constant post-pre pair LTD
const double tau3p = 0.125;//0.114;//0.200; // tau_y:  Time constant triplet LTP
const double tau3m = 0.101;//0.200; // tau_x:  Time constant triple LTD
  
const double dW = 0.005; // synapse update value