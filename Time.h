// Time related values  
const double dt = 0.0002; // time step, .2ms, worked up to .01s
const int tMax = 2;
const int nT = int(tMax / dt);
vector<double> t(nT + 1);