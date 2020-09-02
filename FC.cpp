#include <iostream>
#include <math.h> 
#include <vector>
#include <stdio.h>
// #include "TROOT.h"  
// #include "TGraph.h"
// #include "TCanvas.h"
#include "MPFeldman_Cousins.hh"

using namespace std;

int main() 
{ 
    MPFeldman_Cousins* obj = new MPFeldman_Cousins(3.0, 0.01, 1000, 10, 0.95);   // order for parameters (double _b, double _step, int _rows, int _mu_max, double _CL)

    obj->get_n();
    // obj->draw_upper();
    // cout<< "aaaaa"<< endl;
    // obj->calculate_upper();
    obj->print_A();

  return 0; 
} 
    