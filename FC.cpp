#include <iostream>
#include <math.h> 
#include <vector>
#include <stdio.h>
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "../include/MPFeldman_Cousins.hh"

using namespace std;
R__LOAD_LIBRARY(libMPFC.so);

int FC() 
{ 

    MPFeldman_Cousins* obj = new MPFeldman_Cousins(3, 0.01, 80, 50, 0.9);   // order for parameters (double _b, double _step, int _rows, int _mu_max, double _CL)

    // obj->print_poisson();
    // obj->set_mu(0.5);
    // obj->print_R();
    // obj->print_A();
    // std::vector<std::vector<int>> n_array = get_n();
    // obj->print_n(n_array);
    // obj->calculate_upper();
    // cout<< "aaaaa"<< endl;
    // obj->calculate_lower();
    // obj->get_mu_L_v_b(0);
    // obj->draw_upper();
    obj->draw_lower();
    // obj->print_A();
    // obj->draw_mu_U_v_b(2);
    // obj->correct_mu_U(2);

  return 0; 
} 
    