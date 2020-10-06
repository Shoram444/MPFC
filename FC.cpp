#include <iostream>
#include <math.h> 
#include <vector>
#include <stdio.h>
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "./include/MPFeldman_Cousins.hh"

using namespace std;
R__LOAD_LIBRARY(./lib/libMPFC.so);

int FC() 
{ 

    MPFeldman_Cousins* obj = new MPFeldman_Cousins(6, 0.9);   // order for parameters (double _b, double _step, int _rows, int _mu_max, double _CL)
    obj->print_table  (0.0, 2.5);

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
    // obj->draw_lower();
    // obj->print_n();

    // obj->draw_mu_U_v_b(0);
    // obj->correct_mu_U(2);
    //double b;
    //b = obj->get_b();
    //cout<< "b = " << b<<endl;
    //obj->mu_U_final();
  return 0; 
} 
    
