#include <iostream>
#include <math.h> 
#include <vector>
#include <stdio.h>
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "TMath.h"
#include "./include/MPFeldman_Cousins.hh"

using namespace std;
R__LOAD_LIBRARY(./lib/libMPFC.so);

int FC() 
{ 

    MPFeldman_Cousins* obj = new MPFeldman_Cousins(3, 0.9);   // order for parameters (double _b, double _step, int _rows, int _mu_max, double _CL)
    // obj->print_table  (0.0, 2.5);

    obj->fill_m_table(10);
    // obj->calculate_lim();
    // obj->shift_mu_U(3.0, 10);

    // cout<< "second object: "<<endl;
    // MPFeldman_Cousins* obj2 = new MPFeldman_Cousins(4, 0.9);   // order for parameters (double _b, double _step, int _rows, int _mu_max, double _CL)
    // obj2->calculate_lim();

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
    
    // double p;
    // double p1;
    // double ratio;
    // double max_error = 0;
    // for(int i = 0; i < 100000; i++)
    // {
    //     p = obj->poisson(1,i*0.001);
    //     p1 = TMath::Poisson(1,i*0.001);
    //     ratio = abs((p-p1)/p1);
    //     if(ratio > 0.0001)
    //     {
    //         cout<< "==> p from FC = " << p;

    //         cout<< "<====> p from TMath = "<<p1 << endl;

    //         cout<< " ===== Ratio = "<< ratio << "=====" <<endl;

    //     }
    //     else
    //     {
    //         cout<< " ===== Ratio = "<< ratio << "=====" <<endl;

    //     }
    //     if(max_error < ratio)
    //     {
    //         max_error = ratio;
    //     }

    // }
    // cout<< " Max error = " << max_error << endl;