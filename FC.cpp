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

    MPFeldman_Cousins* obj = new MPFeldman_Cousins(21.0, 0.9);   // generated bmax, CL, 

    obj->print_table(0.0, 10.0, 0, 43, 1.0, true);
    obj->print_table(0.0, 10.0, 0, 43, 1.0, false);
    
	for (int i = 0; i < 20; i++)
	{	
		cout << i << ": " << obj->get_sensitivity(i) << endl;
	}

  return 0; 
} 
