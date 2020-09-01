#include <iostream>
#include <math.h> 
#include <vector>
#include <bits/stdc++.h> 
// #include "TROOT.h"  
// #include "TGraph.h"
// #include "TCanvas.h"
#include "MPFeldman_Cousins.hh"

using namespace std;


MPFeldman_Cousins::MPFeldman_Cousins()
{

}


MPFeldman_Cousins::MPFeldman_Cousins(double _b, double _step, int _rows, double _mu_max, double _CL)
{

	STEP  = _step;  //step at which mu iterates. 
	ROW_N = _rows;	//Number of n.
	mu_max = _mu_max;
	COL_N = _rows/STEP;  //The number of columns dependent on max mu and step chosen.
	b = _b;
	CL = _CL;


	table = (double**)malloc(sizeof(double*)*ROW_N);
    for (int i = 0; i <= ROW_N; i++) 
    {
        table[i] =(double*)malloc(sizeof(double)*COL_N);
    }
    
    
    for (int x = 0; x<=ROW_N;x++)
    {
        for (int y = 0; y<=COL_N;y++)
        {
            int n = x;
            double mu_j = y*STEP;

            
            table[x][y] = poisson(n,mu_j);   
            
        }
        
    }

}

MPFeldman_Cousins::~MPFeldman_Cousins()
{

}

double* MPFeldman_Cousins::get_R()
{	
	// double mu_j = Get_mu();
	
	R = (double*)malloc(sizeof(double)*ROW_N);
	int n_start = 0;
	int n_end = ROW_N;

	for (int i = n_start; i<n_end;i++)
	{
		int mu_best = get_muBest(i, b);

		// cout<< " n_start is : "<< i << "   and int(mu+b)/STEP is: "<< int((mu_j+b)/STEP)<< "   and mub+b is: "<< int((mu_best+b))<< endl;
		// cout<< "table[i][mu+b] "<< table[i][int((mu_j+b)/STEP)]<< "   and table[i][mu_best]:   "<< table[i][int((mu_best+b)/STEP)]<< endl;

		int mu_i = int((mu_j+b)/STEP);
		int mu_b_i = int((mu_best+b)/STEP);

		// cout<<"mu_"<<i<<" = "<< mu_i<<endl;
		// cout<<"mu_b_"<<i<<" = "<< mu_b_i<<endl;
		if(table[i][mu_i] == 0)
		{
			R[i] = 0;
		}
		else if(table[i][mu_b_i] == 0)
		{
			R[i] = 0;
		}
		else
		{
			R[i] = table[i][mu_i]/table[i][mu_b_i];
		}
		// cout<<"R_i is:  "<<R[i]<<endl;

	}
	
	return R;

}

int* MPFeldman_Cousins::get_A(double* R)
{
	int* indexes= (int*)malloc(sizeof(int)*ROW_N);
	// Vector to store element 
	// with respective present index 
	vector<pair<double, int> > vp; 

	// Inserting element in pair vector 
	// to keep track of previous indexes 
	for (int i = 0; i < ROW_N; ++i) 
	{ 
		vp.push_back(make_pair(R[i], i)); 
	} 

	// Sorting pair vector 
	sort(vp.begin(), vp.end(), greater<>()); 

	// Displaying sorted element 
	// with previous indexes 
	// corresponding to each element 

  
	for (int i = 0; i < ROW_N; i++) 
	{ 
		// cout<< "vp first (R)"<< vp[i].first<<"    vp secound (index)"<< vp[i].second<<endl;
		// cout<< "R i "<< R[i]<<endl;
		indexes[i]=vp[i].second;
	} 

	// cout<< "index 4 = "<< indexes[4]<<endl;
	return indexes;

}


double MPFeldman_Cousins::get_muBest(int n, double b)
{

	double mu_bst = double(n) - b;		

	if (mu_bst < 0)
	{
		mu_bst = 0;
	}

	return mu_bst;
}


double MPFeldman_Cousins::poisson(int n, double mu_j)
{

	double const e_to_mu_n  = exp(mu_j / double(n));
	double P = 1.0;
	
	if (n == 0)
	{
		if (mu_j == 0.0)
		{
			P = 1.0;
		}
		else
		{
			P = exp(-mu_j); // exp(-mu_j)*(mu_j^0/0!) = exp(-mu_j)
		}
	}
	else if (mu_j == 0.0 && n != 0)
	{
		P = 0.0;
	}
	else
	{
		int no_multip  = 0;
		int factor_div = 1;
		for (int i = 1; i <= 2*n; i++)	
		{
				if ( P < 1.0 && no_multip < n)
				{ 
					P *= mu_j;
					no_multip++;
				}
				else
				{				
					P /= factor_div*e_to_mu_n;
					factor_div++;
				}
		}
	}
	if (P < 1e-310)
	{
		P = 0;
	}

    return P;
    
}

void MPFeldman_Cousins::print_poisson()
{

	if(ROW_N>50){ROW_N=50;}
	if(COL_N>0){COL_N=10;}

	for(int i = 0; i<COL_N;i++)
	{
		cout<<"\tmu = "<<i*STEP<<" ";
	}
	cout<<endl;

	for (int x = 0; x<=ROW_N;x++)
    {
    	cout<<"n = "<<x<<"\t";
        for (int y = 0; y<=COL_N;y++)
        {
            cout<<table[x][y]<< "\t";
        }
        cout<<endl;
    }

}

void MPFeldman_Cousins::print_R()
{	
	R = get_R();

	if(ROW_N>100){ROW_N=100;}
	for (int i = 0; i<ROW_N;i++)
	{
		cout<< "n = "<< i << "\t R = " << R[i]<<"\tPoisson i = " <<table[i][int((mu_j+b)/STEP)] <<endl;
	}
}

void MPFeldman_Cousins::set_mu(double _mu_j)
{
	mu_j = _mu_j;
}

double MPFeldman_Cousins::get_mu()
{
	return mu_j;
}

void MPFeldman_Cousins::print_A()
{
	R = get_R();
	A = get_A(R);
	if(ROW_N>50){ROW_N=50;}
	// cout<<"mu_j = "<< mu_j << "\tb = "<<b<<endl;
	cout << setw(4) << "n" << "|" <<setw(12) <<"P(n, mu+b)  |" << setw(4)<< " mu |" <<setw(16) <<"  P(n, mu_b+b)  |" <<setw(13)<<"R     |" << setw(5)<< "A|"<<endl;
	cout << "----|------------|----|----------------|------------|----|" << endl;
	for(int i =0; i<ROW_N;i++)
	{
		cout << setw(4) << i 
		     << "|" << setw(12) << table[i][int((mu_j+b)/STEP)]
		     << "|" << setw(4)  << get_muBest(i, b) 
		     << "|" << setw(16) << table[i][int(get_muBest(i, b)+b/STEP)] 
		     << "|" << setw(12) << R[i] 
		     << "|" << setw(4)  << A[i] << "|" <<endl;
	}

}

std::vector<int> MPFeldman_Cousins::CL_check(double mu_j)
{
	R = get_R();
	A = get_A(R);
	
	double P_Sum = 0;
	int i = 0;
	int max_i = 0;
	std::vector<int> B;

	while (P_Sum < CL)	
	{
		if(i == ROW_N - 1)
		{
			cout<< "THE CONDITION OF THE CL WAS NEVER ACHIEVED!! BREAK!"<<endl;

			break;
		}
		// cout<<"mu_j + b = " <<mu_j+b<<endl;
		P_Sum += table[A[i]][int((mu_j+b)/STEP)];
		// cout<<"P_sum = "<< P_Sum<<"i = "<<i <<endl;

		i += 1;
		// cout<<"CL = "<< CL<< endl;
	
	}
	max_i = i;	

	for(int x = 0; x<=max_i; x++)
	{
		B.push_back(A[x]);
	}

	int n_min = B[0];
	int n_max = B[0];
	for (int i = 0; i< max_i; i++)
	{
		if(n_min > B[i])
		{
			n_min = B[i];
		}
	}
	for (int i = 0; i< max_i; i++)
	{
		if(n_max < B[i])
		{
			n_max = B[i];
		}
	}

	n_order.push_back(n_min);
	n_order.push_back(n_max);

	B.clear();
	// cout<< " n_min = "<< n_order[0]<< "\t n_max = "<< n_order[1]<<endl;


	return n_order;
}


std::vector<vector<int>> MPFeldman_Cousins::get_n()
{
	
	for(int col= 0; col< 2*mu_max ; col++)
	{
		set_mu(col*STEP);
		double mu_j = get_mu();
		R = get_R();
		A = get_A(R);
		n_order = CL_check(mu_j);

	
	    std::vector<int> temp;
	    for(int j=0;j<2;j++)
	    {
	        temp.push_back(n_order[j]);
	    }
	    cout<<"mu_j = "<< mu_j <<"\tn_min\t"<< n_order[0]<<"\t"<<"\tn_max\t"<< n_order[1]<<"\t";

	    cout<<endl;
	    n_array.push_back(temp);
		
		n_order.clear();
		
	}
	return n_array;
}	


// TGraph* MPFeldman_Cousins::calculate_upper() ///CHANGE TO GET GRAPH!!!
// {
// 	double mu_array[COL_N];
// 	vector<int> n;
// 	vector<double> mu_U;

// 	for (int i = 0; i<COL_N; i++)
// 	{
// 		mu_array[i] = i*STEP;
		
// 	}

// 	int size = n_array.size();
// 	int n_current = n_array[size - 1][0];
// 	int i = size -1;

// 	int y = 0;
// 	for (int x = size -1 ; x >= 0; x--)
// 	{
		
// 		if (n_current > n_array[x][0])
// 		{
// 			n.push_back(n_current);
		
// 			mu_U.push_back(mu_array[x+1]);
// 			n.push_back(n_current - 1 );
			
		
// 			mu_U.push_back(mu_array[x+1]);
		
// 			n_current = n_array[x][0];
// 			y+=2;
// 		}
		
// 	}


// 	int entries = y;

// 	TGraph *gr  = new TGraph();

// 	for (int i =0; i<y; i++)
// 	{
// 		gr->SetPoint(i, n[i], mu_U[i]);
// 	}
  	

// 	return gr;
// }



void MPFeldman_Cousins::calculate_upper() ///CHANGE TO GET GRAPH!!!
{
	double mu_array[COL_N];
	vector<int> n;
	vector<double> mu_U;

	for (int i = 0; i<COL_N; i++)
	{
		mu_array[i] = i*STEP;
		
	}

	int size = n_array.size();
	int n_current = n_array[size - 1][0];
	// int i = size -1;

	int y = 0;
	for (int x = size -1 ; x >= 0; x--)
	{
		
		if (n_current > n_array[x][0])
		{
			n.push_back(n_current);
		
			mu_U.push_back(mu_array[x+1]);
			n.push_back(n_current - 1 );
			
		
			mu_U.push_back(mu_array[x+1]);
		
			n_current = n_array[x][0];
			y+=2;
		}
		
	}


	// int entries = y;

	for (int i =0 ; i<y; i++)
	{
		cout<< "n = " << n[i] << "\t mu_U = " << mu_U[i] << endl;

	}

	return;
}
// void MPFeldman_Cousins::calculate_lower() ///CHANGE TO GET GRAPH!!!
// {
// 	double mu_array[COL_N];
// 	vector<int> n;
// 	vector<double> mu_L;

// 	for (int i = 0; i<COL_N; i++)
// 	{
// 		mu_array[i] = i*STEP;
		
// 	}

// 	int size = n_array.size();
// 	int n_current = n_array[size - 1][1];
// 	int i = size -1;
// 	// cout<< "n_current" << n_current<<endl;
// 	// cout<< "i " << i <<endl;

// 	int y = 0;
// 	for (int x = size -1 ; x >= 0; x--)
// 	{
		
// 		if (n_current > n_array[x][1])
// 		{
// 			n.push_back(n_current);
// 			cout<< " n = " << n[y]<< "\t";
		
// 			mu_L.push_back(mu_array[x+1]);
// 			cout<< "mu_L = "<< mu_L[y]<<endl;
		
// 			n_current = n_array[x][1];
// 			y+=1;
// 		}
		
// 	}


// 	// int entries = y;

// 	// TGraph *gr  = new TGraph();

// 	// for (int i =0; i<y; i++)
// 	// {
// 	// 	gr->SetPoint(i, n[i], mu_L[i]);
// 	// }
  	

// 	return; //gr;
// }


// void MPFeldman_Cousins::draw_upper()
// {
// 	TGraph *gr  = new calculate_upper();
// 	TCanvas *c1 = new TCanvas("c1","Graph Draw Options",
//                              600,600);

// 	gr->Draw("AC*");

// 	return;
// }										// draws upper limit
