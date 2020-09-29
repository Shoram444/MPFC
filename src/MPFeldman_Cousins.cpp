#include <iostream>
#include <math.h> 
#include <vector>
#include <bits/stdc++.h> 
#include <chrono>
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "MPFeldman_Cousins.hh"
#include <algorithm>

using namespace std;
ClassImp(MPFeldman_Cousins);

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
	auto start = std::chrono::steady_clock::now();
	table = fill_table();
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Time it took to fill table: " << elapsed_seconds.count() << "s\n";
    set_mu(0);
	R = get_R();
	A = get_A(R);
	n_order = CL_check(mu_j, R);
	n_array = get_n();
	// print_n();



}

MPFeldman_Cousins::~MPFeldman_Cousins()
{

}

double** MPFeldman_Cousins::fill_table()
{
	double** fill_table_table = new double*[ROW_N];

	// table = (double**)malloc(sizeof(double*)*ROW_N);
    for (int i = 0; i <= ROW_N; i++) 
    {
        fill_table_table[i] = new double[COL_N];
    }
    
    
    for (int x = 0; x<=ROW_N;x++)
    {
        for (int y = 0; y<=COL_N;y++)
        {
            int n = x;
            double mu_j = y*STEP;

            
            fill_table_table[x][y] = poisson(n,mu_j);   
            
        }
        
    }

    return fill_table_table;
}

void MPFeldman_Cousins::set_b(double _b)
{
	b = _b;
	return;
}

double MPFeldman_Cousins::get_b()
{
	return b;
}

double* MPFeldman_Cousins::get_R()
{	

	mu_j = get_mu();
	
	double* get_R_R = new double[ROW_N];
	int n_start = 0;
	int n_end = ROW_N;

	for (int i = n_start; i<n_end;i++)
	{
		int mu_best = get_muBest(i, b);

		// cout<< " n_start is : "<< i << "   and int(mu+b)/STEP is: "<< int((mu_j+b)/STEP)<< "   and mub+b is: "<< int((mu_best+b))<< endl;
		// cout<< "table[i][mu+b] "<< table[i][int((mu_j+b)/STEP)]<< "   and table[i][mu_best]:   "<< table[i][int((mu_best+b)/STEP)]<< endl;

		int mu_i = int((mu_j+b)/STEP);
		int mu_b_i = int((mu_best+b)/STEP);

		// cout<<"mu_j = "<< mu_j<<endl;
		// cout<<"mu_b_"<<i<<" = "<< mu_b_i<<endl;
		if(table[i][mu_i] == 0)
		{
			get_R_R[i] = 0;
		}
		else if(table[i][mu_b_i] == 0)
		{
			get_R_R[i] = 0;
		}
		else
		{
			get_R_R[i] = table[i][mu_i]/table[i][mu_b_i];
		}
		// cout<<"R_i is:  "<<R[i]<<endl;

	}
	
	return get_R_R;

}

int* MPFeldman_Cousins::get_A(double* R)
{

	int* get_A_A= new int[ROW_N];

	vector<pair<double, int> > vp; 

	for (int i = 0; i < ROW_N; ++i) 
	{ 
		vp.push_back(make_pair(R[i], i)); 
	} 

	// Sorting pair vector from highest (that's what greater stands for)
	sort(vp.begin(), vp.end(), greater<>()); 

	
  
	for (int i = 0; i < ROW_N; i++)  //saves get_A_A to the array
	{ 
		
		get_A_A[i]=vp[i].second;
	} 

	return get_A_A;

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

	if(ROW_N>1000){ROW_N=1000;}
	if(COL_N>10){COL_N=10;}

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

	if(ROW_N>1000){ROW_N=1000;}
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
	if(ROW_N>1000){ROW_N=1000;}
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

std::vector<int> MPFeldman_Cousins::CL_check(double mu_j, double* _R)
{
	// R = get_R();
	// A = get_A(R);
	double CL_check_mu_j = mu_j;
	double* CL_check_R = _R;
	int* CL_check_A = get_A(CL_check_R);
	std::vector<int> CL_check_n_order;
	double P_Sum = 0;
	int i = 0;
	int max_i = 0;
	std::vector<int> B;
	int check = 0;
	while (P_Sum < CL)	
	{
		if(i == ROW_N - 1)
		{
			cout<< "THE CONDITION OF THE CL WAS NEVER ACHIEVED!! YOU HAVE TO CHOSE LARGER ROW_N!"<<endl;
			check = -1;
			break;
		}
		// cout<<"CL_check_mu_j + b = " <<CL_check_mu_j+b<<endl;
		P_Sum += table[CL_check_A[i]][int((CL_check_mu_j+b)/STEP)];
		// cout<<"P_sum = "<< P_Sum<<"i = "<<i <<endl;

		i += 1;
		// cout<<"CL = "<< CL<< endl;
	
	}
	max_i = i;	

	for(int x = 0; x<=max_i; x++)
	{
		B.push_back(CL_check_A[x]);
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

	if(check !=-1)
	{
		CL_check_n_order.push_back(n_min);
		CL_check_n_order.push_back(n_max);
	}
	else
	{
		CL_check_n_order.push_back(10000);
		CL_check_n_order.push_back(10000);
	}

	B.clear();
	// cout<< " n_min = "<< CL_check_n_order[0]<< "\t n_max = "<< CL_check_n_order[1]<<endl;


	return CL_check_n_order;
}


std::vector<vector<int>> MPFeldman_Cousins::get_n()
{
	std::vector<vector<int>> get_n_n_array;
	for(int col= 0; col< mu_max/STEP ; col++)
	{
		set_mu(col*STEP);
		double get_n_mu_j = get_mu();
		double* get_n_R = get_R();
		// int* get_n_A = get_A(R);
		std::vector<int> get_n_n_order = CL_check(get_n_mu_j, get_n_R);

	
	    std::vector<int> temp;
	    for(int j=0;j<2;j++)
	    {
	        temp.push_back(get_n_n_order[j]);
	    }
	    // cout<<"get_n_mu_j = "<< get_n_mu_j <<"\tn_min\t"<< get_n_n_order[0]<<"\t"<<"\tn_max\t"<< get_n_n_order[1]<<"\t";

	    // cout<<endl;
	    get_n_n_array.push_back(temp);
		
		get_n_n_order.clear();
		get_n_n_order.shrink_to_fit();
		delete[] get_n_R;
		
	}

	return get_n_n_array;
}	
void MPFeldman_Cousins::print_n()
{
	int size = n_array.size();
	for (int i = 0; i < size - 1 ; i++)
	{
		cout<<"mu_j = "<< i*STEP <<"\tn_min\t"<< n_array[i][0]<<"\t"<<"\tn_max\t"<< n_array[i][1]<<"\t";
	    cout<<endl;
	}

}


TGraph* MPFeldman_Cousins::TGraph_calculate_upper() 
{
	auto start = std::chrono::steady_clock::now();
	double mu_array[int(mu_max/STEP)];
	vector<int> n;
	vector<double> mu_U;
	std::vector<vector<int>> calculate_upper_n_array = get_n();

	for (int i = 0; i<mu_max/STEP; i++)
	{
		mu_array[i] = i*STEP;
		
	}


	int size = calculate_upper_n_array.size();
	int n_current = calculate_upper_n_array[size - 1][0];
	// cout<< "n_current = "<< n_current<< endl;
	// int i = size -1;

	int y = 0;
	for (int x = size -1 ; x > 0; x--)  			//algortihm checks the n_array generated in get_n from back to front. 
													//Logs only jumps by one, so that there are no repetitions.
	{
		
		if (n_current > calculate_upper_n_array[x][0])
		{

			if(calculate_upper_n_array[x][0] == 10000)
			{
				continue;

			}
			else
			{
				n.push_back(n_current);
			
				mu_U.push_back(mu_array[x+1]);
				n.push_back(n_current - 1 );
				
			
				mu_U.push_back(mu_array[x+1]);
			
				n_current = calculate_upper_n_array[x][0];
				y+=2;
			}
			
		}
		
	}
	n.push_back(0);					//Adding point (0,0) so that the graph starts at 0. 
	mu_U.push_back(0);

	int mu_U_size = mu_U.size();
	reverse(mu_U.begin(),mu_U.end());
	reverse(n.begin(),n.end());

	TGraph *gr  = new TGraph();


	for (int i =0; i<mu_U_size; i++)
	{
		gr->SetPoint(i, n[i], mu_U[i]);
		// cout<< "n_" << i << " = " << n[i] << "\t mu_U[" << i << "]  = "<< mu_U[i] << endl;
	}
  	
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Time it took to calculate_upper: " << elapsed_seconds.count() << "s\n";
	// delete[] mu_array;
	return gr;
}

void MPFeldman_Cousins::draw_upper()
{

	TGraph* gr = new TGraph();
	gr = TGraph_calculate_upper();

	stringstream graph_title;
    graph_title  << "Upper Limit for bkg = " << b << " with C.L. = " << CL*100 << "%";
    string strname  = graph_title.str();

 	gr->SetTitle(strname.c_str());
	gr->SetFillStyle(1000);
	gr->SetLineColor(2);
	gr->SetLineWidth(4);
	gr->SetMarkerStyle(0);
	// gr->SetMaximum(20);
	TAxis* xaxis;
	xaxis = gr->GetXaxis();
	// xaxis->SetLimits(0,20);
	xaxis->CenterTitle(true);
	xaxis->SetTitle("n [Counts]");

	TAxis* yaxis;
	yaxis = gr->GetYaxis();
	yaxis->SetTitle("\\mu_U  [Counts]");
	yaxis->CenterTitle(true);


	TCanvas *c1 = new TCanvas("c1","Graph Draw Options",600,600);


	c1->SetFrameFillStyle(3013);
	c1->SetFrameLineWidth(2);
	c1->SetFrameBorderMode(0);
	c1->SetFrameBorderSize(4);
	c1->SetGridx();
	c1->SetGridy();


	gr->Draw("ALP");

	return;
}	



vector<double> MPFeldman_Cousins::calculate_upper() ///CHANGE TO GET GRAPH!!!
{
	auto start = std::chrono::steady_clock::now();
	double* mu_array = new double[COL_N];
	vector<int> n;
	vector<double> mu_U;
	std::vector<vector<int>> calculate_upper_n_array = get_n();

	for (int i = 0; i<COL_N; i++)
	{
		mu_array[i] = i*STEP;
		
	}


	int size = calculate_upper_n_array.size();
	int n_current = calculate_upper_n_array[size - 1][0];
	// cout<< "n_current = "<< n_current<< endl;
	// int i = size -1;
				
	int y = 0;
	for (int x = size -1 ; x >= 0; x--)
	{
		
		if (n_current > calculate_upper_n_array[x][0])
		{

			if(calculate_upper_n_array[x][0] == 10000)
			{
				continue;

			}
			else
			{
				n.push_back(n_current);
			
				mu_U.push_back(mu_array[x+1]);
				n.push_back(n_current - 1 );
				
			
				mu_U.push_back(mu_array[x+1]);
			
				n_current = calculate_upper_n_array[x][0];
				y+=2;
			}
			
		}
		
	}


	// int entries = y;

	// for (int i =0 ; i<y; i++)
	// {
	// 	cout<< "n = " << n[i] << "\t mu_U = " << mu_U[i] << endl;
	// 	// cout<< "Got to printing"<<endl;


	// }
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Time it took to calculate_upper: " << elapsed_seconds.count() << "s\n";
	delete[] mu_array;

	return mu_U;
}

vector<double>  MPFeldman_Cousins::calculate_lower() 
{
	
	double mu_array[COL_N];
	vector<int> n;
	vector<double> mu_L;
	std::vector<vector<int>> calculate_lower_n_array = get_n();


	for (int i = 0; i<COL_N; i++)
	{
		mu_array[i] = i*STEP;
		
	}

	int size = calculate_lower_n_array.size();
	int n_current = calculate_lower_n_array[size - 1][1];
	// int i = size -1;
	// cout<< "n_current" << n_current<<endl;
	// cout<< "i " << i <<endl;

	int y = 0;
	
	for (int x = size -1 ; x >= 0; x--)
	{
		
		if (n_current > calculate_lower_n_array[x][1])
		{
			n.push_back(n_current);
		
			mu_L.push_back(mu_array[x+1]);
			n.push_back(n_current - 1 );
			
		
			mu_L.push_back(mu_array[x+1]);
		
			n_current = calculate_lower_n_array[x][1];
			y+=2;
		}
		
	}


	// int entries = y;

	// TGraph *gr  = new TGraph();

	// for (int i =0; i<y; i++)
	// {
	// 	gr->SetPoint(i, n[i], mu_L[i]);
	// }
  	for (int i =0 ; i<y; i++)
	{
		cout<< "n = " << n[i] << "\t mu_L = " << mu_L[i] << endl;

	}

	return mu_L;
}



double* MPFeldman_Cousins::get_mu_U_v_b(int n)
{

	double* bkg      = new double[ROW_N];
	double* mu_U_v_b = new double[ROW_N];

	for (int i = 0; i<ROW_N; i++)
	{
		// cout<< "started loop get mu v b"<< endl;
		auto start = std::chrono::steady_clock::now();

		set_b(i*STEP*20);

		bkg[i] = get_b();

		// cout<< "got to bkg = "<< bkg[i]<< endl;
		vector<double> mu_U;
		mu_U = calculate_upper();

		// cout<< "got to calc upper "<< endl;
		// print_upper();

		if (n == 0)
		{
			int index = mu_U.size() - 1;
			mu_U_v_b[i] = mu_U[index];
		}
		else
		{
			int index = mu_U.size() - 1;

			mu_U_v_b[i] = mu_U[index-2*n];
		}
		cout<< "b = " << bkg[i] << "\t mu_u = "<< mu_U_v_b[i]<<endl;

		// vector<double>().swap(mu_U);
		// mu_U.shrink_to_fit();
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
	    std::cout << "Time it took to find mu vs b: " << elapsed_seconds.count() << "s\n";
	}

	return mu_U_v_b;
} 

void MPFeldman_Cousins::get_mu_L_v_b(int n)
{

	double* bkg      = new double[ROW_N];
	double* mu_L_v_b = new double[ROW_N];


	for (int i = 0; i<ROW_N; i++)
	{
		// cout<< "started loop get mu v b"<< endl;
		auto start = std::chrono::steady_clock::now();

		set_b(i*STEP);

		bkg[i] = get_b();

		// cout<< "got to bkg = "<< bkg[i]<< endl;
		vector<double> mu_L;
		mu_L = calculate_upper();

		// cout<< "got to calc upper "<< endl;
		// print_upper();

		if (n == 0)
		{
			int index = mu_L.size() - 1;
			mu_L_v_b[i] = mu_L[index];
		}
		else
		{
			int index = mu_L.size() - 1;

			mu_L_v_b[i] = mu_L[index-2*n];	
		}
		cout<< "b = " << bkg[i] << "\t mu_L = "<< mu_L_v_b[i]<<endl;

		mu_L.clear();
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
	    std::cout << "Time it took to find mu vs b: " << elapsed_seconds.count() << "s\n";
	}
	return;
} 

void MPFeldman_Cousins::draw_mu_U_v_b(int n)
{
	double* mu_U_v_b = new double[ROW_N];
	double* bkg      = new double[ROW_N];

	mu_U_v_b = get_mu_U_v_b(n);
	for (int i = 0; i < ROW_N; ++i)
	{
		bkg[i] = i*STEP*20;
	}

	TGraph* gr = new TGraph(ROW_N,bkg,mu_U_v_b);

	stringstream graph_title;
    graph_title  << "Upper Limit for n = " << n << " with C.L. = " << CL*100 << "%";
    string strname  = graph_title.str();

 	gr->SetTitle(strname.c_str());
	gr->SetFillStyle(1000);
	gr->SetLineColor(2);
	gr->SetLineWidth(4);
	gr->SetMarkerStyle(0);
	// gr->SetMaximum(20);
	TAxis* xaxis;
	xaxis = gr->GetXaxis();
	// xaxis->SetLimits(0,20);
	xaxis->CenterTitle(true);
	xaxis->SetTitle("b [Counts]");

	TAxis* yaxis;
	yaxis = gr->GetYaxis();
	yaxis->SetTitle("\\mu_U  [Counts]");
	yaxis->CenterTitle(true);


	TCanvas *c1 = new TCanvas("c1","Mu_u vs background",600,600);


	c1->SetFrameFillStyle(3013);
	c1->SetFrameLineWidth(2);
	c1->SetFrameBorderMode(0);
	c1->SetFrameBorderSize(4);
	c1->SetGridx();
	c1->SetGridy();


	gr->Draw("ALP");
	return;
}
