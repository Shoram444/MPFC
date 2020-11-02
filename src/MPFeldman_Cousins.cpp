#include <iostream>
#include <math.h> 
#include <vector>
#include <bits/stdc++.h> 
#include <chrono>
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"
#include "MPFeldman_Cousins.hh"
#include <algorithm>

using namespace std;
ClassImp(MPFeldman_Cousins);

MPFeldman_Cousins::MPFeldman_Cousins()
{
}

MPFeldman_Cousins::MPFeldman_Cousins(double _b, double _CL)
{
	table_set = false;
	b         = _b;
	CL        = _CL;

	auto start = chrono::steady_clock::now();
		fill_table();
	auto end   = chrono::steady_clock::now();
	chrono::duration<double> elapsed_seconds = end-start;
    
	cout << "Time it took to fill table: " << elapsed_seconds.count() << "s\n";
/*

	b_const = _b;
	     table = fill_table();
    	
	set_mu(0);

	R       = get_R();
	A       = get_A(R);
	n_order = CL_check(mu_j, R);
	n_array = get_n();*/
	// print_n();
}

MPFeldman_Cousins::~MPFeldman_Cousins()
{

}

double MPFeldman_Cousins::calculate_lim()
{
	int 	cycles = 2000; 
	double* mu_U_temp = new double[cycles];
	double* mu_L_temp = new double[cycles];
	int*	n_U_temp  = new int[cycles];
	int* 	n_L_temp  = new int[cycles];




	for (int m = 0; m<cycles; m ++)
	{
		double mu 			= m*STEP;
		double P_Sum 		= 0.0;
		int    n_min        = 0;
		int    n_max        = 0;
		int	   checker 		= 0;
		int    iter;
		bool   n_3_used     = false;

		int 	n_top	 	= ceil( mu + b ) ; 

		double  mu_bst; 
		double  mub_highest;

		double* buffer  =	 new double[7]; //Prepare buffer for R values of size 7
		double* buf_n   =    new double[7]; //Prepare buffer of corresponding n values associated to buffer

		for(int i = 0 ; i < 7; i++) //Filling out buffer around the n_top value. 
		{
			if (n_top - 3 + i < 0)  //Makes sure that the buffer of n values doesn't go to negatives 
			{
				buf_n[i]  = 	 0;
				buffer[i] = 	-1; 
			}
			else
			{
				buf_n[i]  = 	n_top - 3 + i;
				mu_bst    = 	get_muBest ( buf_n[i] , b );
				buffer[i] = 	get_R(buf_n[i] , mu + b , mu_bst + b );
			}

		}

		double max 			= 	buffer[0]; //max is used to search for highest R value
		int    n_highest 	= 	buf_n[0];  //n_highest is the corresponding n for highest R

		for (int j = 0 ; j <7 ; j++)       //find the highest R value from buffer and it's corresponding n
		{
			if(max <= buffer[j])
			{
				max 		= 	buffer[j];
				n_highest 	= 	buf_n[j];
				mub_highest =	get_muBest(n_highest , b);
				iter = j;					// iter is used to make sure that the R value that has been used once will not be used again
			}
		}

		n_min 	= 	n_highest;			//In first round of searching for n_min and n_max both are set to default first highest value
		n_max 	= 	n_highest; 
		buffer[iter] = -1;
		P_Sum += poisson(n_highest , mu + b);



		if( buf_n[0] != 0 && (poisson(buf_n[0] - 1, mu ) > poisson(buf_n[6] + 1, mu )) ) //this condition makes sure that when moving buffer toward lower numbers, we do not get out of bounds.
		{
			buf_n[0] 	 = buf_n[0] - 1;
			mu_bst    	 = get_muBest ( buf_n[0] , b );
			buffer[0] 	 = get_R(buf_n[0] , mu + b , mu_bst + b );


			for (int i = 1 ; i<4 ; i++)
			{

				if(buf_n[i] - 1 <0 )
				{
					buf_n[i]  = 	0;
					mu_bst    = 	get_muBest ( buf_n[i] , b );
					buffer[i] = 	get_R(buf_n[i] , mu + b , mu_bst + b );
				}
				else
				{
					buf_n[i]  = 	buf_n[i] - 1;
					mu_bst    = 	get_muBest ( buf_n[i] , b );
					buffer[i] = 	get_R(buf_n[i] , mu + b , mu_bst + b );
				}
			}
		}
		else 
		{
			buf_n[3] 	 = buf_n[4];
			mu_bst    	 = get_muBest ( buf_n[3] , b );
			buffer[3] 	 = get_R(buf_n[3] , mu + b , mu_bst + b );

			for (int i = 4 ; i<7 ; i++)
			{
				if(buffer[i] == -1)
				{
					buf_n[i]  = 	buf_n[i] + 1;
					buffer[i] = -1;
				}
				else
				{
					buf_n[i]  = 	buf_n[i] + 1;
					mu_bst    = 	get_muBest ( buf_n[i] , b );
					buffer[i] = 	get_R(buf_n[i] , mu + b , mu_bst + b );
				}

			}
		}
		max 		= 	buffer[0];
		n_highest 	= 	buf_n[0];

		for (int j = 0 ; j <7 ; j++) //find the highest R value from buffer and it's corresponding n
		{
			if(max <= buffer[j])
			{
				max 		= 	buffer[j];
				n_highest 	= 	buf_n[j];
				mub_highest =	get_muBest(n_highest , b);
				iter = j;
			}
		}

		if( n_highest <= n_min )
		{
			n_min = n_highest;
			buffer[iter] = -1;
		}
		else
		{
			n_max = n_highest;
			buffer[iter] = -1;
		}
		
		int n_index = iter; 		//the index is used to remember at which place in buffer the last R has been used
		P_Sum += poisson(n_highest , mu + b);


		while(P_Sum < CL)			// buffer moves up or down depending on whether the last used n was above or below n_top respectively
		{
			if( n_highest > n_top )
			{
				buf_n[n_index] 		= buf_n[n_index] +1 ;
				mu_bst    	 		= get_muBest ( buf_n[n_index] , b );
				buffer[n_index] 	= get_R(buf_n[n_index] , mu + b , mu_bst + b );

				for (int i = n_index + 1 ; i<7 ; i++)
				{
					if(buffer[i] == -1)
					{
						buf_n[i]  = 	buf_n[i] + 1;
						buffer[i] = 	-1;
					}
					else
					{
						buf_n[i]  = 	buf_n[i] + 1;
						mu_bst    = 	get_muBest ( buf_n[i] , b );
						buffer[i] = 	get_R(buf_n[i] , mu + b , mu_bst + b );
					}
				}
			}
			else if ( buf_n[0] != 0 )
			{
				buf_n[n_index] 	= buf_n[n_index] - 1;
				mu_bst    	 	= get_muBest ( buf_n[n_index] , b );
				buffer[n_index] = get_R(buf_n[n_index] , mu + b , mu_bst + b );


				for (int i = n_index - 1 ; i>=0 ; i--)
				{
					if(buf_n[i] - 1 <0 )
					{
						buf_n[i]  = 	 0;
						buffer[i] = 	-1; 
					}
					else
					{
						buf_n[i]  = 	buf_n[i] - 1;
						mu_bst    = 	get_muBest ( buf_n[i] , b );
						buffer[i] = 	get_R(buf_n[i] , mu + b , mu_bst + b );
					}
				}
			}
			
			max 		= 	buffer[0];
			n_highest 	= 	buf_n[0];

			for (int j = 0 ; j <7 ; j++) //find the highest R value from buffer and it's corresponding n
			{
				if(max <= buffer[j])
				{
					max 		= 	buffer[j];
					n_highest 	= 	buf_n[j];
					mub_highest =	get_muBest(n_highest , b);
					iter = j;
				}
			}

			if( n_highest <= n_min )
			{
				n_min = n_highest;
				buffer[iter] = -1;
			}
			else
			{
				n_max = n_highest;
				buffer[iter] = -1;
			}
			
			n_index = iter;
			P_Sum += poisson(n_highest , mu + b);

		}
		mu_U_temp[m] = mu;
		mu_L_temp[m] = mu;
		n_U_temp[m]  = n_min;
		n_L_temp[m]  = n_max;


		// cout<< "mu = " << mu_U_temp[m] << " \t nmin = " << n_U_temp[m] <<  "\t nmax = " << n_L_temp[m] << endl;


	}
	
	vector<int> n_U;
	vector<int> n_L; 

	vector<double> mu_U;
	vector<double> mu_L;


	int size = cycles;
	int n_current = n_U_temp[size - 1];
	// cout<< "n_current = "<< n_current<< endl;
	// int i = size -1;

	int y = 0;
	for (int x = size -1 ; x > 0; x--)  			//algortihm checks the n_array generated in get_n from back to front. 
													//Logs only jumps by one, so that there are no repetitions.
	{
		
		if (n_current > n_U_temp[x])
		{

			n_U.push_back(n_current);
		
			mu_U.push_back(mu_U_temp[x+1]);
			n_U.push_back(n_current - 1 );
			
		
			mu_U.push_back(mu_U_temp[x+1]);
		
			n_current = n_U_temp[x];
			y+=2;
			
		}
		
	}
	n_U.push_back(0);					//Adding point (0,0) so that the graph starts at 0. 
	mu_U.push_back(0);

	int mu_U_size = mu_U.size();
	reverse(mu_U.begin(),mu_U.end());
	reverse(n_U.begin(),n_U.end());

	for(int i = 0; i<mu_U_size; i++)
	{
		cout<< "n_U = " << n_U[i] << "\t mu_U = "<< mu_U[i] << endl;
	}

	delete[] mu_U_temp;
	delete[] mu_L_temp;
	delete[] n_U_temp;
	delete[] n_L_temp;

	return 0.0;
}

double MPFeldman_Cousins::get_R(int _n, double _lam1, double _lam2)
{
	double p_numer = poisson(_n, _lam1);
	double p_denom = poisson(_n, _lam2);
	
	if (p_denom != 0.0)
	{
		return poisson(_n, _lam1)/poisson(_n, _lam2);
	}
	else
	{
		cout << "ERROR: MPFeldman_Cousins::get_R(int, double, double): division by zero!" << endl;
		
		return -1.0; 
	}	
}

void MPFeldman_Cousins::fill_table()
{
	table = new double*[ROW_N];

    	for (int i = 0; i < ROW_N; i++) 
    	{
    	    table[i] = new double[COL_N];
    	}
    
    	for (int m = 0; m < COL_N; m++)
    	{
    	    for (int n = 0; n < ROW_N; n++)
    	    {
    	        double mu = double(m) * STEP;
	
		if (m == 0)
		{
			if (n == 0)
			{
				table[0][0] = 1.0;
			}
			else
			{
				table[n][0] = 0.0;
			}
		}	
		else
		{
			if (n == 0)
			{
				table[0][m] = exp(-mu);
			}
			else
			{
				table[n][m] = (mu / n) * table[n-1][m];
			}
		}	

		
		if ( table[n][m] < 1e-310 &&
		     table[n][m] > 0.0 )
		{
			// cout << "WARNING: MPFeldman_Cousins::fill_table(): value P(" << n << ", " << mu << ") lower than 1e-310, value was set to zero!" << endl;  
			table[n][m] = 0.0;
		}
            }
       	}

	table_set = true;

    	return;
}

void   MPFeldman_Cousins::print_table  (double _mumin, double _mumax)
{
	if (table_set)
	{
		cout << "*  mu|";
		
		for (int i = ceil(_mumin/STEP); double(i)*STEP < _mumax; i++)
		{
			cout << "        |";
		}

		cout << endl << " *   |";
		
		for (int i = ceil(_mumin/STEP); double(i)*STEP < _mumax; i++)
		{
			cout << fixed << setprecision(5) << setw(7) << double(i)*STEP << " |";
		}
	
		cout << endl << "n *  |";
		
		for (int i = ceil(_mumin/STEP); double(i)*STEP < _mumax; i++)
		{
			cout << "        |";
		}

		cout << endl << "-----";

		for (int i = ceil(_mumin/STEP); double(i)*STEP < _mumax; i++)
		{
			cout << "---------";
		}

		cout << endl;
		cout << showpoint;
		for (int n = 0; n < ROW_N; n++)
		{
			cout << setw(4) << n << " |";

			for (int m = floor(_mumin/STEP); double(m)*STEP < _mumax; m++)
			{
				cout << setw(6) << table[n][m] << " |";
			}

			cout << endl;
		}
	}
	else
	{
		cout << "ERROR: MPFeldman_Cousins::print_table(int, double): table with values was not set! Please call MPFeldman_Cousins::fill_table() to generate values first." << endl; 
	}
		
	return;
}

double MPFeldman_Cousins::poisson(int _n, double _mu)
{
	if (table_set)
	{
		int m1 = floor(_mu / STEP);
		int m2 = m1 + 1;

		double mu1 = double(m1) * STEP;
		double mu2 = double(m2) * STEP;
		double P1  = table[_n][m1];
		double P2  = table[_n][m2];
	
		double a   = (P1 - P2) / (mu1 - mu2);
		double b   = P1 - a * mu1; 
		
		return	a*_mu + b;	
	}
	else
	{
		cout << "ERROR: MPFeldman_Cousins::poisson(int, double): table with values was not set! Please call MPFeldman_Cousins::fill_table() to generate values first." << endl; 
		return NULL;
	}
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



























































/*


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

	mu_j = get_mu(); //???
	
	double* get_R_R = new double[ROW_N];
	int     n_start = 0;
	int      n_end  = ROW_N;

	for (int i = n_start; i<n_end;i++)
	{
		int mu_best = get_muBest(i, b);

		// cout<< " n_start is : "<< i << "   and int(mu+b)/STEP is: "<< int((mu_j+b)/STEP)<< "   and mub+b is: "<< int((mu_best+b))<< endl;
		// cout<< "table[i][mu+b] "<< table[i][int((mu_j+b)/STEP)]<< "   and table[i][mu_best]:   "<< table[i][int((mu_best+b)/STEP)]<< endl;

		int mu_i   = int( (mu_j   + b)/STEP );
		int mu_b_i = int( (mu_best+ b)/STEP );

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

TGraph* MPFeldman_Cousins::TGraph_calculate_lower()
{
	auto start = std::chrono::steady_clock::now();
	double mu_array[int(mu_max/STEP)];
	vector<int> n;
	vector<double> mu_L;
	std::vector<vector<int>> calculate_lower_n_array = get_n();

	for (int i = 0; i<mu_max/STEP; i++)
	{
		mu_array[i] = i*STEP;
		
	}


	int size = calculate_lower_n_array.size();
	int n_current = calculate_lower_n_array[size - 1][1];
	// cout<< "n_current = "<< n_current<< endl;
	// int i = size -1;

	int y = 0;
	for (int x = size -1 ; x > 0; x--)  			//algortihm checks the n_array generated in get_n from back to front. 
													//Logs only jumps by one, so that there are no repetitions.
	{
		
		if (n_current > calculate_lower_n_array[x][1])
		{

			if(calculate_lower_n_array[x][1] == 10000)
			{
				continue;

			}
			else
			{
				n.push_back(n_current);
			
				mu_L.push_back(mu_array[x+1]);
				n.push_back(n_current - 1 );
				
			
				mu_L.push_back(mu_array[x+1]);
			
				n_current = calculate_lower_n_array[x][1];
				y+=2;
			}
			
		}
		
	}

	// for (int i =0 ; i<y; i++)
	// {
	// 	cout<< "n = " << n[i] << "\t mu_L = " << mu_L[i] << endl;
	// 	// cout<< "Got to printing"<<endl;
	// }

	n.push_back(calculate_lower_n_array[0][1]);					//Adding first point so that the graph starts from bottom. 
	mu_L.push_back(mu_array[0]);

	int mu_L_size = mu_L.size();
	reverse(mu_L.begin(),mu_L.end());
	reverse(n.begin(),n.end());

	TGraph *gr  = new TGraph();


	for (int i =0; i<mu_L_size; i++)
	{
		gr->SetPoint(i, n[i], mu_L[i]);
	}
  	
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Time it took to calculate_lower: " << elapsed_seconds.count() << "s\n";
	// delete[] mu_array;
	return gr;
}

void MPFeldman_Cousins::draw_lower()
{

	TGraph* gr2 = new TGraph();
	gr2 = TGraph_calculate_lower();

	stringstream graph_title;
    graph_title  << "Lower Limit for bkg = " << b << " with C.L. = " << CL*100 << "%";
    string strname  = graph_title.str();

 	gr2->SetTitle(strname.c_str());
	gr2->SetFillStyle(1000);
	gr2->SetLineColor(kBlue);
	gr2->SetLineWidth(4);
	gr2->SetMarkerStyle(0);

	TAxis* xaxis;
	xaxis = gr2->GetXaxis();
	xaxis->CenterTitle(true);
	xaxis->SetTitle("n [Counts]");

	TAxis* yaxis;
	yaxis = gr2->GetYaxis();
	yaxis->SetTitle("\\mu_L  [Counts]");
	yaxis->CenterTitle(true);


	TCanvas *c2 = new TCanvas("c2","Graph Draw Options",600,600);


	c2->SetFrameFillStyle(3013);
	c2->SetFrameLineWidth(2);
	c2->SetFrameBorderMode(0);
	c2->SetFrameBorderSize(4);
	c2->SetGridx();
	c2->SetGridy();


	gr2->Draw("ALP");

	return;
}

vector<double> MPFeldman_Cousins::calculate_upper() ///CHANGE TO GET GRAPH!!!
{
	auto start = std::chrono::steady_clock::now();
	double* mu_array = new double[int(mu_max/STEP)];
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
	for (int x = size-1 ; x > 0; x--)
	{
		
		if (n_current > calculate_upper_n_array[x][0])
		{

			if(calculate_upper_n_array[x][0] == 10000)
			{
				// cout<<"watch out!" <<endl;
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
	// n.push_back(0);					//Adding point (0,0) so that the graph starts at 0. 
	// mu_U.push_back(0);

	// int entries = y;
	// cout<< "  n|  mu_U"<<endl;
	// cout<< "----------"<<endl;
	// for (int i =0 ; i<y+1; i++)
	// {
	// 	// cout<< "n = " << n[i] << "\t mu_U = " << mu_U[i] << endl;
	// 	cout<<setw(3)<< n[i]<<"|"<<setw(6) << mu_U[i] << endl;

	// 	// cout<< "Got to printing"<<endl;


	// }
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
    // std::cout << "Time it took to calculate_upper: " << elapsed_seconds.count() << "s\n";
	delete[] mu_array;

	return mu_U;
}

vector<double>  MPFeldman_Cousins::calculate_lower() 
{
	auto start = std::chrono::steady_clock::now();
	
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


  	
	auto end = std::chrono::steady_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;
    std::cout << "Time it took to calculate_upper: " << elapsed_seconds.count() << "s\n";
	// delete[] mu_array;

	return mu_L;
}



double* MPFeldman_Cousins::get_mu_U_v_b(int n)
{
	int range = 42;

	double* bkg      = new double[range];
	double* mu_U_v_b = new double[range];

	for (int x = 0; x < range; x++)
	{
		bkg[x] = b+x*0.5;
	}


	for (int i = 0; i<range; i++)
	{
		// cout<< "started loop get mu v b"<< endl;
		auto start = std::chrono::steady_clock::now();

		set_b(bkg[i]);

		// bkg[i] = get_b();

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
		// cout<< "b = " << bkg[i] << "\t mu_u = "<< mu_U_v_b[i]<<endl;

		// vector<double>().swap(mu_U);
		// mu_U.shrink_to_fit();
		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
	    // std::cout << "Time it took to find mu vs b: " << elapsed_seconds.count() << "s\n";

	}
	set_b(b_const);

	return mu_U_v_b;
} 

double* MPFeldman_Cousins::get_mu_L_v_b(int n)
{

	double* bkg      = new double[ROW_N];
	double* mu_L_v_b = new double[ROW_N];


	for (int i = 0; i<ROW_N; i++)
	{
		// cout<< "started loop get mu v b"<< endl;
		auto start = std::chrono::steady_clock::now();

		set_b(i*0.2); //checking b in 0.2 steps.

		bkg[i] = get_b();

		// cout<< "got to bkg = "<< bkg[i]<< endl;
		vector<double> mu_L;
		mu_L = calculate_lower();

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

		auto end = std::chrono::steady_clock::now();
		std::chrono::duration<double> elapsed_seconds = end-start;
	    std::cout << "Time it took to find mu vs b: " << elapsed_seconds.count() << "s\n";
	}
	return mu_L_v_b;
} 

void MPFeldman_Cousins::draw_mu_U_v_b(int n)
{
	int range = 42;

	double* mu_U_v_b = new double[range];
	double* bkg      = new double[range];
	double* mu_corrected = new double[range];


	TCanvas *c1 = new TCanvas("c1","Mu_u vs background",600,600);


	c1->SetFrameFillStyle(3013);
	c1->SetFrameLineWidth(2);
	c1->SetFrameBorderMode(0);
	c1->SetFrameBorderSize(4);
	c1->SetGridx();
	c1->SetGridy();


	mu_U_v_b = get_mu_U_v_b(n);
	for (int i = 0; i < range; ++i)
	{
		bkg[i] = i*0.2;
	}

	TGraph* gr = new TGraph(range,bkg,mu_U_v_b);

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
	gr->Draw("ALP");

	

	mu_corrected = correct_mu_U(n);

	TGraph* gr2 = new TGraph(range,bkg,mu_corrected);
	gr2->SetLineColor(kBlue);
	gr2->SetLineStyle(7);
	gr2->SetLineWidth(4);

	gr2->Draw("SAME");

	TLegend*  legend = new TLegend(0.6,0.6,0.9,0.7);
   	legend->AddEntry(gr,"Original \t \\mu_U","l");
   	legend->AddEntry(gr2,"Corrected \t \\mu_U","l");
   	legend->Draw();


	// gr->Draw("ALP");
	return;
}

void MPFeldman_Cousins::draw_mu_L_v_b(int n)
{
	double* mu_L_v_b = new double[ROW_N];
	double* bkg      = new double[ROW_N];
	double* mu_corrected = new double[ROW_N];


	TCanvas *c1 = new TCanvas("c1","Mu_L vs background",600,600);


	c1->SetFrameFillStyle(3013);
	c1->SetFrameLineWidth(2);
	c1->SetFrameBorderMode(0);
	c1->SetFrameBorderSize(4);
	c1->SetGridx();
	c1->SetGridy();


	mu_L_v_b = get_mu_L_v_b(n);
	for (int i = 0; i < ROW_N; ++i)
	{
		bkg[i] = i*0.2;
	}

	TGraph* gr = new TGraph(ROW_N,bkg,mu_L_v_b);

	stringstream graph_title;
    graph_title  << "Lower Limit for n = " << n << " with C.L. = " << CL*100 << "%";
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
	yaxis->SetTitle("\\mu_L  [Counts]");
	yaxis->CenterTitle(true);
	gr->Draw("ALP");

	

	mu_corrected = correct_mu_L(n);

	TGraph* gr2 = new TGraph(ROW_N,bkg,mu_corrected);
	gr2->SetLineColor(kBlue);
	gr2->SetLineStyle(7);
	gr2->SetLineWidth(4);

	gr2->Draw("SAME");

	TLegend*  legend = new TLegend(0.6,0.6,0.9,0.7);
   	legend->AddEntry(gr,"Original \t \\mu_U","l");
   	legend->AddEntry(gr2,"Corrected \t \\mu_U","l");
   	legend->Draw();


	// gr->Draw("ALP");
	return;
}

double* MPFeldman_Cousins::correct_mu_U(int n)
{
	int range = 42;

	double* mu_U_corrected = new double[range];
	mu_U_corrected = get_mu_U_v_b(n);
	for (int i = range-1; i>=0; i--)
	{
		if(i == 0)
		{	
			continue;
		}

		else if(mu_U_corrected[i]<=mu_U_corrected[i-1])
		{
			mu_U_corrected[i-1] = mu_U_corrected[i-1];
		}
		else if(mu_U_corrected[i]>mu_U_corrected[i-1])
		{
			mu_U_corrected[i-1] = mu_U_corrected[i];
		}
	}
	// for (int x = 0; x<range;x++)
	// {
	// 	cout<< mu_U_corrected[x]<<endl;
	// }
	return mu_U_corrected;
}

double* MPFeldman_Cousins::correct_mu_L(int n)
{
	double* mu_L_corrected = new double[ROW_N];
	mu_L_corrected = get_mu_L_v_b(n);
	for (int i = ROW_N-1; i>=0; i--)
	{
		if(i == 0)
		{	
			continue;
		}

		else if(mu_L_corrected[i]<=mu_L_corrected[i-1])
		{
			mu_L_corrected[i-1] = mu_L_corrected[i-1];
		}
		else if(mu_L_corrected[i]>mu_L_corrected[i-1])
		{
			mu_L_corrected[i-1] = mu_L_corrected[i];
		}
	}
	// for (int x = 0; x<ROW_N;x++)
	// {
	// 	cout<< mu_L_corrected[x]<<endl;
	// }
	return mu_L_corrected;
}

double* MPFeldman_Cousins::mu_U_final()
{
	int range = 42;
	double* mu_U_final = new double[range];

	cout<< "  n|  mu_U"<<endl;
	cout<< "----------"<<endl;
	for (int i = 0; i< range/2; i++)
	{
		double* mu_U_corrected = new double[range];
		mu_U_corrected = correct_mu_U(i);
		mu_U_final[i] = mu_U_corrected[0];

		cout<<setw(3)<< i<<"|"<<setw(6) << mu_U_final[i] << endl;
		delete[] mu_U_corrected;
	}

	return mu_U_final;
}*/
