// Standard headers
#include <iostream>
#include <math.h> 
#include <vector>
#include <bits/stdc++.h> 
#include <chrono>
#include <algorithm>

// ROOT headers
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TLegend.h"

// MPFeldman_Cousins header
#include "MPFeldman_Cousins.hh"

ClassImp(MPFeldman_Cousins);

double** MPFeldman_Cousins::p_table       =   NULL;
bool     MPFeldman_Cousins::p_table_set   =   false;

double** MPFeldman_Cousins::m_table       =   NULL;
bool     MPFeldman_Cousins::m_table_set   =   false;

using namespace std;
//***************************************************************
// PUBLIC METHODS
//***************************************************************
//////////////////////////////////////////////////////
// Constructors and destructor

MPFeldman_Cousins::MPFeldman_Cousins()
{
}

MPFeldman_Cousins::MPFeldman_Cousins(double _b, double _CL)
{
    b         = _b;
    CL        = _CL;

	COLS_M = 3*_b/STEP_M;

    if(!p_table_set)
    {
        auto start = chrono::steady_clock::now();
        fill_p_table();
        auto end   = chrono::steady_clock::now();
        chrono::duration<double> elapsed_seconds = end-start;
        
        cout << "Time it took to fill table: " << elapsed_seconds.count() << "s\n";
    }
    
    if(!m_table_set)
    {
    	auto start = chrono::steady_clock::now();
    	fill_m_table();
    	auto end   = chrono::steady_clock::now();
    	chrono::duration<double> elapsed_seconds = end-start;

    	cout << "Time it took to fill table: " << elapsed_seconds.count() << "s\n";
    }
}

MPFeldman_Cousins::~MPFeldman_Cousins()
{
}

//////////////////////////////////////////////////////
// Mathematical methods

double MPFeldman_Cousins::get_muBest(int n, double b, bool _warn)
{

    double mu_bst = double(n) - b;

    if(n < 0 || b < 0)
    {
        if (_warn)
        {
            cout << "WARNING: MPFeldman_Cousins::get_muBest(int n, double b, bool): negative values given as argument! Please check input values!" << endl; 
        }

        mu_bst = -1.0;
    }
    else if (n == 0 && b == 0)
    {
        mu_bst = 0;
    }
    else if(mu_bst < 0)
    {
        mu_bst = 0;
    }

    return mu_bst;
}

double MPFeldman_Cousins::get_R(int _n, double _lam1, double _lam2, bool _warn)
{
    double p_numer = poisson(_n, _lam1, _warn);
    double p_denom = poisson(_n, _lam2, _warn);
    
    if (_n < 0 || _lam1 < 0.0 || _lam2 < 0.0)
    {
        if (_warn)
        {
            cout << "WARNING: MPFeldman_Cousins::get_R(int, double, double, bool): one of the arguments is negative!" << endl;
        }       

        return -1.0; 
    }
    else if (p_denom != 0.0)
    {
        return poisson(_n, _lam1, _warn)/poisson(_n, _lam2, _warn);
    }
    else
    {
        cout << "ERROR: MPFeldman_Cousins::get_R(int, double, double, bool): division by zero!" << endl;
        
        return -1.0; 
    }   
}

double MPFeldman_Cousins::poisson(int _n, double _mu, bool _warn)
{
    if (_n < 0 || _mu < 0)
    {
        if(_warn)
        {
            cout << "WARNING: MPFeldman_Cousins::poisson(int, double, bool): negative values given as argument! Please check input values!" << endl; 
        }
        
        return -1.0;
    }
    
    if (p_table_set)
    {
        int m1 = floor(_mu / STEP_P);
        int m2 = m1 + 1;

        double mu1 = double(m1) * STEP_P;
        double mu2 = double(m2) * STEP_P;
        double P1  = p_table[_n][m1];
        double P2  = p_table[_n][m2];
    
        double a   = (P1 - P2) / (mu1 - mu2);
        double b   = P1 - a * mu1; 
        
        return  a*_mu + b;  
    }
    else
    {
        cout << "ERROR: MPFeldman_Cousins::poisson(int, double, bool): table with values was not set! Please call MPFeldman_Cousins::fill_p_table() to generate values first." << endl; 
        return NULL;
    }
}

//////////////////////////////////////////////////////
// Printing and visualization

void MPFeldman_Cousins::draw_upper()
{
    vector<belt> bt;                  

    int m    =    0;
    int _n_0 =   10;

    bt.push_back( calculate_limit( m, b ) );

    m = 1;

    while( bt.back().n_min < _n_0 + 1 )

    {
        bt.push_back( calculate_limit( double( m * STEP_M ) , b ) );
        m++;
    }

    double* v  = shift_mu_U( bt);
    cout<< "v fine " << endl;

    int*           n =       new int[int(2*_n_0)];
    double*       mu =   new double[int(2*_n_0)];
    int            y =                         0;

    TGraph *gr  = new TGraph();

    
    for (int i = 0; i<2*_n_0 - 1; i+=2)
    {
        n[i]    = y;
        n[i+1]  = y;
        y++;
    }

    cout<< " n fine " << endl;

    mu[0] = 0 ;

    y     = 0 ;

    for (int i = 1; i<2*_n_0 -1 ; i+=2)
    {
        mu[i]   = v[y];
        mu[i+1] = v[y];
        cout<< "mu[i] = " << mu[i] << "\t mu[i+1] = " <<mu[i+1] <<endl;
        y++;
    }

    cout<< " mu fine " << endl;


    for (int i =0; i<2*_n_0 - 1; i++)
    {
        gr->SetPoint(i, n[i], mu[i]);
    }

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

void   MPFeldman_Cousins::print_table  (double _colmin, double _colmax, int _nmin, int _nmax, bool _isp)
{
	// _isp = True - p_table, _isp = False -> m_table
	double** tbl;

	if (_isp)
	{
        tbl = p_table;
        cout << "*******************************************" << endl;
        cout << "* Table of chosen values of Poisson(n,mu) *" << endl;
        cout << "*******************************************" << endl;
	}
	else
	{
    	tbl = m_table;
        cout << "*******************************************" << endl;
        cout << "*  Table of chosen values of Mu_max(n,b)  *" << endl;
        cout << "*******************************************" << endl;
	}

	double stp = _isp ? STEP_P : STEP_M;

	if ( ( _isp && p_table_set ) ||
    	 (!_isp && m_table_set ) )
    {
    	int min_order = 0;
    	int max_order = 0;

    	for (int n = _nmin; n < _nmax; n++)
        {
        	for (int m = floor(_colmin/stp); double(m)*stp <= _colmax; m++)
            {
            	if( ( min_order < ceil(-log10(_isp ? tbl[n][m] : tbl[m][n])) ) && 
            		(_isp ? tbl[n][m] : tbl[m][n] > 1e-10 ) )
            	{
            		min_order = ceil(-log10(_isp ? tbl[n][m] : tbl[m][n]));
            	}

            	if( max_order < ceil( log10(_isp ? tbl[n][m] : tbl[m][n])) )
            	{
            		max_order = ceil( log10(_isp ? tbl[n][m] : tbl[m][n]));
            	} 
            }
        }

        if (min_order > 9)
        {
        	min_order = 9;
        }

        if (min_order < 2)
        {
        	min_order = 2;
        }

        print_symb("-", 5);

        for (int i = ceil(_colmin/stp); double(i)*stp <= _colmax; i++)
        {
        	print_symb("-", max_order + min_order + 4);
        }

        cout << endl;

    	if (_isp)
    	{
        	cout << "*  mu|";
        }
        else
        {
        	cout << "*  b |";	
        }
        
        for (int i = ceil(_colmin/stp); double(i)*stp <= _colmax; i++)
        {
        	print_symb(" ", max_order + min_order + 3);
            cout << "|";
        }

        cout << endl << " *   |";
        
        for (int i = ceil(_colmin/stp); double(i)*stp <= _colmax; i++)
        {
            cout << fixed << setprecision(min_order) << setw(max_order + min_order + 2) << double(i)*stp << " |";
        }
    
        cout << endl << "n *  |";
        
        for (int i = ceil(_colmin/stp); double(i)*stp <= _colmax; i++)
        {
        	print_symb(" ", max_order + min_order + 3);
            cout << "|";
        }

        cout << endl;
        print_symb("-", 5);

        for (int i = ceil(_colmin/stp); double(i)*stp <= _colmax; i++)
        {
        	print_symb("-", max_order + min_order + 4);
        }

        cout << endl;
        cout << showpoint;

        for (int n = _nmin; n < _nmax; n++)
        {
            cout << setw(4) << n << " |";

            for (int m = floor(_colmin/stp); double(m)*stp <= _colmax; m++)
            {
                	cout << setprecision(min_order) << setw(max_order + min_order + 2) << (_isp ? tbl[n][m] : tbl[m][n])  << " |";
            }

            cout << endl;
        }

        print_symb("-", 5);

        for (int i = ceil(_colmin/stp); double(i)*stp <= _colmax; i++)
        {
        	print_symb("-", max_order + min_order + 4);
        }
        cout << endl;
    }
    else
    {
        cout << "ERROR: MPFeldman_Cousins::print_table(double, double, int, int, bool): table with values was not set! Please call generate the values first." << endl; 
    }
}

void MPFeldman_Cousins::print_symb(string _s, int _rep)
{
	for (int i = 0; i < _rep; i++)
        	{
            	cout << _s;
            }
}

//***************************************************************
// PRIVATE METHODS
//***************************************************************
//////////////////////////////////////////////////////
// Calculation algorithm methods

belt MPFeldman_Cousins::calculate_limit(double _mu, double _b)
{
    belt*     r = new belt;
    
    double  bkg =        _b;
    double mu   =      _mu;
    int n_min   =        0;
    int n_max   =        0;

    double mu_bst;

    double* buf_R  = new double[7]; //Prepare buffer for R values of size 7
    double* buf_n  = new double[7]; //Prepare buffer of corresponding n values associated to buffer
    
    // cout<< "mu = " << mu << endl << " n_top : "<< ceil( mu + bkg ) << endl << " First round of buffering" <<endl;

    for(int i = 0 ; i < 7; i++) //Filling out buffer around the n_top value. 
    {
        if (ceil( mu + bkg ) - 3 + i < 0)  //Makes sure that the buffer of n values doesn't go to negatives 
        {
            buf_n[i]  =      0;
            buf_R[i]  =     -1; 
        }
        else
        {
            buf_n[i]  =     ceil( mu + bkg ) - 3 + i;
            buf_R[i]  =     get_R(buf_n[i], mu + bkg, get_muBest(buf_n[i], bkg, false) + bkg, false);
        }

    }

    double  P_Sum       = 0.0;

    n_min               = ceil( mu + bkg ); //this initiates the n_min to the "n_top" value (middle of the buffer)
    n_max               = ceil( mu + bkg ); //this initiates the n_max to the "n_top" value (middle of the buffer)


    while(P_Sum < CL)
    {   
        int    cho_n = 0;   //cho_n is the corresponding n for highest R
        double cho_R = 0.0; //cho_R is used to search for highest R value
        double cho_P = 0;   //cho_n is the corresponding n for highest R
        int    cho_j = 0;   //cho_n is the corresponding n for highest R

        // cout << "PO| n |   P    |    P_bst    |     R   " << endl;
        // for(int i = 0 ; i < 7; i++)
        //     {
        //         cout << i << " | " << buf_n[i] << " | " << poisson(buf_n[i] , mu + bkg, false) << " | " << poisson(buf_n[i] , get_muBest(buf_n[i], bkg, false) + bkg, false) << " | " << buf_R[i] << endl;
        //     }
        // cout << endl << endl;

        for (int j = 0 ; j < 7 ; j++)       //find the highest R value from buffer and it's corresponding n
        {
            if(cho_R <= buf_R[j])
            {
                cho_n       =   buf_n[j];
                cho_R       =   buf_R[j];
                cho_P       =   poisson(cho_n , mu + bkg, false);
                cho_j       =   j;      // cho_j is used to make sure that the R value that has been used once will not be used again
            }

        }


        if(n_min >= cho_n)
        {
            n_min       =   cho_n; 
        }
        else
        {
            n_max       =   cho_n; 
        }
                         
        P_Sum       +=  cho_P;

        // cout << "*****CHOSEN = " << cho_n << " cho_j " << cho_j << ", P(cho_n) = " << poisson(cho_n , mu + bkg, false) << "\t P_Sum = " << P_Sum << endl << endl;

        if( cho_n == 0 || 
                get_R(buf_n[0] - 1, mu + bkg, get_muBest(buf_n[0] - 1, bkg, false) + bkg, false ) 
              > get_R(buf_n[6] + 1, mu + bkg, get_muBest(buf_n[6] + 1, bkg, false) + bkg, false ) ) //this condition makes sure that when moving buffer toward lower numbers, we do not get out of bounds.
        {                   
            if (cho_j == 0)
            {
                buf_n[0]    =  buf_n[0] - 1;
                buf_R[0]    =  get_R(buf_n[0], mu + bkg, get_muBest(buf_n[0], bkg, false) + bkg, false );                   
            }
            else 
            {
                while (cho_j >= 1)
                {
                    buf_R[cho_j]    =  buf_R[cho_j - 1];
                    buf_n[cho_j]    =  buf_n[cho_j - 1];

                    cho_j--;
                }

                buf_n[0]    =  buf_n[1] - 1;
                buf_R[0]    =  get_R(buf_n[0], mu + bkg, get_muBest(buf_n[0], bkg, false) + bkg, false );   
            }                   
        }
        else 
        {
            if (cho_j == 6)
            {
                buf_n[6]    =  buf_n[6] + 1;
                buf_R[6]    =  get_R(buf_n[6], mu + bkg, get_muBest(buf_n[6], bkg, false) + bkg, false );   
            }
            else
            {
                while (cho_j <= 5)
                {
                    buf_R[cho_j]    =  buf_R[cho_j + 1];
                    buf_n[cho_j]    =  buf_n[cho_j + 1];

                    cho_j++;
                }

                buf_n[6]    =  buf_n[5] + 1;
                buf_R[6]    =  get_R(buf_n[6], mu + bkg, get_muBest(buf_n[6], bkg, false) + bkg, false );   
            }
        }    
// 
        // cout<< " Cho_n =  "<< cho_n << " \t P_0 = " << cho_P << "\t P_Sum = " << P_Sum << endl;  
    }

    r->mu        =     mu;
    r->n_min     =  n_min;  
    r->n_max     =  n_max;

    return *r;
}

void MPFeldman_Cousins::fill_m_table()
{   
	vector<belt> mu;          
	m_table = new double*[COLS_M];

    for (int bkg = 0; bkg <= COLS_M; bkg++)
	{
    	double m = 0.0;
		//cout << "bkg " << bkg << endl;

        do
        {
            mu.push_back( calculate_limit(m, bkg*STEP_M ) );
            m += MUPREC;
        }
        while( mu.back().n_min <= NLST_M + 1 );

        m_table[bkg] = shift_mu_U(mu);

        mu.clear();
	}
    
    m_table_set = true;
    

    for(int n = 0; n <= NLST_M; n++)
    {
        for(int bkg = COLS_M - 1; bkg > 0; bkg--)
        {
            if( m_table[bkg-1][n] <= m_table[bkg][n])
            {
                m_table[bkg-1][n] = m_table[bkg][n];
            }
        }
    }


    return;
}

void MPFeldman_Cousins::fill_p_table()
{
    p_table = new double*[ROWS_P];

        for (int i = 0; i < ROWS_P; i++) 
        {
            p_table[i] = new double[COLS_P];
        }
    
        for (int m = 0; m < COLS_P; m++)
        {
            for (int n = 0; n < ROWS_P; n++)
            {
                double mu = double(m) * STEP_P;
    
                if (m == 0)
                {
                    if (n == 0)
                    {
                        p_table[0][0] = 1.0;
                    }
                    else
                    {
                        p_table[n][0] = 0.0;
                    }
                }   
                else
                {
                    if (n == 0)
                    {
                        p_table[0][m] = exp(-mu);
                    }
                    else
                    {
                        p_table[n][m] = (mu / n) * p_table[n-1][m];
                    }
                }   

                if ( p_table[n][m] < 1e-310 &&
                     p_table[n][m] > 0.0 )
                {
                    p_table[n][m] = 0.0;
                }
            }
        }

    p_table_set = true;

    return;
}

double* MPFeldman_Cousins::shift_mu_U(vector<belt> _bt)
{
	int n_last = _bt.back().n_min;
	double* mu = new double[n_last];
	
	for (int i = 0; i < n_last; i++)
    {
        mu[i] = -1.0;
    }

    int j = _bt.size() - 1;
    int n = _bt[j].n_min;
    do
    {   
        if ( mu[_bt[j-1].n_min] == -1           &&      
        	    _bt[j-1].n_min   < _bt[j].n_min )
        {
            mu[_bt[j-1].n_min] = _bt[j].mu;
        }
        j--;
    }
    while (j > 0);              //algortihm checks the n_array generated in get_n from back to front. 
                                                          //Logs only jumps by one, so that there are no repetitions.
    return mu;
}