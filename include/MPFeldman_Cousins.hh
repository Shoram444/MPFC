#include <vector>
// ROOT headers
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "TObject.h"

class MPFeldman_Cousins: public TObject 
{ 
	public:    
    	MPFeldman_Cousins();
    	MPFeldman_Cousins(double _b, double _step, int _rows, double _mu_max, double _CL);
    	~MPFeldman_Cousins();
        double** fill_table();
	    double poisson(int n, double mu_j);  					// Generates a table of ROW_N rows and COL_N columns and fills it with Poisson values with n = ROW_N[i] and mu = COL_N[i]*STEP
	    double* get_R();   										// Generates an array of R ratios where R = P(n,mu+b)/P(n, mu_best+b)
	    double get_muBest(int n, double b); 					// Returns mu_best
	    int* get_A(double* R);									// Generates an array of indexes, which point to R values from highest to lowest. (index = 0 is the highest R , index = 1 is second Highest R and so on.)
	    void print_poisson();									// Prints out the table of P values
	    void print_R();											// Prints array of R values
	    void set_mu(double _mu_j);								// Sets mu for priting R array
	    double get_mu();										// Outputs mu
	    void print_A();											// Prints array of n and their according ordering number
	    std::vector<int> CL_check(double mu_j, double* _R);					// Checks for the CL condition and outputs array of indexes for n
	    std::vector<std::vector<int> > get_n();					// Outputs array of n_min and n_max associated with each mu step
        std::vector<double> calculate_upper();									// Calculates TGraph for upper limit
        TGraph* TGraph_calculate_upper();
        std::vector<double>  calculate_lower();									// Calculates TGraph for upper limit
        void draw_upper();										// draws upper limit
        void print_n();                                         // 
        void set_b(double _b);
        double get_b();
        double* get_mu_U_v_b(int n);
        void get_mu_L_v_b(int n);
        void draw_mu_U_v_b(int n);


    private:

        double b;
        double STEP;  											//step at which mu iterates. 
		int ROW_N;												//Number of n.
		int COL_N;  											//Number of columns
		double mu_max;
        double** table;											//2D Array generated in constructor - holds values of P. Row -> n ; Column -> mu/STEP. 
        double* R;												//1D Array generated by get_R, provides array of Ratios
        double mu_j;											//mu
        int* A;													//1D Array generated by get_A, provides array of indexes from R in order from highest to lowest. (e.g. highest ratio is R[A[0]]. )
        std::vector<int> n_order;								//Provides n_min and n_max 
        double CL; 												//Desired confidence level. 
        std::vector<std::vector<int> > n_array;  				// holds table of n_min and n_max in columns and mu_j/step in rows
        // TGraph* calculate_upper();									// Calculates TGraph for upper limit
	    // TGraph* calculate_lower();									// Calculates TGraph for lower limit

        ClassDef(MPFeldman_Cousins,1);
};