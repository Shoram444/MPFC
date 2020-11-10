// Standard C++ headers
#include <vector>

// ROOT headers
#include "TROOT.h"  
#include "TGraph.h"
#include "TCanvas.h"
#include "TObject.h"

using namespace std;

struct belt
{
    double  mu;
    int     n_min;
    int     n_max;
};



class MPFeldman_Cousins: public TObject 
{ 
    public:    
        // Constructors and destructor
        MPFeldman_Cousins             ();
            MPFeldman_Cousins         (double _b, double _CL);
            ~MPFeldman_Cousins        ();

        double  calculate_lim         ();
        double  get_R                 (int _n, double _lam1, double _lam2, bool _warn = true);   // Generates an array of R ratios where R = P(n,mu+b)/P(n, mu_best+b)
        void    fill_table            ();
        void    print_table           (double _mumin, double _mumax);
        double  print_poisson         (int _n, double _mu);         
        double  get_muBest            (int n, double b, bool _warn = true);
        double  poisson               (int _n, double _mu, bool _warn = true);
        void    fill_m_table          (int _n_0);
        void    draw_upper            ();


        belt                calculate_limit       (double _mu, double _b);
        double*              shift_mu_U             (std::vector<belt> _bt , int _n_0);
        

        static double** p_table;             // 2D Array generated in constructor - holds values of P. Row -> n ; Column -> mu/STEP. 
        static bool     p_table_set;

        static double** m_table;             // 2D Array generated in constructor - holds values of P. Row -> n ; Column -> mu/STEP. 
        static bool     m_table_set;

        // vector<belt>  b;













/*



        
    
        
        // "Calculate" methods
            vector<double>  calculate_lower();                  // Calculates TGraph for upper limit
            vector<double>  calculate_upper();          // Calculates TGraph for upper limit
            TGraph*         TGraph_calculate_upper();               // Same as calculate upper, but returns TGraph for visualisation
            TGraph*         TGraph_calculate_lower();

        // "Draw" methods
            void draw_upper   ();                                   // draws upper limit
            void draw_lower   ();                       // draws upper limit
            void draw_mu_U_v_b(int n);                              // draws graph of mu vs b, with original mu from get_mu_v_b and correct mu
            void draw_mu_L_v_b(int n);                              // draws graph of mu vs b, with original mu from get_mu_v_b and correct mu
    
    // "Get" methods
        
            double           get_muBest  (int n, double b);     // Returns mu_best
            int*             get_A       (double* R);       // Generates an array of indexes, which point to R values from highest to lowest. (index = 0 is the highest R , index = 1 is second Highest R and so on.)
            double           get_mu      ();            // Outputs mu
            double           get_b       ();                    // gets b
            double*              get_mu_U_v_b(int n);               // outputs array of mu values based on different backgrounds
            double*              get_mu_L_v_b(int n);               // outputs array of mu values based on different backgrounds
            vector<vector<int> > get_n       ();            // Outputs array of n_min and n_max associated with each mu step

    // "Print" methods
        void print_poisson();                   // Prints out the table of P values
            void print_R();                     // Prints array of R values
        void print_A();                     // Prints array of n and their according ordering number
            void print_n();                                     // 

        // "Set" methods
            void set_b (double _b);                                 // sets b
        void set_mu(double _mu_j);              // Sets mu for priting R array

        // Other methods
        vector<int> CL_check(double mu_j, double* _R);      // Checks for the CL condition and outputs array of indexes for n
            double*     correct_mu_U(int n);                        // outputs array of corrected mu values based on get mu -> correction is done so that higher b doesnt give lower mu
            double*     correct_mu_L(int n);
            double*     mu_U_final();                               // outputs array of mu for different n and constant b. The array is corrected based on correct_mu()
*/      

    private:

        const double STEP   = 0.01;         // step at which mu iterates. 
        const int    ROW_N  = 1000;     // Number of n.
        const int    COL_N  = ROW_N/STEP;            // The number of columns dependent on max mu and step chosen.
        double       b;                         // background b, variable through get and set
        double       CL;            // Desired confidence level. 
        
        vector<belt> mu_U;
                    
                     // 2D Array generated in constructor - holds values of P. Row -> n ; Column -> mu/STEP. 

        // bool     table_set;















        // double b_const;                                         // constant b, keeps same as in constructor
        
            
        // double mu_j;                        //mu
        // int* A;                         //1D Array generated by get_A, provides array of indexes from R in order from highest to lowest. (e.g. highest ratio is R[A[0]]. )
        // vector<int> n_order;                    //Provides n_min and n_max 
        // vector<vector<int> > n_array;               // holds table of n_min and n_max in columns and mu_j/step in rows
        // TGraph* calculate_upper();               // Calculates TGraph for upper limit
        // TGraph* calculate_lower();               // Calculates TGraph for lower limit

        ClassDef(MPFeldman_Cousins,1);
};
