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
        MPFeldman_Cousins         ();
        MPFeldman_Cousins         (double _b, double _CL);
        ~MPFeldman_Cousins        ();

        // Mathematical methods
        double  get_muBest        (int n, double b, bool _warn = true); 
        double  get_R             (int _n, double _lam1, double _lam2, bool _warn = true); 
        double  poisson           (int _n, double _mu, bool _warn = true);
        double  get_mu_U          (int _n, double _b);
        double  get_sensitivity   ();

        // Printing and visualization
        void    draw_upper        ();
        void    print_table       (double _colmin, double _colmax, int _nmin, int _nmax, bool _isp = false);
        void    print_symb        (string _s, int _rep);

        // Extend method
        void    m_extend          ();
        // Static variables to save table of poisson values and limits for mu
        static double** m_table;            
        static bool     m_table_set;

        static double** p_table;              
        static bool     p_table_set;

    private:

        // Background and confidence level values given by user
        double       b;             
        double       CL;        

        // Constants given by the design of the algorithm
        const double STEP_P  = 0.01;      
        const int    ROWS_P  = 1000;             
        const int    COLS_P  = ROWS_P/STEP_P;    

        const double STEP_M = 0.1;    
        const int    NLST_M = 10; 
              int    COLS_M; 
        const double MUPREC = 0.01;              

        // Calculation algorithm methods
        belt    calculate_limit       (double _mu, double _b);
        void    fill_m_table          ();
        void    fill_p_table          ();  
        double* shift_mu_U            (std::vector<belt> _bt);
                            
        ClassDef(MPFeldman_Cousins,1);
};
