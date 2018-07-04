// Copyright Alessandro Fabbri 2015

#include <tba.hpp>
#include <utils.hpp>

/*
 *  Lee-Yang TBA class with full features
 */
enum
{
  TBA_LEEYANG,
  TBA_SINEGORDON
};

template<int TBA_MODEL>
class tba_single : public tba
{
public:
  // Functions
  double *eps;                               // pseudoenergy
  double *L;                                 // log(1 + e^-epsilon)
  double *phi;                               // kernel
  double *nu;                                // forcing term

  // Constructors
  tba_single() {}                            // default constructor

  ~tba_single() {}                           // default destructor

  tba_single(int Nstep_)                     // step size constructor
  {
    Nstep = Nstep_;
    x     = new double [Nstep];
    eps   = new double [Nstep];
    L     = new double [Nstep];
    nu    = new double [Nstep];
    phi   = new double [2*Nstep-1];
  }

  void operator=(const tba_single<TBA_MODEL> &tba_ly)
  {
    Nstep = tba_ly.Nstep;
    iter  = tba_ly.iter;
    error = tba_ly.error;
    x     = tba_ly.x;
    eps   = tba_ly.eps;
    L     = tba_ly.L;
    nu    = tba_ly.nu;
    phi   = tba_ly.phi;
  }

  tba_single(const tba_single &tba_ly)
  {
    Nstep = tba_ly.Nstep;
    iter  = tba_ly.iter;
    error = tba_ly.error;
    x     = tba_ly.x;
    eps   = tba_ly.eps;
    L     = tba_ly.L;
    nu    = tba_ly.nu;
    phi   = tba_ly.phi;
  }

  // Methods
  void init_all(double xmax, double r_)      // initial operations, divided into
  {
    iter = 0;
    error = 1;
    r = r_;
    c = 0;
    x_max= ( r < 1e-2 ? 2.*abs(log(2/r)) : xmax );
    x_min = -x_max;
    dx = (x_max-x_min)/(Nstep-1.);
    for(int i=0; i<Nstep; i++) x[i] = x_min+i*dx;

    setup << "Lee-Yang TBA setup" << std::endl
          << "Nstep   = " << Nstep << std::endl
          << "X Range = [ "  << x_min << " , " << x_max << " ] " << std::endl
          << "dx      = " << dx << std::endl
          << "r       = " << r << std::endl;

    init_kernels();
    init_guess();
  }

  void init_guess()                          // initial guess for epsilon
  {
    for(int i=0; i<Nstep; i++)
    {
      nu[i] = r*cosh(x[i]);
      eps[i] = nu[i];
      L[i] = log(1+exp(-eps[i]));
    }
  }

  void init_kernels()                        // convolution kernel init
  {
    double xi;
    for(int i=0; i<2*Nstep-1; i++)
    {
      xi = 2.*x_min + i*dx;
      if( abs(xi) < 1E-25 ) xi = 1E-10;
      phi[i] = sqrt(3.)/M_PI*sinh(2.*xi)/sinh(3.*xi);
    }
  }

  void solve()                               // solution loop
  {
    double eps_old_0=eps[Nstep/2];
    double temp_sum=0;
    do{
      for(int i=0; i<Nstep; i++)
      {
        temp_sum=0;
        for(int j=0; j<Nstep; j++)
        {
          int ij = Nstep-1+abs(i-j);
          temp_sum += phi[ij]*L[j];
        }
        eps[i] = nu[i] + dx*temp_sum;
        L[i] = log(1+exp(-eps[i]));
      }
      error = abs(eps[Nstep/2]-eps_old_0);
      eps_old_0 = eps[Nstep/2];
      iter++;
    } while(error > PRECISION);
  }

  void evaluate_c()                          // evaluation of the c-function
  {
    for(int i=0;i<Nstep;i++)
      c += cosh(x[i])*L[i];
    c*=dx*6/M_PI/M_PI*r;

    report << "Lee-Yang TBA report" << std::endl
           << "Iter  = " << iter << std::endl
           << "Error = " << error << std::endl
           << "c     = " << c << std::endl;
  }

  void save_kernels(std::string filename)
  {
    std::ofstream out(filename);
    if(!out)
    {
      std::cerr << "Error creating kernel file : " << filename << std::endl;
      exit(2);
    }
    double xi;
    out << "index;x;phi" << std::endl;
    for(int i=0; i<2*Nstep-1; i++)
    {
      xi = 2.*x_min + i*dx;
      out << i << ";" << xi << ";" << phi[i] << std::endl;
    }

    out.close();
  }

  void save_results(std::string filename)    // save results to file
  {
    std::ofstream out(filename);
    out << "index;x;eps;L" << std::endl;
    for(int i=0; i<Nstep; i++) out << i << ";" << x[i] << ";" << eps[i] << ";" << L[i] << std::endl;
    out.close();
  }
};

/*
 *  C-Function Lee-Yang class for automated evaluation of c(r)
 */
template<int TBA_MODEL>
class cfunc_single{
public:
  // Variables
  int Nr;                                                      // number of r values
  double *r;                                                   // r values
  double *c;                                                   // c(r) function
  tba_single<TBA_MODEL> *tba_ly;                               // tba, one for each r
  // Constructors
  cfunc_single() {};                                           // default constructor

  ~cfunc_single() {};                                          // default destructor

  cfunc_single(int Ns, double r_min, double r_max, int Nr_)    // constructor range and steps
  {
    Nr = Nr_;
    r = new double[Nr];
    c = new double[Nr];
    tba_ly = new tba_single<TBA_MODEL>[Nr];
    double dr = (r_max-r_min)/(Nr-1);
    for(int i=0; i<Nr-1; i++) r[i] = r_min + i*dr;
    r[Nr-1] = r_max;
    for(int i=0; i<Nr; i++)
    {
      tba_ly[i] = tba_single<TBA_MODEL>(Ns);
      tba_ly[i].init_all(30,r[i]);
    }
  }

  cfunc_single(int Ns, double r_min, double r_max, double dr)  // constructor range and dr
  {
    Nr = (int) (r_max-r_min)/dr;
    r = new double[Nr];
    c = new double[Nr];
    tba_ly = new tba_single<TBA_MODEL>[Nr];
    for(int i=0; i<Nr-1; i++){
      r[i] = r_min + i*dr;
    }
    r[Nr-1] = r_max;
    for(int i=0; i<Nr-1; i++){
      tba_ly[i].init_all(30,r[i]);
    }
  }

  // Methods
  void evaluate_cfun()                                         // solves the tba's and evaluates c(r)
  {
    for(int i=0; i<Nr; i++)
    {
      tba_ly[i].solve();
      tba_ly[i].evaluate_c();
      c[i] = tba_ly[i].c;
    }
  }

  void save_results(std::string filename)                      // save c(r) to file
  {
    std::ofstream out(filename);
    out << "index;r;c" << std::endl;
    for(int i=0; i<Nr; i++) out << i << ";" << r[i] << ";" << c[i] << std::endl;
    out.close();
  }

};
