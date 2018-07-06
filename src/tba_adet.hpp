// Copyright Alessandro Fabbri 2015

#include <cmath>
#include <tba.hpp>

/*
 *  ADET TBA class with full features
 */
class tba_adet : public tba
{
public:
  // Additional parameters
  char alg_type;                                                             // A,D,E,T type switch
  int Nrank;                                                                 // number of nodes
  int mass_pos;                                                              // index of massive node
  // Functions
  double ** eps;                                                             // pseudoenergy
  double ** L;                                                               // log(1 + e^-epsilon)
  double ** nu;                                                              // forcing term
  double * phi;                                                              // kernel
  double * z;                                                                // counterterms
  // Constructors
  tba_adet() {}                                                              // default constructor

  ~tba_adet() {}                                                             // default destructor

  tba_adet(const tba_adet &adet)                                             // copy constructor
  {
    Nstep = adet.Nstep;
    Nrank = adet.Nrank;
    iter  = adet.iter;
    error = adet.error;
    x     = adet.x;
    eps   = adet.eps;
    L     = adet.L;
    nu    = adet.nu;
    phi   = adet.phi;
  }

  tba_adet(int Ns, char Alg, int Rank, int Index)                            // step size constructor
  {
    alg_type = Alg;
    Nstep = Ns;
    Nrank = Rank;
    if ( Index != 0 )
    {
      std::cerr << "WARNING: Mass index != 0 coming soon..." << std::endl;
      Index = 0;
    }
    mass_pos = Index;
    x   = new double [Nstep];
    phi = new double [2*Nstep-1];
    z   = new double [Nrank];
    eps = new double* [Nrank];
    L   = new double* [Nrank];
    nu  = new double* [Nrank];
    for(int a=0; a<Nrank; ++a)
    {
      eps[a] = new double [Nstep];
      L[a]   = new double [Nstep];
      nu[a]  = new double [Nstep];
    }
  }

  void operator=(const tba_adet &adet)                                       // assignement operator
  {
    Nstep = adet.Nstep;
    iter  = adet.iter;
    error = adet.error;
    x     = adet.x;
    eps   = adet.eps;
    L     = adet.L;
    nu    = adet.nu;
    phi   = adet.phi;
  }

  // Methods
  void init_all(double xmax, double r_)                                      // initial operations, divided into:
  {
    iter  = 0;
    error = 1;
    r     = r_;
    c     = 0;
    x_max = ( r < 1e-2 ? 2.*std::abs(std::log(2/r)) : xmax );
    x_min = -x_max;
    dx = (x_max-x_min)/(Nstep-1.);
    for(int i=0; i<Nstep; ++i) x[i] = x_min+i*dx;

    setup << alg_type << Nrank << " TBA setup" << std::endl
          << "Algebra    = " << alg_type << "_" << Nrank << std::endl
          << "Mass index = " << mass_pos << std::endl
          << "Nstep      = " << Nstep << std::endl
          << "X Range    = [ "  << x_min << " , " << x_max << " ]" << std::endl
          << "dx         = " << dx << std::endl
          << "r          = " << r <<  std::endl;

    init_kernels();
    init_guess();
  }

  void init_kernels()                                                        // convolution kernel init and saving to file
  {
    double h=0;
    switch (alg_type)
    {
      case 'A':
      {
        h = Nrank+1.;
        double alpha=M_PI/(Nrank+2.);
        for(int a=0; a<Nrank; ++a)
        {
          if(a==0) z[a] = 0;
          else     z[a] = std::sin((a+2)*alpha)*std::sin(a*alpha)/std::sin((Nrank+1)*alpha)/std::sin(alpha);
        }
        break;
      }
      case 'D':
      {
        std::cerr << "Coming soon" << std::endl;
        break;
      }
      case 'E':
      {
        std::cerr << "Coming soon" << std::endl;
        break;
      }
      case 'T':
      {
        std::cerr << "Coming soon" << std::endl;
        break;
      }
      default:
      {
        std::cerr << "Coming soon" << std::endl;
        break;
      }
    }

    double xi;
    for(int i=0; i<2*Nstep-1; ++i){
      xi = 2.*x_min + i*dx;
      phi[i] = h/4/M_PI/std::cosh(h*xi/2);
    }
  }

  void init_guess()                                                          // initial guess for epsilon
  {
    for(int a=0; a<Nrank; ++a)
    {
      for(int i=0; i<Nstep; ++i)
      {
        if(a==mass_pos) nu[a][i] = r*std::cosh(x[i]);
        else            nu[a][i] = 0;
        eps[a][i] = nu[a][i];
        L[a][i] = std::log(1+std::exp(-eps[a][i]));
      }
    }
  }

  void solve()                                                               // solution loop
  {
    double *temp_sum, temp0; temp_sum = new double [Nrank];
    do{
      for(int i=0; i<Nstep; ++i)
      {
        for(int a=0; a<Nrank; ++a) temp_sum[a]=0;
        for(int j=0; j<Nstep; ++j)
        {
          int ij = Nstep-1+std::abs(i-j);
          temp_sum[0] += phi[ij]*(L[1][j]-std::log(1+z[1]));
          for(int a=1; a<Nrank-1; ++a)
          {
            temp_sum[a] += phi[ij]*(L[a-1][j]+L[a+1][j]-std::log(1+z[a-1])-std::log(1+z[a+1]));
          }
          temp_sum[Nrank-1] += phi[ij]*(L[Nrank-2][j]-std::log(1+z[Nrank-2]));
        }
        eps[0][i] =.5*(eps[0][i] + nu[0][i] - dx*temp_sum[0]-.5*std::log(1+z[1]));
        for(int a=1; a<Nrank-1; ++a)
        {
          eps[a][i] =.5*(eps[a][i] + nu[a][i] - dx*temp_sum[a]-.5*std::log(1+z[a-1])-.5*std::log(1+z[a+1]));
        }
        eps[Nrank-1][i] =.5*(eps[Nrank-1][i] + nu[Nrank-1][i] - dx*temp_sum[Nrank-1]-.5*std::log(1+z[Nrank-2]));
        for(int a=0; a<Nrank; ++a) L[a][i] = std::log(1+std::exp(-eps[a][i]));
      }
      error = std::abs(eps[0][(Nstep+1)/2]-temp0);
      temp0 = eps[0][(Nstep+1)/2];
      ++iter;
    } while( error > PRECISION );
  }

  void evaluate_c()                                                          // c-function evaluation
  {
    for(int i=0; i<Nstep-1; ++i)
      c += .5*(nu[mass_pos][i]*L[mass_pos][i]+nu[mass_pos][i+1]*L[mass_pos][i+1]);
    c *= dx*3/(M_PI*M_PI);

    report << alg_type << Nrank << " TBA report" << std::endl
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
    for(int i=0; i<2*Nstep-1; ++i)
    {
      xi = 2.*x_min + i*dx;
      out << i << ";" << xi << ";" << phi[i] << std::endl;
    }
    out.close();
  }

  void save_results(std::string filename)                                    // save results to file
  {
    std::ofstream out(filename);
    out << "index;x;";
    for(int a=0; a<Nrank; ++a) out << "eps" << a << ";L" << a << ( (a != Nrank-1) ? ";" : "");
    out << std::endl;
    for(int i=0; i<Nstep; ++i)
    {
      out << i << ";" << x[i] << ";";
      for(int a=0; a<Nrank; ++a) out << eps[a][i] << ";" << L[a][i] << ( (a != Nrank-1) ? ";" : "");
      out << std::endl;
    }
    out.close();
  }
};

/*
 *  C-Function ADET class for automated evaluation of c(r)
 */
class cfunc_adet{
public:
  // Variables
  int r_sample;
  double * r;
  double * c;
  tba_adet* tba_ly;
  // Constructors
  cfunc_adet();                                                         // default constructor
  ~cfunc_adet();                                                        // default destructor
  cfunc_adet(int, int, double, double, int);                            // constructor range and steps
  cfunc_adet(int, int, double, double, double);                         // constructor range and dr
  // Methods
  void evaluate_cfun();                                                 // solves the tba's and evaluates c(r)
  void save_results(std::string filename); // save c(r) to file
};