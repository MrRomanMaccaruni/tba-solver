// Copyright Alessandro Fabbri 2015

#include <fstream>
#include <sstream>

#define PRECISION 1e-15

#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define DEFAULT_NSTEP  1000
#define DEFAULT_XMAX   10.
#define DEFAULT_R      1E-8

/*
 *  Generic TBA class containing the minimal amount of
 *  parameters which can belong to any TBA model and no methods
 */
class tba{
public:
  int Nstep, iter;
  double x_min, x_max, dx;
  double * x;
  double error;
  double r, c;
  std::stringstream setup;
  std::stringstream report;

  void show_setup()                          // show setup values on stdout
  {
    std::cout << setup.str();
  }

  void save_setup(std::string filename)      // save setup values to file
  {
    std::ofstream out(filename);
    out << setup.str();
    out.close();
  }

  void show_report()                         // show report on stdout
  {
    std::cout << report.str();
  }

  void save_report(std::string filename)     // save report to file
  {
    std::ofstream out(filename);
    out << report.str();
    out.close();
  }
};
