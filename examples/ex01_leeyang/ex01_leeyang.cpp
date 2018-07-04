#include <iostream>
#include <tba_single.hpp>

using namespace std;

int main(){
	cout << "*********** LEE-YANG TBA EXAMPLE ***********" << endl;

	tba_single<TBA_LEEYANG> tba(1000);
	tba.init_all(10.,1e-8);
	tba.show_setup();
	tba.show_setup("setup.txt");
	tba.save_kernels("kernel.csv");
	tba.solve();
	tba.evaluate_c();
	tba.show_report();
	tba.show_report("report.txt");
	tba.save_results("results.csv");

	cfunc_single<TBA_LEEYANG> cfunc(1000, 1e-8, 1, 20);
	cfunc.evaluate_cfun();
	cfunc.save_results("cfunc.csv");

	return 0;
}
