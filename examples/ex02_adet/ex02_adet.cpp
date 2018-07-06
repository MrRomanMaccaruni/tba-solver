#include <iostream>
#include <tba_adet.hpp>

using namespace std;

int main(){
	cout << "*********** ADET TBA EXAMPLE ***********" << endl;

	tba_adet adet(1000, 'A', 5, 1);
	adet.init_all(30, 1e-8);
	adet.show_setup();
	adet.save_setup("setup.txt");
	adet.save_kernels("kernel.csv");
	adet.solve();
	adet.evaluate_c();
	adet.show_report();
	adet.save_report("report.txt");
	adet.save_results("results.csv");

	return 0;
}
