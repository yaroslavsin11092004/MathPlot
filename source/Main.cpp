#include "RungeKutt.h"
#include "interpolation.h"
double func_1(double var, matrix<double> params)
{
	return var;
}
double func_2(double var, matrix<double> params)
{
	return params(0,0) + var / 2.0;
}
double f(double x) { return pow(x,3) / 6.0 + pow(x,2) / 4.0;}
using namespace std;
int main()
{
	vector<function<double(double, matrix<double>)>> eq = {func_1, func_2};
	vector<double> x0({0.0,0.0});
	rk_params task;
	task.param_func = nullptr;
	task.step = 0.01;
	task.equations = make_shared<vector<function<double(double, matrix<double>)>>>(std::move(eq));
	task.index_var = 0;
	task.u0 = make_shared<vector<double>>(std::move(x0));
	task.dir = true;
	task.t1 = 2;
	task.dim_var = 1;
	task.t0 = 0.0;
	matrix<double> ans = runge_kutt(ref(task));
	matrix<double> f0(ans.size_row(), 2);
	f0.set_column(0,ans.make_column_acceptor(0));
	f0.set_column(1, ans.make_column_acceptor(2));
	auto ratios = spline_interpolation(f0);
	interpol_param ip; 
	ip.ratios = std::move(ratios);
	ip.mesh = std::move(f0);
	ip.dir = true;
	cout << value(ip.ratios, ip.mesh, 1.015,ip.dir) << endl << f(1.015);
	return 0;
}
