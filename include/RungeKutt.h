#ifndef RUNGEKUTT_H
#define RUNGEKUTT_H
#include "interpolation.h"
#include "Matrix.h"
#include <optional>
#include <functional>
#include <memory>
struct rk_params
{
	shared_ptr<vector<function<double(valarray<double>)>>> equations = nullptr;
	shared_ptr<vector<interpol_param>> param_func = nullptr;
	double t0 = 0.0;
	double t1 = 0.0;
	shared_ptr<valarray<double>> u0 = nullptr;
	double step = 0.0;
	shared_ptr<vector<size_t>> indeces = nullptr;
	bool dir = true;
	size_t size_var;
	rk_params() = default;
};
using RkRef = optional<reference_wrapper<rk_params>>;
matrix<double> runge_kutt(RkRef args);
#endif
