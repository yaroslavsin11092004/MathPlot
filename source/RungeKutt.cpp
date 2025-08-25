#include "RungeKutt.h"
matrix<double> runge_kutt(RkRef args)
{
	auto cond = [args](double ti)->bool
	{
		if (args->get().dir) return ti <= args->get().t1 + args->get().step ? true : false;
		else return ti >= args->get().t0 - args->get().step ? true : false;
	};
	double mul = args->get().dir == true ? 1 : -1;
	size_t dim_vec = args->get().u0->size();
	std::vector<std::vector<double>> var(dim_vec + 1);
	std::vector<double> cur(dim_vec + 1);
	cur[0] = args->get().t0;
	for (size_t i = 0;i < dim_vec; i++)
		cur[i + 1] = (*args->get().u0)[i];
	do 
	{
		for (size_t i = 0; i < dim_vec + 1;i++)
			var[i].push_back(cur[i]);
		matrix<double> params(dim_vec, args->get().dim_var);
		std::vector<double> start_value(dim_vec);
		double t = cur[0];
		auto calc_param_func = [args, &params](double t)->void 
		{
			if (args->get().param_func != nullptr)
			{
				for (size_t i = 0; i < args->get().param_func->size(); i++)
					params(args->get().param_func_indeces[i].first, args->get().param_func_indeces[i].second) = value((*args->get().param_func)[i].ratios, (*args->get().param_func)[i].mesh, t, (*args->get().param_func)[i].dir);
			}
		};
		for (size_t i = 0; i < dim_vec; i++)
			start_value[i] = cur[i + 1];
		params.set_column(args->get().index_var, start_value);
		calc_param_func(t);
		matrix<double> k(4, dim_vec);
		for (size_t i = 0;i < dim_vec;i++)
			k(0,i) = (*args->get().equations)[i](t, params);
		t = cur[0] + mul * args->get().step / 2;
		for (size_t i = 0; i < dim_vec; i++)
			params(i, args->get().index_var) = start_value[i] + mul * args->get().step / 2 * k(0,i);
		calc_param_func(t);
		for (size_t i = 0; i < dim_vec; i++)
			k(1,i) = (*args->get().equations)[i](t,params);
		for (size_t i = 0; i < dim_vec; i++)
			params(i, args->get().index_var) = start_value[i] + mul * args->get().step / 2 * k(1,i);
		for (size_t i = 0; i < dim_vec; i++)
			k(2,i) = (*args->get().equations)[i](t, params);
		t = cur[0] + mul * args->get().step;
		for (size_t i = 0; i < dim_vec;i++)
			params(i, args->get().index_var) = start_value[i] + mul * args->get().step * k(2,i);
		calc_param_func(t);
		for (size_t i = 0; i < dim_vec; i++)
			k(3,i) = (*args->get().equations)[i](t, params);
		cur[0] += args->get().step * mul;
		for (size_t i = 0; i < dim_vec; i++)
			cur[i + 1] += mul * args->get().step / 6 * (k(0,i) + 2 * (k(1,i) + k(2,i)) + k(3,i));
	}
	while(cond(cur[0]));
	matrix<double> ans(var[0].size(), dim_vec + 1);
	for (size_t i = 0; i < dim_vec + 1; i++)
		ans.set_column(i, var[i]);
	return ans;
}
