#include "RungeKutt.h"

matrix<double> runge_kutt(RkRef args)
{
	auto cond = [args](double ti)->bool
	{
		if (args->get().dir) return ti <= args->get().t1 ? true : false;
		else return ti >= args->get().t0 ? true : false;
	};
	double mul = args->get().dir == true ? 1 : -1;
	vector<vector<double>> var(args->get().u0->size() + 1);
	valarray<double> cur(args->get().u0->size() + 1);
	cur[0] = args->get().t0;
	for (size_t i = 0;i < args->get().u0->size(); i++)
		cur[i + 1] = (*args->get().u0)[i];
	do 
	{
		for (size_t i = 0; i < cur.size();i++)
			var[i].push_back(cur[i]);
		valarray<double> var1(args->get().size_var);
		var1[(*args->get().indeces)[0]] = cur[0];
		size_t offset = args->get().param_func->size();
		for (size_t i = 0; i < offset; i++)
			var1[(*args->get().indeces)[i + 1]] = value((*args->get().param_func)[i].ratios, (*args->get().param_func)[i].mesh, var1[(*args->get().indeces)[0]], (*args->get().param_func)[i].dir);
		for (size_t i = 0; i< args->get().u0->size();i++)
			var1[(*args->get().indeces)[i + 1 + offset]] = cur[i + 1];
		matrix<double> k(4, args->get().u0->size());
		for (size_t i = 0;i < args->get().u0->size();i++)
			k[0][i] = (*args->get().equations)[i](var1);
		valarray<double> var2(var1);
		var2[(*args->get().indeces)[0]] += mul * args->get().step / 2;
		for (size_t i = 0; i < offset; i++)
			var2[(*args->get().indeces)[i + 1]] = value((*args->get().param_func)[i].ratios, (*args->get().param_func)[i].mesh, var2[(*args->get().indeces)[0]], (*args->get().param_func)[i].dir);
		for (size_t i = 0; i < args->get().u0->size();i++)
			var2[(*args->get().indeces)[i + 1 + offset]] += mul * args->get().step / 2 * k[0][i];
		for (size_t i = 0; i < args->get().u0->size();i++)
			k[1][i] = (*args->get().equations)[i](var2);
		valarray<double> var3(var1);
		var3[(*args->get().indeces)[0]] += mul * args->get().step / 2;
		for (size_t i = 0; i < offset;i++)
			var3[(*args->get().indeces)[i + 1]] = value((*args->get().param_func)[i].ratios, (*args->get().param_func)[i].mesh, var3[(*args->get().indeces)[0]], (*args->get().param_func)[i].dir);
		for (size_t i = 0; i < args->get().u0->size(); i++)
			var3[(*args->get().indeces)[i + 1 + offset]] += mul * args->get().step / 2 * k[1][i];
		for (size_t i = 0; i < args->get().u0->size(); i++)
			k[2][i] = (*args->get().equations)[i](var3);
		valarray<double> var4(var1);
		var4[(*args->get().indeces)[0]] += mul * args->get().step;
		for (size_t i = 0; i < offset; i++)
			var4[(*args->get().indeces)[i + 1]] = value((*args->get().param_func)[i].ratios, (*args->get().param_func)[i].mesh, var4[(*args->get().indeces)[0]], (*args->get().param_func)[i].dir);
		for (size_t i = 0; i < args->get().u0->size();i++)
			var4[(*args->get().indeces)[i + 1 + offset]] += mul * args->get().step * k[2][i];
		for (size_t i = 0; i < args->get().u0->size();i++)
			k[3][i] = (*args->get().equations)[i](var4);
		cur[0] += args->get().step * mul;
		for (size_t i = 0; i < args->get().u0->size();i++)
			cur[i + 1] += mul * args->get().step / 6 * (k[0][i] + 2 * k[1][i] + 2 * k[2][i] + k[3][i]);
	}
	while(cond(cur[0]));
	matrix<double> res(var[0].size(), args->get().u0->size() + 1);
	for (size_t i = 0; i < args->get().u0->size() * 2; i++)
	{
		for (size_t j = 0; j < var[i].size();j++)
			res[j][i] = var[i][j];
	}
	return res;
}
