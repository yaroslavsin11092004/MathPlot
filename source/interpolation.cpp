#include "interpolation.h"
matrix<double> running(matrix<double>& A, matrix<double>& D)
{
	int size = static_cast<int>(A.size_row());
	matrix<double> ktx(size + 1, 3);
	for (int i = 0;i < size;i++)
	{
		double a,b,c,d;
		d = D(i,0);
		if (i == 0)
		{
			a = 0;
			b = -A(i,i);
			c = A(i,i + 1);
		}
		else if (i == size - 1)
		{
			a = A(i,i-1);
			b = -A(i,i);
			c = 0;
		}
		else
		{
			a = A(i,i-1);
			b = -A(i,i);
			c = A(i,i + 1);
		}
		ktx(i + 1, 0) = c / (b - ktx(i,0) * a);
		ktx(i + 1, 1) = (a * ktx(i,1) - d) / (b - a * ktx(i,0));
	}
	matrix<double> res(size,1);
	int index = 0;
	for (int i = size - 1; i >= 0; --i)
	{
		ktx(i,2) = ktx(i + 1, 0) * ktx(i + 1,2) + ktx(i + 1,1);
		res(index,0) = ktx(i,2);
		index++;
	}
	return res;
}
matrix<double> spline_interpolation(matrix<double>& grid)
{
	int N = static_cast<int>(grid.size_row());
	matrix<double>A(N,N);
	matrix<double>D(N,1);
	A(0,0) = -1;
	for (int i = 2; i < N;i++)
	{
		double h0 = grid(i-1,0) - grid(i-2,0);
		double h1 = grid(i,0) - grid(i-1,0);
		A(i-1,i-1) = 2 * (h0 + h1);
		A(i-1,i) = h1;
		A(i-1,i-2) = h0;
		D(i-1,0) = 3 * ((grid(i,1) - grid(i-1,1)) / h1 - (grid(i-1,1) - grid(i-2,1))/h0);
	}
	A(N-1,N-1) = 1;
	D(N-1,0) = 0;
	matrix<double> C = running(A,D);
	matrix<double> result(N-1, 4);
	for (int i = 0;i < N-1;i++)
	{
		double h = grid(i + 1,0) - grid(i,0);
		result(i,0) = grid(i,1);
		result(i,1) = (grid(i + 1,1) - grid(i,1)) / h - h / 3 * (C(i + 1,0) + 2 * C(i,0));
		result(i,2) = C(i,0);
		result(i,3) = (C(i + 1,0) - C(i,0)) / (3 * h);
	}
	return result;
}
size_t get_range(matrix<double>& grid, double x, bool dir)
{
	size_t start, end, center;
	if (dir)
	{
		start = 0;
		end = grid.size_row() - 1;
		if (x >= grid(grid.size_row()-1,0)) return grid.size_row()-2;
		if (x <= grid(0,0)) return 0;
		while(start <= end)
		{
			center = (start + end) / 2;
			if (x <= grid(center,0) && x >= grid(center-1,0))
				return center - 1;
			else if (x > grid(center,0))
				start = center + 1;
			else if (x < grid(center-1,0))
				end = center - 1;
		}
		return -1;
	}
	else 
	{
		if (x <= grid(grid.size_row()-1,0))return grid.size_row()-2;
		if (x >= grid(0,0)) return 0;
		start = grid.size_row()-1;
		end = 0;
		while(end <= start)
		{
			center = (start + end) / 2;
			if (x >= grid(center,0) && x <= grid(center-1,0))
				return center-1;
			else if (x < grid(center,0))
				end = center + 1;
			else if (x > grid(center-1,0))
				start = center - 1;
		}
		return -1;
	}
}
double value(matrix<double> & ratios, matrix<double> & grid, double x, bool dir)
{
	size_t index = get_range(grid, x, dir);
	if (index == -1) return -1;
	else 
	{
		double z = x - grid(index,0);
		return ratios(index,0) + ratios(index,1) * z + ratios(index,2) * pow(z,2) + ratios(index,3) * pow(z,3);
	}
}
matrix<double> make_draw_matrix(matrix<double> & ratios, matrix<double> & grid, double step, bool dir)
{
	double x0 = grid(0,0);
	double xn = grid(grid.size_row()-1,0);
	std::vector<double> argX, argY;
	while(x0 < xn)
	{
		argX.push_back(x0);
		argY.push_back(value(ratios, grid, x0, dir));
		x0 += step;
	}
	argX.push_back(xn);
	argY.push_back(value(ratios, grid, xn, dir));
	matrix<double> res(argX.size(),2);
	res.set_column(0, argX);
	res.set_column(1, argY);
	return res;
}
