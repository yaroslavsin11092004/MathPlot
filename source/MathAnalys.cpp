#include "MathAnalys.h"
double DownRectIntegral(matrix<double> &mesh)
{
	size_t n = mesh.size_row();
	double ans = 0;
	for (size_t i = 0; i < n-1; i++)
	{
		double h = mesh(i + 1,0)- mesh(i,0);
		ans += mesh(i,1) * h;
	}
	return ans;
}
double UpRectIntegral(matrix<double> &mesh)
{
	size_t n = mesh.size_row();
	double ans = 0;
	for (size_t i = 0; i < n - 1; i++)
	{
		double h = mesh(i + 1,0) - mesh(i,0);
		ans += mesh(i + 1,1) * h;
	}
	return ans;
}
double CentrRectIntegral(matrix<double> &mesh)
{
	size_t n = mesh.size_row();
	double ans = 0;
	for (size_t i = 0; i < n - 1;i++)
	{
		double mWidth = mesh(i + 1,1) >= mesh(i,1) ? mesh(i,1) : mesh(i + 1,1);
		double sWidth = mWidth == mesh(i + 1,1) ? mesh(i,1) : mesh(i + 1,1);
		double delta = (sWidth - mWidth) / 2;
		double h = mesh(i + 1,0) - mesh(i,0);
		ans += (sWidth - delta) * h;
	}
	return ans;
}
double TrapezeIntegral(matrix<double> &mesh)
{
	size_t n = mesh.size_row();
	double ans = 0;
	for (size_t i = 0; i < n - 1; i++)
	{
		double h = mesh(i + 1,0) - mesh(i,0);
		ans += (mesh(i + 1,1) + mesh(i,1)) * h;
	}
	return ans / 2;
}
double SimpsonIntegral(matrix<double> &mesh)
{
	size_t n = mesh.size_row();
	size_t re_size = (n%2 == 0) ? n / 2 : n / 2 + 1;
	matrix<double> _mesh(re_size, 2);
	for (size_t i = 0;i < re_size;i++)
		_mesh.set_row(i, mesh.make_row_acceptor(2 * i));
	double val1 = TrapezeIntegral(mesh);
	double val2 = TrapezeIntegral(_mesh);
	return val1 + (val1 - val2) / 3;
}
double MonteCarloIntegral(matrix<double> &mesh)
{
	size_t n = mesh.size_row();
	double h = (mesh(n-1,0) - mesh(0,0)) / n;
	return mesh.sum_column(1) * h;
}
size_t BinarySearch(double value, matrix<double> &mesh)
{
	size_t start = 0;
	size_t end = mesh.size_row();
	size_t ans = end + 1;
	int data = static_cast<int>(value);
	while(start <= end)
	{
		size_t c = (start + end) / 2;
		if (mesh(c,0) >= data)
		{
			ans = c;
			end = c - 1;
		}
		else 
			start = c + 1;
	}
	return ans;
}
double DifferentSubtract(matrix<double> &mesh)
{
	size_t n = mesh.size_row();
	double sum = 0;
	for (size_t i = 0; i < n; i++)
	{
		double mul = 1;
		for (size_t j = 0; j < n; j++)
			if (j != i) mul *= 1.0 / (mesh(i,0) - mesh(j,0));
		sum += mesh(i,1) * mul;
	}
	return sum;
}
int Factorial(int n)
{
	return (n == 1 || n == 0) ? 1 : Factorial(n-1) * n;
}
double Differential(double x, int order, int offset, matrix<double> &mesh)
{
	size_t n = mesh.size_row();
	double answer = 0;
	size_t index = BinarySearch(x,mesh);
	if (index - order * offset < 0 || index + order* offset > n)
		throw std::runtime_error("Error of order or offset");
	else
	{
		matrix<double> m_right(order + 1,2);
		matrix<double> m_left(order + 1,2);
		size_t l_index = 0, r_index = 0;
		for (size_t i = index - order * offset; i < index + order* offset + 1; i += offset)
		{
			if (i < index)
			{
				m_left.set_row(l_index, mesh.make_row_acceptor(i));
				l_index++;
			}
			else if (i > index)
			{
				m_right.set_row(r_index, mesh.make_row_acceptor(i));
				r_index++;
			}
			else 
			{
				m_left.set_row(l_index, mesh.make_row_acceptor(i));
				m_right.set_row(r_index, mesh.make_row_acceptor(i));
				r_index++;
			}
		}
		double left_dif = DifferentSubtract(m_left);
		double right_dif = DifferentSubtract(m_right);
		return (left_dif + right_dif) / 2 * Factorial(order);
	}
}
double ExactDifferential(double x,int order, matrix<double> &mesh)
{
	try 
	{
		double dif1 = Differential(x, order, 1,mesh);
		double dif2 = Differential(x, order, 2, mesh);
		return dif1 + (dif1 - dif2) / 3;
	} 
	catch(const std::exception& e)
	{
		throw std::runtime_error(e.what());
	}
}
double DifferentBack(double x, matrix<double> &mesh)
{
	size_t index = BinarySearch(x, mesh);
	if (index == -1 || index - 1 < 0)
		throw std::runtime_error("Error of search value");
	return (mesh(index,1) - mesh(index - 1,1)) / (mesh(index, 0) - mesh(index - 1,0));
}
double DifferentForward(double x, matrix<double> &mesh)
{
	size_t index = BinarySearch(x,mesh);
	if (index == - 1 || index + 1 > mesh.size_row() - 1)
		throw std::runtime_error("Error of search value");
	return (mesh(index + 1,1) - mesh(index,1)) / (mesh(index + 1,0) - mesh(index, 0));
}
double DifferentCenter(double x, matrix<double> &mesh)
{
	size_t index = BinarySearch(x, mesh);
	if (index == -1 || index + 1 > mesh.size_row() - 1 || index - 1 < 0)
		throw std::runtime_error("Error of search value");
	return (mesh(index + 1,1) - mesh(index - 1,1)) / (mesh(index + 1,0) - mesh(index - 1,0));
}
