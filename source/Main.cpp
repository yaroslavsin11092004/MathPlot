#include "matrix.h"
using namespace std;
int main()
{
	matrix<double> m(2,2);
	matrix<double> l(2,2);
	m.set_row(0, {1,2});
	m.set_row(1, {3,4});
	l.set_column(0, {5,6});
	l.set_column(1, {7,8});
	cout << m << endl << l << endl;
	m.set_column(0, l.make_column_acceptor(0));
	cout << m << endl;
	return 0;
}
