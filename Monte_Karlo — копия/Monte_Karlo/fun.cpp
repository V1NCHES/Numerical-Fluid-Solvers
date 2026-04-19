#include"include.h"

double normal_extractionem_0_1(int & a1_, int &a2_, int & a3_ ){
	int ic15 = 32768, ic10 = 1024;
	int mz = 710, my = 17784, mx = 11973;
	double xi = 9.0949470177292824E-13, c = 1.073741824E9;
	double b;
	int i13, i12, i11, ii;
	i13 = mz * a1_ + my * a2_ + mx * a3_;
	i12 = my * a1_ + mx * a2_;
	i11 = mx * a1_;
	ii = i11 / ic15;
	i12 = i12 + ii;
	a1_ = i11 - ic15 * ii;
	ii = i12 / ic15;
	i13 = i13 + ii;
	a2_ = i12 - ic15 * ii;
	a3_ = i13 % ic10;
	b = xi * (c * a3_ + ic15 * a2_ + a1_);
	return b;
}
double extractionem_1_1()
{
//	double result = -1 +2*(double)rand() / RAND_MAX;
//	while (is_equal(-1, result) || is_equal(1, result))
	{
		///result = -1 + 2 *( (double)rand() / RAND_MAX);
	}
	return -1 + 2 * (double)rand() / RAND_MAX;
}
double extractionem_0_1_()
{
	return (double)rand() / RAND_MAX;
}
double extractionem_0_1()
{
	double result = (double)rand() / RAND_MAX;
	while (is_equal(0, result) || is_equal(1, result))
	{
		result = (double)rand() / RAND_MAX;
	}
	return result;
}
bool is_equal(double x, double y) {
	return std::fabs(x - y) < std::numeric_limits<double>::epsilon();
}

double sqrt_check(double x, int index)
{
	
	if (x < 0)
	{
		cout << index<< " ERROR sqrt x < 0\n";
		return 0;
	}
	return sqrt(x);
}
bool check_inf(double V)
{
	return -inf < V && V < inf;
}
double log_check(double x)
{
	if (!(check_inf(x)) || is_equal(0, x))
	{
		cout << "ERROR log : x = 0  || x = -+inf\n";
		return 0;
	}
	return log(x);
}

