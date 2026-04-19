#ifndef distribution_function_h
#define distribution_function_h 
#include"include.h"
#include "velocity_at_border.h"

class Interval
{
public:
	double a, b;
	int n;
	Interval() { a = 0, b = 0; n = 0; }
	Interval(double a1, double b1, int n1) { a = a1; b = b1; n = n1; }
	void calculete_n(double d){ n = (int)((b - a) / d); }
	friend ostream& operator<<(ostream& cout, Interval& P);//îïðàòîð âûâîäà
};
class dataFH
{
public:
	double t;
	long long int N;
	dataFH() { t = 0, N = 0;}
	dataFH(double t1, long long  int N1) {t = t1; N = N1; }	
};
class ÑMatrix
{
public:
	//vector<vector<vector<vector<double>>>> vec;
	vector<vector<double>> vec;
	//vector<long long int> kolL;
	Interval x,vx,vy,vz;
	double dx,dy, dl;

	ÑMatrix()
	{
		dx = 0.1; dy = 0.1; dl = 0.05;
		x.a = -L, x.b = 0; x.calculete_n(dl);
		//vx.a = -3, vx.b = 3;  vx.calculete_n(d);
		//vy.a = -3, vy.b = 3;  vy.calculete_n(d);
		//vz.a = -3, vz.b = 3;  vz.calculete_n(d);
		vx.a = -6, vx.b = 6;  vx.calculete_n(dx);
		vy.a = -3, vy.b = 3;  vy.calculete_n(dy);
		vz.a = -3, vz.b = 3;  vz.calculete_n(dy);
		memory_allocation();
	}
/*	~ÑMatrix()
	{ 
		//cout << "HIIIIIIIIIII\n";
		if (x.n > 0)
		{
			x.n = 0; vx.n = 0; vy.n = 0; vz.n = 0;
			for (int i = 0; i < vz.n; ++i)
			{
				
				for (int j = 0; j < vy.n; ++j)
				{
					
					for (int k = 0; k < vx.n; ++k)
					{
						vec[i][j][k].clear();
						
					}
					vec[i][j].clear();
				}
				
				vec[i].clear();
			}
			vec.clear();
		}
		//cout << "1 HIIIIIIIIIII\n";
	}//*/
	ÑMatrix(Interval x1, Interval vx1, Interval vy1, Interval vz1, double d1)
	{
		x = x1; dx = d1, dy = dl; dl = d1; x.calculete_n(dl); vx.calculete_n(dx); //vx = vx1; vy = vy1; vz = vz1;   vy.calculete_n(dy); vz.calculete_n(dy);
	}
//	double integrate_vy_vz_Matrix(int k, int l);
	void memory_allocation()
	{
		//kolL.resize(x.n,0);
		//int kol = vy.n * vz.n;
	//	
		//for (int i = 0; i < x.n*vx.n; ++i)
		//for (int i = 0; i < vz.n; ++i)
		{
			//vector <double> martrixR3(kol, 0);
		//	vector <vector <vector<double>>> martrixR3;
			//cout << "i = " << i << endl;
		//	for (int j = 0; j < vy.n; ++j)
			{
			//	vector <vector<double>> martrixR2;
				for (int k = 0; k < x.n+1; ++k)
				{
					
					vec.push_back(vector<double> (vx.n, 0));
					
				}
				//martrixR3.push_back(martrixR2);
			}
			//*/
			//vec.push_back(martrixR3);
		}
	}
	void matrix_vx_x_fH(vector<vector<double>> &Matrix);

//	void distribution_function_vx(size_t N);
	void filling_matrix(double x, double x0, double t, Velociti_Atom& V, double v, int N, vector<vector<double>> &graficV,double Le, int &err, int &KOL);
	bool check_beloning(Velociti_Atom& V);
	
	void printN()
	{ 
		cout <<  vz.n << " ; " << vy.n <<" ; "<< vx.n << " ; " << x.n <<"\n"; 
		cout << "vz = {" << vz.a << " ; " << vz.b << " } " << "vx = {" << vx.a << " ; " << vx.b << " }\n";
	}
};


double A();
//double A(double min_lim, double  max_lim);

//double func(double a);
//double integrate(double(*func)(double), double min_lim, double  max_lim, int n);
//double method_of_rectangles(double(*func)(double), double min_lim, double max_lim, double delta);

void OutPut_distribution_function_vx(vector<vector<double>>& grafic, ÑMatrix& M, int N, string& name, string& name1);
void nameFile(string& name, int index, int i, int j,int k);
#endif
