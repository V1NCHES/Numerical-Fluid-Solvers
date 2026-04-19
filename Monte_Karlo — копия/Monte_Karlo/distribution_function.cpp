#include"distribution_function.h"
#include "include.h"
namespace fs = std::experimental::filesystem;
ostream& operator<<(ostream& cout, Interval& P)
{
	cout << P.a << "\n" << P.b;
	return cout;
}

bool СMatrix::check_beloning(Velociti_Atom& V)
{
	return (V.x > vx.a && V.x < vx.b) && (V.y > vy.a && V.y < vy.b) && (V.z > vz.a && V.z < vz.b);
}
void СMatrix::filling_matrix(double x, double x0, double t, Velociti_Atom& V, double v, int N, vector<vector<double>> &graficV,double Le,int &err,int &KOL)
{
//	double C = A() / (N * dx * dy * dy * dl);
	int l = 0, p = 0; //  индексы границы L
	//double T = dl / v; // время находжения в полной ячейке
	double modV = A()  /(v * N * dx);
	//double modV = 1. /v;
	//СMatrix A;

	bool end = true; // для проверки где находится x 
	//int kol = vy.n * vz.n;
	int k = (int)(fabs(this->vx.a - V.x) / dx); // высчитываем индексы по скоростям
	//int j = (int)(fabs(this->vy.a - V.y) / dy);
	//int i = (int)(fabs(this->vz.a - V.z) / dy);
	//double tl = 0; // время кусочка клетки
	//cout << i << " " << j << " " << k << endl;
	

	//if (Le < dl)  cout << "p = " << p << "; l =" << l << "; k =" << k << "length < M.dl\n";

	if (Le < dl )
	{
	//	cout << "x0 = " << x0 << "; x = " << x << endl;
		if (x > 0)
		{
		//	cout << "x>0\n\n";
			vec[this->x.n][k] += modV;
			//graficV[[this->x.n][k] += modV;
			return;
		}
		else if (x < -L)
		{
		//	cout << "x<-L\n\n";
			vec[0][k] += modV;
			//graficV[[0][k] += modV;
			return;
		}
		else if (fabs(x0) < 1e-13 && fabs(x) < dl)
		{
			++err;
			vec[this->x.n][k] += modV;
			return;
		}
		else if ( (int)(fabs(this->x.a - x) / dl) == (int)(fabs(this->x.a - x0) / dl) )
		{
			//cout << "else\n\n";
			//p = (int)(fabs(this->x.a - x) / dl);
			//vec[p][k] += Le / v;
			//graficV[[p][k] += modV;
			return;
		}
		
		//cout << "end\n\n";
		//if(x - x0 > 0 )
		//cout << "Le < dl\n";
		//p = (int)(fabs(this->x.a - x) / dl);
		//if (p==0) ++Kol[p][];
		
		
	}//*/
	

	if (x - x0 > 0) // случай когда летит с лева на право
	{

		if (x > 0) { x = 0., vec[this->x.n][k] += modV, end = false; p = this->x.n - 1; /*cout << "x > 0\n";*/ } //условия на концы 
		else
		{
			p = (int)(fabs(this->x.a - x) / dl);
		}

		l = (int)(fabs(this->x.a - x0) / dl) + 1;
	}
	else // случай когда летит с права  на лево
	{
		if (fabs(x0) < 1e-13)
		{
			p = this->x.n;
			++KOL;
		}
		else
		{
			p = (int)(fabs(this->x.a - x0) / dl);
		}

		if (x < -L) { x = -L, vec[0][k] += modV, end = false; l = 1;/*cout << "x < -L\n";*/ }
		else
		{
			l = (int)(fabs(this->x.a - x) / dl) +1;
		}	
	}
		
	for (; l <= p; l++) // пробегаемся по всех клеткам в которые атом пролетел полностью 
	{
		vec[l][k] += modV;	
	}
}
void nameFile(string& name, int index, int i, int j, int k)
{
	if (index < 10)
	{
		name[k] = (int)(48);
		name[i] = (int)(48);
		name[j] = (int)(48 + index);

	}
	else if (index < 100 && index > 9)
	{
		name[k] = (int)(48);
		name[i] = (int)(48 + index / 10);
		index = index - (index / 10) * 10;
		name[j] = (int)(48 + index);

	}
	else if (index > 99 && index < 1000)
	{
		name[k] = (int)(48 + index / 100);
		index = index - (index / 100) * 100;
		name[i] = (int)(48 + index / 10);
		index = index - (index / 10) * 10;
		name[j] = (int)(48 + index);
	}
}



/*
void СMatrix::distribution_function_vx(size_t N)
{
	double C = A() / (N * dx * dy * dy * dl);


	//for (int i = 0; i < x.n; ++i)
	for (int i = 0; i < vz.n; ++i)
	{

		//for (int j = 0; j < vx.n*vy.n*vz.n; ++j)
		for (int j = 0; j < vy.n; ++j)
		{
			for (int k = 0; k < vx.n; ++k)
			{
				for (int l = 0; l < x.n; ++l)
				{

					vec[i][j][k][l] *= C;
					//vec[i][j] *= C;
					//vec[i * vy.n + j][k * this->x.n + l] *= C;
					//vec[i][j][k][l]*= (C);
				}
			}
		}
	}
}*/

double A()
{
	return  (exp(-U_H * U_H) / sqrt(PI) - U_H * erfc(U_H)) / 2;
}




void OutPut_distribution_function_vx(vector<vector<double>>& grafic, СMatrix& M, int N, string &name, string& name1)
{
	//string name = (string)("data/dat2/000.txt");
	//char name1[30] = "data/dat2/0.txt";
	ofstream out;
	out.open(name1);
	cout << "name = " << name << endl;
	//      0 1             2                 3             4 5            6               7                8 9          10                 11   12          13
	out << M.vx << "\n" << M.vx.n << "\n" << M.dx << "\n" << M.x << "\n" << M.x.n << "\n" << M.dl << "\n" << M.vy << "\n" << M.vy.n << "\n" << M.vz << "\n" << M.vz.n << "\n";
	//    14          15             16            17              18           19              20                 21
	out << N << "\n" << L << "\n" << U_H << "\n" << U_p << "\n" << Cp << "\n" << CH;
	out.close();
	for (size_t i = 0; i < grafic.size(); ++i)
	{
		nameFile(name, i, 11, 12, 10);
		out.open(name);
		if (out.is_open())
		{
			for (size_t j = 0; j < grafic[i].size(); ++j)
			{
				out << grafic[i][j] << endl;
			}
		}
		out.close();
	}
}
/*
double СMatrix::integrate_vy_vz_Matrix(int k, int l)
{
	double integral = 0;
	int i, j;
	double D = dy * dy;
	//int kol = vy.n * vz.n;
	for (i = 0; i < vz.n; ++i)
	{
		for (j = 0; j < vy.n; ++j)
		{
			//integral += vec[i][j][k][l].t;
			integral += D * vec[i][j][k][l];
			//integral += D * vec[i * vy.n + j][k * this->x.n + l] ;
			//integral += D * vec[l][k * kol + (j * vz.n + i)];
		}
	}
	return  integral;

}

void СMatrix::matrix_vx_x_fH(vector<vector<double>>& Matrix)
{
	for (size_t i = 0; i < x.n; ++i)
	{
		for (size_t j = 0; j < vx.n; ++j)
		{
			Matrix[i][j] = integrate_vy_vz_Matrix(j, i);
		}
	}
}

*/


/*#include"distribution_function.h"
#include "include.h"
#include <omp.h>
namespace fs = std::experimental::filesystem;
ostream& operator<<(ostream& cout, Interval& P)
{
	cout << P.a <<"\n" << P.b;
	return cout;
}

bool СMatrix::check_beloning(Velociti_Atom& V)
{
	return (V.x > vx.a && V.x < vx.b) && (V.y > vy.a && V.y < vy.b )&&(V.z > vz.a && V.z < vz.b);
}
void СMatrix::filling_matrix(double x, double x0, double t, Velociti_Atom& V,double v)
{
	int l = 0, p = 0; //  индексы границы L
	double T =  dl/v; // время находжения в полной ячейке
	bool end = true; // для проверки где находится x 
	//int kol = vy.n * vz.n;
	int j;
	int index_x = (int)(fabs(this->vx.a - V.x) / dx); // высчитываем индексы по скоростям
	int index_y = (int)(fabs(this->vy.a - V.y) / dy);
	int index_z = (int)(fabs(this->vz.a - V.z) / dy);
	double tl = 0; // время кусочка клетки
	//cout << i << " " << j << " " << k << endl;
	if (x > 0) x = 0,  end = false; //условия на концы 
	if (x < -L) x = -L,  end = false;
	
	if (x - x0 > 0) // случай когда летит с лева на право
	{		
		p = (int)(fabs(this->x.a - x) / dl);	l = (int)(fabs(this->x.a - x0) / dl);	// высчитываем индексы по границы L	с учетом полета
		//cout <<"* " << p << " " << l << endl;
		l++;  tl = fabs(this->x.a + l * dl - x0) / v;		
		//vec[index_x * this->x.n + l - 1][index_z*vy.n +index_y]+= tl;
		//vec[l - 1 ][k*kol +(j*vz.n+i)] += tl;
		vec[index_z][index_y][index_x][l-1] += tl; p--;
		if (end)
		{
			tl = fabs(this->x.a + p * dl - x) / v;
			vec[index_z][index_y][index_x][p] += tl; p--;
			//vec[index_x * this->x.n + p][index_z * vy.n + index_y] += tl; p--;
			//vec[p][k * kol + (j * vz.n + index_x)] += tl; p--;
		}
	} 
	else // случай когда летит с права  на лево
	{	
		p = (int)(fabs(this->x.a - x0) / dl);  l = (int)(fabs(this->x.a - x) / dl);  // высчитываем индексы по границы L	с учетом полета
		//cout << "** " << p << " " << l << endl;
		if (x0 < 0)
		{
			tl = fabs(this->x.a + p * dl - x0) / v;
			vec[index_z][index_y][index_x][p] += tl;  p--;	
			//vec[index_x * this->x.n + p][index_z * vy.n + index_y] += tl;  p--;
			//vec[p][k * kol + (j * vz.n + index_x)] += tl; p--;
		}
		if (end)
		{
		
			l++;   tl = fabs(this->x.a + l * dl - x) / v;
			vec[index_z][index_y][index_x][l - 1] += tl;		
			//vec[index_x * this->x.n + l - 1][index_z * vy.n + index_y] += tl;
			//vec[l - 1][k * kol + (j * vz.n + i)] += tl;
		}
	}
	//cout <<"*** " <<  p <<" " << l << endl;
	
	
	// устанавливаем перемножаемые и результирующую матрицу общими для всех потоков;
	 //переменные i, j, k делаем частными, т.е. для каждого потока своя;
	 //устанавливаем статическое планирование распределения итераций по потокам (делим общее число итераций на потоки),
	 //каждый поток будет выполнять примерно одинаковое количество итераций
	for (; l <= p && l <this->x.n; l++) // пробегаемся по всех клеткам в которые атом пролетел полностью 
		vec[index_z][index_y][index_x][l] += T;	 // vec[index_x * this->x.n + l][index_z * vy.n + index_y]+= T;          ///vec[l][k * kol + (j * vz.n + i)] += T;  
}
void nameFile(string& name, int index, int i, int j,int k)
{
	if (index < 10)
	{
		name[k] = (int)(48);
		name[i] = (int)(48);
		name[j] = (int)(48 + index);
		
	}
	else if (index < 100 && index > 9)
	{
		name[k] = (int)(48);
		name[i] = (int)(48 + index / 10);
		index = index -  (index / 10) * 10;
		name[j] = (int)(48 + index);
		
	}
	else if (index > 99 && index < 1000)
	{
		name[k] = (int)(48 + index / 100);
		index = index - (index / 100) * 100;
		name[i] = (int)(48 + index/10);
		index = index -(index / 10) * 10;
		name[j] = (int)(48 + index);
	}
}




void СMatrix::distribution_function_vx(size_t N)
{
	double C = A() / (N*dx*dy*dy*dl);

	//#pragma omp parallel for shared(result, matrix_1, matrix_2) private(i, j, k) schedule(static, N / omp_get_num_threds()) 

	//for (int i = 0; i < x.n*vx.n; ++i)
	for (int i = 0; i < vz.n; ++i)
	{
		
		//for (int j = 0; j < vy.n*vz.n; ++j)
		for (int j = 0; j < vy.n; ++j)
		{
			for (int k = 0; k < vx.n; ++k)
			{
				for (int l = 0; l < x.n; ++l)
				{

					vec[i][j][k][l] *= C;
					//vec[i][j] *= C;
					//vec[i][k * this->vy.n + j] *= C;
					//vec[i][j][k][l]*= (C);
				}
			}
		}
	}
}

double A()
{
	return  ( exp(- U_H*U_H) / sqrt(PI) - U_H*erfc(U_H) ) /2;
}




void OutPut_distribution_function_vx(vector<vector<double>> & grafic, СMatrix& M,int N)
{
	string name = (string)("data/dat3/000.txt");
	char name1[30] = "data/dat3/0.txt";
	ofstream out;
	out.open(name1);
	cout << "name = " << name << endl;
	//      0 1             2                 3             4 5            6               7                8 9          10                 11   12          13
	out << M.vx << "\n" << M.vx.n <<"\n" << M.dx <<"\n" << M.x << "\n" << M.x.n << "\n" << M.dl << "\n" << M.vy << "\n" << M.vy.n << "\n" << M.vz << "\n" << M.vz.n << "\n";
	//    14          15             16            17              18           19    
	out <<N << "\n" << L << "\n" << U_H << "\n" << U_p << "\n" << Cp << "\n" << CH;
	out.close();
	for (size_t i = 0; i < M.x.n; ++i)
	{
		nameFile(name, i, 11, 12,10);
		out.open(name);
		if (out.is_open())
		{
			for (size_t j = 0; j < M.vx.n; ++j)
			{
				out << grafic[i][j] << endl;
			}
		}
		out.close();
	}
}
double СMatrix::integrate_vy_vz_Matrix(int k,int l)
{
	double integral = 0;
	long long int i, j;
	double D = dy * dy;
	int kol = vy.n * vz.n;
	for (i = 0; i < vz.n; ++i)
	{
		for (j = 0; j < vy.n; ++j)
		{		
				//integral += vec[i][j][k][l].t;
				integral += D*vec[i][j][k][l];
				//integral += D * vec[k * this->x.n + l][i * vy.n + j] ;
				//integral += D * vec[l][k * kol + (j * vz.n + i)];
		}
	}
	return  integral;

}

void СMatrix::matrix_vx_x_fH(vector<vector<double>>& Matrix)
{
	for (size_t i = 0; i < x.n; ++i)
	{
			for (size_t j = 0; j < vx.n; ++j)
			{
					//Matrix[j*x.n +i] = integrate_vy_vz_Matrix(j,i);	
					Matrix[i][j] = integrate_vy_vz_Matrix(j, i);
			}
	}
}




/*
double func(double a)
{
	return fabs(a)*exp(-(a - U_H) *(a - U_H))/sqrt(PI);
}
double integrate(double(*func)(double), double min_lim, double  max_lim, int n)
{
	double integral = 0;
	long long int i;
	double   step = (max_lim - min_lim) / n;
	double c, d;
	c = min_lim; d = c + step;
	for (i = 0; i < n; ++i)
	{
		integral += step * func((c + d) / 2);
		c = d; d += step;
	}
	return integral;
}
double method_of_rectangles(double(*func)(double), double min_lim, double max_lim, double delta)
{
	long long int n = 100;
	double sum1 = 0, sum2 = 0;

	sum2 = integrate(func, min_lim, max_lim, n);
	do
	{
		sum1 = sum2; sum2 = 0;
		n *= 2;
		sum2 = integrate(func, min_lim, max_lim, n);
		//cout << "sum1 = " << sum1 << "; sum2 = " << sum2 << endl;
	} while ((fabs(sum1 - sum2) > delta) && (n < N_));
	if (n > N_)
	{
		cout << "ERROR n > N\n";
	}
	return sum2;
	//cout << "Rectangles:" << n << " ; " <<sum2 << endl;
}
*/ 