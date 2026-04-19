#include"include.h"
#include"velocity_at_border.h"
#include "particle_collision.h"
#include "free_path_length.h"
#include"distribution_function.h"
#include "flight_atom.h"
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <random>
#include <string>

double min_(double& x, double& y)
{
	if (x < y)
		return x;
	return y;
}
double max_(double& x, double& y)
{
	if (x > y)
		return x;
	return y;
}
int main()
{
	srand(time(0));
	
	double startTime, endTime;
	
	startTime = getCPUTime();
	СMatrix M;
	endTime = getCPUTime();
	cout<<"Created M: CPU time:"<< endTime - startTime <<"\n";
	M.printN();

	vector<vector<double>> graficV;
	//for (size_t i = 0; i <= M.x.n; ++i)
	//	graficV.push_back(vector<double>(M.vx.n, 0));
	//cout << "M.x.n = " << M.x.n << "\n";

	size_t N = 1e6;// количество атомов на границе
	cout << "N = " << N << endl; //3/(M.d* M.d* M.d *M.dl

	int ptr = 0;
	startTime = getCPUTime();
	flight_atom(M, ptr, N, graficV);
	endTime = getCPUTime();
	cout << "\nM.filling_matrix: CPU time:" << endTime - startTime << "\n";
	
	//startTime = getCPUTime();
	//M.distribution_function_vx(N);
	//endTime = getCPUTime();
	//cout << "\nM.distribution_function_vx: CPU time:" << endTime - startTime << "\n";// 

	cout << "A() = " << A() << endl;

	//startTime = getCPUTime();
	//vector<vector<double>> grafic;
	//for (size_t i = 0; i < M.x.n; ++i)
	//	grafic.push_back(vector<double>(M.vx.n,0));
	//endTime = getCPUTime();
	//cout << "\nGrafic: CPU time:" << endTime - startTime << "\n";


	//startTime = getCPUTime();
	//M.matrix_vx_x_fH(grafic);
	//endTime = getCPUTime();
	//cout << "\nM.matrix_vx_x_fH:  CPU time:" << endTime - startTime << "\n"; //
	
	cout << "ptr = " << ptr <<" N /ptr =  "<< 100 *(double)ptr  / N <<"%" << endl;
	//
	//OutPut(M);
	//OutPut_distribution_function_vx(M);
	string name = (string)("data/dat5/000.txt"), name1 = (string)"data/dat5/0.txt";

	startTime = getCPUTime();
	OutPut_distribution_function_vx(M.vec, M,N, name, name1);
	endTime = getCPUTime();
	cout << "\nOutPut_distribution_function_vx:  CPU time:" << endTime - startTime << "\n";

	cout << "Ve = " << Ve <<"C_H = "<<C_H << "\n";
	//string name2 = (string)("dat1/dat3/000.txt");
	//startTime = getCPUTime();
	//OutPut_distribution_function_vx(graficV, M, N, name2, name1);
	//endTime = getCPUTime();
	//cout << "\nOutPut_distribution_function_vx:  CPU time:" << endTime - startTime << "\n";
	//OutPut(FH,M);
	//*/
	//for (size_t i = 0; i < M.kolL.size(); ++i)
	//	cout << M.kolL[i] << endl;
	//*/
	//distagramma_normal_extractionem();
	
	//distagramma_velociti_at_border();
	// distagramma_free_path_length_hard_balls();
	//distagramma_velociti_proton();

   // distagramma_free_path_length_Reloads();

	//cout << "C_H = " << C_H << "\n";
	//distagramma_new_velociti();
	// 
	//for (int i = 0; i < 100; ++i)
	//	cout << acos(extractionem_1_1()) << endl;
	//vector <double> min_lim = { -8.,-5. , -5.};
	//vector <double> max_lim = { 8.,5. , 5. };
	//cout << method_of_rectangles_R3(f_P, min_lim, max_lim, 1e-3, -2.5, 1.25, -0.9) << endl;
	//cout << A();
	/*
	int a1 = 1,  a2 = 1,  a3 =1 ;
	long long int n = 0;
	double a = 0., b = 1., d = 0.01;
	int size_ = (int)((b - a) / d);
	vector <long long int> distagramma(size_, 0);
	double cur = 0.;
	for (int i = 0; i < 1e8; ++i)
	{
		cur = normal_extractionem_0_1(a1, a2, a3);
		if (cur > a && cur < b)
		{
			++distagramma[(int)(fabs(a - cur) / d)], n++;
		}
		else
		{
			cout << cur << "\n";
		}
	}
	string name = "data/0.txt";
	OutPut(distagramma, name);
	//*/
	cout << "end\n";
	return 0;
}