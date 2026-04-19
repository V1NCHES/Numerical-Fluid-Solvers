#include "flight_atom.h"

void distagramma_normal_extractionem()
{
	cout << "distagramma_normal_extractionem\n";
	string name = "check1/chek/000.txt";
	size_t N = 1e8;
	long long int n = 0;
	double Gamma;
	double a = 0, b = 1, d = 0.01;
	int a1 = 0, a2 = 1, a3 = 0;
	int size_ = (int)((b - a) / d);
	vector <long long int> distagramma_Vx(size_, 0);
	for (size_t i = 0; i < N; ++i)
	{
		Gamma = normal_extractionem_0_1(a1,a2,a3);
		if (Gamma > a && Gamma < b)
			++distagramma_Vx[(int)(fabs(a - Gamma) / d)], n++;
	}
	vector<double> check_ = { a,b,d };
	distagramma_Vx.push_back(n);
	OutPut(distagramma_Vx, name);
	name = "check1/chek/000_.txt";
	OutPut(check_, name);
}

void distagramma_velociti_at_border()
{
	cout << "distagramma_velociti_at_border\n";
	string name = "check1/000.txt";
	size_t N = 1e8;
	long long int n = 0;
	Velociti_Atom VH;
	double a = -7, b = 0, d = 0.025;
	int size_ = (int)((b - a) / d);
	vector <long long int> distagramma_Vx(size_,0);
	for (size_t i = 0; i < N; ++i)
	{
		VH.extractionem_velociti_at_border();
		if(VH.x > a && VH.x < b)
			++distagramma_Vx[(int)(fabs(a - VH.x) / d)], n++;
	}
	vector<double> check_ = {a,b,d , U_H, Cp};
	distagramma_Vx.push_back(n);
	OutPut(distagramma_Vx, name);
	name = "check1/000_.txt";
	OutPut(check_, name);
}

void distagramma_free_path_length_hard_balls()
{
	cout << "distagramma_free_path_length_hard_balls\n";
	Velociti_Atom VH, W_p(U_p/Cp, 0, 0);
	VH.extractionem_velociti_at_border();
	string name = "check1/0011.txt";
	size_t N = 1e8;
	long long int n = 0;
	double length, U_, nu, gamma;
	double a = 0, b = 10, d = 0.05;
	int size_ = (int)((b - a) / d);
	vector <long long int> distagramma_free_path_length(size_, 0);	
	//double g_x = (VH.x - W_p.x) / Cp, g_y = (VH.y - W_p.y) / Cp, g_z = (VH.z - W_p.z) / Cp;
	VH.set_V(VH.x/Cp, VH.y/Cp, VH.z/Cp);
	double g_x = (VH.x - W_p.x) , g_y = (VH.y - W_p.y) , g_z = (VH.z - W_p.z);

	double x = sqrt(g_x * g_x + g_y * g_y + g_z * g_z);
	double  modelVH = VH.calculete_length();
	for (size_t i = 0; i < N; ++i)
	{		
		gamma = extractionem_0_1();	
		U_ = exp(-x * x) / sqrt(PI) + (x + 1 / (2 * x)) * erf(x);
		nu = U_ *CH/ Cp;
		if (i % 10000000 == 0) cout << i << endl;
		length = -modelVH * log(gamma) / nu;
		if (length > a && length < b)
			++distagramma_free_path_length[(int)(fabs(a - length) / d)], n++;
	}
	cout << "n = " << n << endl;	
	vector<double> check_ = { a,b,d , VH.x , VH.y , VH.z , nu / VH.calculete_length(), Cp };
	distagramma_free_path_length.push_back(n);
	OutPut(distagramma_free_path_length, name);
	name = "check1/0011_.txt";
	OutPut(check_, name);
}
double sigma_(double v)
{
	return pow(a1 - a2 * log(v),2);
}
void distagramma_free_path_length_Reloads()
{
	cout << "distagramma_free_path_length_Reloads\n";
	Velociti_Atom VH(-3.1,-0.4, 0.4), W_p(U_p, 0, 0);
	//VH.extractionem_velociti_at_border();
	string name = "data/1.txt";
	size_t N = 1e8;
	long long int n = 0;
	double length, U_, nu, gamma;
	double a = 0, b = 10, d = 0.05;
	int size_ = (int)((b - a) / d);
	vector <long long int> distagramma_free_path_length(size_, 0);
	/*
	const double k = 1.38064852 * 1e-16;//  протона в граммах
	const double mp = 1.6735575e-24; // г
	const double T_H = 6530; //K
	const double k_mh = 13806.49 / 1.6735575;
	const double a1 = 2.2835 * 1e-7;
	const double a2 = 1.062 * 1e-8;
	CH = 1.0
	Cp = sqrt(2 * k * T_H / mp)/1e5;
*/
	double E = 1.6 * 1e-12;
	double sigma = 0, w = sqrt(4*E/mp) ;
	double sqrt_sigma_0 = sigma_(w);
	double x = 0.;
	double cp = sqrt(2 * k * T_H / mp);

	for (int i = 0; i < 1e7; ++i)
	{

		VH.extractionem_velociti_at_border();
		b = 1 - a2 * log(VH.calculete_length()*C_H/w)/ sigma_Ve;
		sigma = b*b;	
		x = sqrt((VH.x - W_p.x) * (VH.x - W_p.x) + (VH.y - W_p.y) * (VH.y - W_p.y) + (VH.z - W_p.z) * (VH.z - W_p.z)) /Cp;
		U_ = exp(-x * x) / sqrt(PI) + (x + 1 / (2 * x)) * erf(x);
		nu = U_ * sigma * Cp;
		gamma = extractionem_0_1();
		length = -VH.calculete_length() * log(gamma) / nu;
		if (length > a && length < b)
			++distagramma_free_path_length[(int)(fabs(a - length) / d)], n++;
		//nu = 
		//cout << "length = " << length <<  "\n\n";
	}//*/
	OutPut(distagramma_free_path_length, name);
	return;
}//*/


void distagramma_velociti_proton()
{
	cout << "distagramma_velociti_proton\n";
	Velociti_Atom VH(-2.5, 1.25, -0.9);
	string name = "check1/0022.txt";
	size_t N = 1e8;

	long long int n = 0;
	double a = -7, b = 1, d = 0.01;
	int size_ = (int)((b - a) / d);
	vector <long long int> distagramma_W(size_, 0);

	double P4, Xi1, h1, Xi2, Xi3, Xi4, Xi5, Xi6;
	double w_x, w_y, w_z;
	Velociti_Atom u_p(U_p, 0, 0);
	Velociti_Atom W, U, V;
	VH.set_V(VH.x , VH.y, VH.z );
	//Velociti_Atom V_H(V_H1.x, V_H1.y, V_H1.z);
	double x = sqrt((VH.x - u_p.x) * (VH.x - u_p.x) + (VH.y - u_p.y) * (VH.y - u_p.y) + (VH.z - u_p.z) * (VH.z - u_p.z));
	double C1, C2;
	double sigma_max = 1;
	P4 = calculatе_P4(sqrt(PI) * x / 2);
	for (size_t i = 0; i < N; ++i)
	{ 
		do {
			
			Xi1 = extractionem_0_1();

			if (P4 < Xi1)
			{
				Xi2 = extractionem_0_1(), Xi3 = extractionem_0_1();
				Xi4 = extractionem_0_1(), Xi5 = 2 * PI * extractionem_0_1();
				C1 = sqrt(-log(Xi2 * Xi3));
				C2 = sqrt(1 - calculatе_omega_x(Xi4) * calculatе_omega_x(Xi4));
				w_x = C1 * calculatе_omega_x(Xi4);
				w_y = C1 * calculatе_omega_y(C2, Xi5);
				w_z = C1 * calculatе_omega_z(C2, Xi5);
			}
			else {
				Xi2 = extractionem_0_1(), Xi3 = extractionem_0_1();
				Xi4 = extractionem_0_1(), Xi5 = 2 * PI * extractionem_0_1();
				Xi4 = sqrt(-log(Xi4));
				w_x = calculatе_omega_y(sqrt(-log(Xi2)), PI * Xi3);
				w_y = calculatе_omega_y(Xi4, Xi5);
				w_z = calculatе_omega_z(Xi4, Xi5);
			}
			W.set_V(w_x, w_y, w_z);
			V = W + u_p;
			U = VH - V;
			h1 = calculatе_h1(sigma_max, W, U, x);
			Xi6 = extractionem_0_1();
			//if (h1 > 1)
			//	cout << "h1 > 1\n";
		} while (!(h1 > Xi6));

		if (i % 1000000 == 0) cout << i << endl;
		if (V.x > a && V.x < b)
			++distagramma_W[(int)(fabs(a - V.x) / d)], n++;
	}

	cout << " n = " << n << endl;
	vector<double> check_ = { a,b,d , VH.x, VH.y, VH.z, U_p, Cp, exp(-x * x) / sqrt(PI) + (x + 1 / (2 * x)) * erf(x) };
	distagramma_W.push_back(n);
	OutPut(distagramma_W, name);
	name = "check1/0022_.txt";
	OutPut(check_, name);
}

void distagramma_new_velociti()
{
	cout << "distagramma_new_velociti\n";
	Velociti_Atom VH(-3, 0.75, -0.8);
	string name = "check1/003.txt";
	size_t N = 1e7;
	Velociti_Atom VH_new, V;
	long long int n = 0, nxi =0;
	double a = -7, b = 1, d = 0.01;
	int size_ = (int)((b - a) / d);
	vector <long long int> distagramma_W(size_, 0);
	//a = 0, b =PI, d = 0.01;
	//size_ = (int)((PI) / d) +1;
	//vector <long long int> distagramma_Xi(size_, 0);

	//a = 0, b = 2 * PI, d = 0.01;
	//size_ = (int)((2*PI) / d) +1;
	//vector <long long int> distagramma_eps(size_, 0);
	//double Xi , eps_, l,C1,C2,C3;
	//Velociti_Atom g,e1,e ,e2,e3 , g_new, g_new_;
	//for (size_t i = 0; i < N; ++i)
	{
		V.extractionem_velociti_P(VH);
		VH_new = calculete_new_velocity(VH, V);
		/*
		Xi = extractionem_angle_xi(), eps_ = extractionem_angle_eps();
		g.set_V(VH.x - (VH.x + V.x) / 2, VH.y - (VH.y + V.y) / 2, VH.z - (VH.z + V.z) / 2);

		//Velociti_Atom g((V.x - W.x)/2 ,  (V.y - W.y)/2 , (V.z  - W.z)/2);
		l = g.calculete_length();
		C1 = cos(Xi), C2 = sin(Xi) * cos(eps_), C3 = sin(Xi) * sin(eps_);

		e1.set_V(g.x / l, g.y / l, g.z / l), e.set_V(1, 0, 0);
		if (is_equal(g.y, 0) && is_equal(g.z, 0))
			e.set_V(0, 1, 0);
		
		e2 = e.vector_product(e1);
		e2.normalization();
		//Le3 = e2.calculete_length();
		//e2.set_V(e2.x / Le3, e2.y / Le3, e2.z / Le3);
		e3 = e1.vector_product(e2);
		

		//Le3 = e3.calculete_length();
		//e3.set_V(e3.x / Le3, e3.y / Le3, e3.z / Le3);

	//	if ( fabs(e1.calculete_length() - 1) > 1e-5)
	//		cout << e1.calculete_length() <<"    fabs(e1.calculete_length() - 1)>1e-5\n";
		if (fabs(e2.calculete_length() - 1) > 1e-5)
			cout << e2.calculete_length()<< "    fabs(e2.calculete_length() - 1) > 1e-5\n";
		if (fabs(e3.calculete_length() - 1) > 1e-5)
			cout << e3.calculete_length()<<"    fabs(e3.calculete_length() - 1) > 1e-5\n";
		//if (!is_equal(e1.scalar_product(e2), 0))
		//	cout << "e1.scalar_product(e2)\n";
		
		g_new.set_V( l*(e1.x * C1 + e2.x * C2 + e3.x * C3), l*(e1.y * C1 + e2.y * C2 + e3.y * C3), l*(e1.z * C1 + e2.z * C2 + e3.z * C3));

		//if (fabs(g_new.calculete_length() -l) > 1e-5)
		//	cout << g_new.calculete_length() << "    fabs(g_new.calculete_length() -l >1e-5\n";
		//cout << "VH = " << VH << " ; V =" <<V <<" ; g = "<<g<<" ; l = "<<l << " ; g_new = " << g_new << endl;
		g_new_.set_V(e1.x * g_new.x + e2.x * g_new.y + e3.x * g_new.z, e1.y * g_new.x + e2.y * g_new.y + e3.y * g_new.z, e1.z * g_new.x + e2.z * g_new.y + e3.z * g_new.z);
		//cout << "e1 = " << e1 << " ; e2 = " << e2 << " ; e3 = " << e3 <<" ; g_new_ ="<< g_new_ << endl<<endl;
		//g_new_.set_V(e1.x * g_new.x + e2.x * g_new.y + e3.x * g_new.z, e1.y * g_new.x + e2.y * g_new.y + e3.y * g_new.z, e1.z * g_new.x + e2.z * g_new.y + e3.z * g_new.z);

		//VH_new.set_V(g_new_.x + (VH.x + V.x) / 2, g_new_.y + (VH.y + V.y) / 2, g_new_.z + (VH.z + V.z) / 2);
		VH_new.set_V(g_new.x + (VH.x + V.x) / 2, g_new.y + (VH.y + V.y) / 2, g_new.z + (VH.z + V.z) / 2);
		*/
		//if (i % 1000000 == 0) cout << i << endl;

		//if (VH_new.x > a && VH_new.x < b)
		//	++distagramma_W[(int)(fabs(a - VH_new.x) / d)], n++;

		//if (Xi > 0 && Xi < PI)
		//	++distagramma_Xi[(int)(Xi / d)], nxi++;

		//if (eps_ > 0 && eps_ < 2* PI)
		//	++distagramma_eps[(int)(fabs(- eps_) / d)], neps++;
	}
	cout << "VH = "<< VH<<"V= "<<V <<"VH_new = " <<VH_new;
	/*cout << " n = " << n << endl;
	vector<double> check_ = { a,b,d , VH.x, VH.y, VH.z, U_p, Cp};
	distagramma_W.push_back(n);
	OutPut(distagramma_W, name);
	name = "check1/003_.txt";
	OutPut(check_, name);
	*/
	//distagramma_Xi.push_back(nxi);
	//name = "check1/004.txt";
	//OutPut(distagramma_Xi, name);
	//distagramma_eps.push_back(neps);
	//name = "check1/005.txt";
	//OutPut(distagramma_eps, name);
}


void flight_atom(СMatrix& M, int& ptr, int N,vector<vector<double>>& graficV)
{
	//Velociti_Atom W, V, W_p(U_p, 0, 0); //скорость протона при столкновении, и скрость нужна для подсчета длины свободного пробега
	Velociti_Atom W, V, W_p(U_p, 0, 0); //скорость протона при столкновении, и скрость нужна для подсчета длины свободного пробега

	//int ptr = 0;
	//vector<int> Error(M.x.n, 0);
	//vector<int>Kol(M.x.n, 0);
	int err = 0, kol = 0;
	///vector<vector<int>> Kol;
	//for (size_t i = 0; i < M.x.n; ++i)
	//	Kol.push_back(vector<int>(M.x.n, 0));
	
///	vector<int> Kol;
	double modelV = 0;
	int k = 0, j = 0, l = 0; // счетчики колчества улетевших вправо и количество столкновений
	//double dL =  0.01; //размер ячейки, препологаю что она настолько мала что в ней не могут столкнуться,	
	double length = 0; // длина свободного пробега для тукущего частицы i в текущий момент времени t
	double t = 0, x0 = 0, x; // время, и координата текущей частицы
	//ofstream out;
	//out.open("data/OutPut_Error.txt" );

	for (size_t i = 0; i < (size_t)N; ++i)
	{
		if (i % 1000000 == 0)
			cout << i << endl;
		//cout << "!!!!!!!!!\n";
		do{ V.extractionem_velociti_at_border();} while (!M.check_beloning(V));

		t = 0, x0 = 0; x = 0; //начальные данные: время, текущее положение, положение в момент столкновения	
	//	cout << "------------------------------\n";
		do
		{
		//	cout << "^^^^^^\n";
			/*if (fabs(V.x)<1e-2)
			{
				//cout << "ERROR 1 nan(inf)\n";
				break;
			}*/

		    length = hard_balls(V, W_p); // высчитываем длину свободного пробега до столкновения для упругих шаров
			// 
			//length = reload(V, W_p); // высчитываем длину свободного пробега до столкновения для перезарядки
			//if (length < M.dl) cout << "length < M.dl\n";


			modelV = V.calculete_length();
			//if (is_equal(0, modelV))
			//{
			//	cout << "ERROR 1 nan(inf)\n";
			//	break;
			//}
			 
			t = length / modelV; // считаем время полета до столкновения 
			x = x0 + V.x * t; // высчитываем координату столкновения
			//cout << "V = " << V <<"; length = "<< length<<"; t = "<<t << "; x = " << x << "; x0 = " << x0 << endl;
			if (M.check_beloning(V))
			{
				//if (length < M.dl) 
			//	cout << "M.check_beloning(V)\n";

				M.filling_matrix(x, x0, t, V, fabs(V.x),N,graficV, length, err, kol);
				x0 = x;
				//cout << "vec[i] = " << vec[i] << "\nt = " << t << " ; length = " << length << " ; x = " << x << endl;
				if (x > 0) { k++; break; }
				if (x > -L )
				{
					
					W.extractionem_velociti_P(V); //разыгрываем скорость протона при столновении с атомом водорода
					//W.extractionem_velociti_reload_P(V); //разыгрываем скорость протона при перезарядке с атомом водорода
					V = calculete_new_velocity(V, W); j++;// вычесляем новую скорость после упругого столкновения	
					//V = W; 
					j++; 
				}
				else { l++;  break; }
			}
			else
			{
				 break;//cout << "check_beloning(V)" << vec[i] << "\n"; 
			}
		} while (true); //проверяем пока шар не вылетит из границы
	}
	cout << "KOL = " << kol <<"; Err = "<<err<<" sum = "<<kol+err << "\n";
	cout << "N = " << N << "\npoints flying from the right = " << k << "\nnumber of collisions = " << j << endl;
	cout << "points flying from the left = " << l << endl;
	//cout << "number of discarded trajectories = " << ptr << endl;
	//out.close();
	
	cout << "x = " << -L << "; kol_Error = " << err <<"; kol = "<<kol << "\n";

}
double G(double wx, double wy, double wz, double vx, double vy, double vz)
{
	return sqrt((vx-wx)* (vx - wx) + (vy - wy) * (vy - wy) + (vz - wz) * (vz - wz));
}
double f_P(double wx, double wy, double wz, double vx, double vy, double vz, double U)
{
	return G(wx, wy, wz, vx, vy, vz) * exp(- ((wx - U)*(wx - U) - wy*wy - wz*wz ) / (Cp *Cp) );
}



double method_of_rectangles_R3(double(*f)(double, double, double, double, double, double, double), vector <double> min_lim, vector <double> max_lim, double delta, double v_x, double v_y, double v_z)
{
	long long int n = 1000;
	double sum1 = 0, sum2 = 0;

	sum2 = integrate_R3(f, min_lim, max_lim, n, v_x, v_y, v_z); //вычесляем интеграл для разбиения n
	do
	{
		sum1 = sum2; sum2 = 0;   n *= 2;
		sum2 = integrate_R3(f, min_lim, max_lim, n, v_x, v_y, v_z); //вычесляем интеграл для разбиения n*2

	    cout << "sum1 = " << sum1 << "; sum2 = " << sum2 << " ;  n*2 = " << n << endl;
	} while ((fabs(sum1 - sum2) > delta) && (n < 1e5));// поверка на разность сумм меньше eps  + для того чтобы была сходимость ограничение по разбиению

	if (n > 1e5)
	{
		cout << "ERROR n > N\n";
	}

	// cout << "Rectangles:" << n << " ;\n";// << sum2 << endl;
	return sum2;
}

double integrate_R3(double(*f)(double, double, double, double, double, double, double), vector <double> min_lim, vector <double> max_lim, int n, double v_x, double v_y, double v_z)
{
	double integral = 0;
	long long int i, j, k;
	vector <double> delte; // вектор разбиения отрезков, то есть  сдив по каждой координате запоминаю для текущего разбиения n
	for (size_t p = 0; p < min_lim.size(); ++p)
		delte.push_back((double)(max_lim[p] - min_lim[p]) / n), cout<< delte[p]<<endl;

	double V = delte[0] * delte[1] * delte[2]; // вычесляем обЪем клетки разбиения 

	 // текущие x_i и x_i+1
	double cx = min_lim[0], cy = min_lim[1], cz = min_lim[2], dx = min_lim[0] + delte[0], dy = min_lim[1] + delte[1], dz = min_lim[2] + delte[2];


	for (i = 0; i < n; ++i)
	{
		cy = min_lim[1]; dy = min_lim[1] + delte[1];
		for (j = 0; j < n; ++j)
		{
			cz = min_lim[2]; dz = min_lim[2] + delte[2];
			for (k = 0; k < n; ++k)
			{
				//integral += delte[0] * delte[1] * delte[2] * f((cx + dx) / 2, (cy + dy) / 2, (cz + dz) / 2, v_x, v_y, v_z, x, U_H);
				integral +=  f((cx + dx) / 2, (cy + dy) / 2, (cz + dz) / 2, v_x, v_y, v_z, U_p);
				cz = dz;  dz += delte[2]; //c[1] = d[1]; d[1] += delte[1];
			}
			cz = 0;  dz = 0;
			cy = dy; dy += delte[1];
		}
		cx = min_lim[0], dx += delte[0];
		cy = 0; dy = 0;
	}
	return V * integral;
}