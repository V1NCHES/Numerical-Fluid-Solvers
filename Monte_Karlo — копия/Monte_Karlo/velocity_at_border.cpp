#include"velocity_at_border.h"

ostream& operator<<(ostream& cout, Velociti_Atom& V)
{
	cout << V.x << " " << V.y << " " << V.z;
	return cout;
}

void Velociti_Atom::extractionem_velociti_at_border()
{

	double P4, Xi1 , z, h3, Xi4, Xi5, Xi6;
	do {
		P4 = calculatå_P4(sqrt(PI) * fabs(U_H));
		Xi1 = extractionem_0_1();
		z = calculatå_z(P4, Xi1);
		h3 = calculatå_h3(z);
		Xi4 = extractionem_0_1();
		//this->x = calculatå_V_x(z);
	} while (!(h3 > Xi4));
	this->x = calculatå_V_x(z);	

	Xi5 = sqrt(-log(extractionem_0_1())), Xi6 = 2 * PI * extractionem_0_1();
	
	this->y = calculatå_V_y(Xi5, Xi6);
	this->z = calculatå_V_z(Xi5, Xi6);
}
bool departure_check(Velociti_Atom& V, double length, double x)
{
	return true;
}


void Velociti_Atom::extractionem_velociti_P(Velociti_Atom& V_H1)
{
	double P4, Xi1, h1, Xi2, Xi3, Xi4, Xi5, Xi6;
	double w_x,  w_y,  w_z;
	Velociti_Atom u_p(U_p/Cp, 0, 0);
	//Velociti_Atom u_p(0.0, 0.0, 0.0);
	Velociti_Atom V, U;
	Velociti_Atom V_H(V_H1.x / Cp, V_H1.y / Cp, V_H1.z / Cp);
	//Velociti_Atom V_H(V_H1.x , V_H1.y , V_H1.z );

	double x = sqrt((V_H.x - u_p.x)* (V_H.x - u_p.x) + (V_H.y - u_p.y)* (V_H.y - u_p.y) + (V_H.z - u_p.z)* (V_H.z - u_p.z));
	//double x = sqrt(V_H.x* V_H.x + V_H.y * V_H.y + V_H.z * V_H.z);

	double C1, C2;
	double sigma_max = 1;
	do {
		P4 = calculatå_P4(sqrt(PI) * x / 2);
		Xi1 = extractionem_0_1();
		
		if (P4 < Xi1) 
		{			
			Xi2 = extractionem_0_1(), Xi3 = extractionem_0_1();
			Xi4 = extractionem_0_1(), Xi5 = 2 * PI * extractionem_0_1();
			C1 = sqrt(-log(Xi2 * Xi3));
			double omega_x = calculatå_omega_x(Xi4);
			C2 = sqrt(1 - omega_x * omega_x);
			w_x = C1* omega_x;
			w_y = C1 * calculatå_omega_y(C2, Xi5);
			w_z = C1 * calculatå_omega_z(C2, Xi5);
		}
		else {

			Xi2 = extractionem_0_1(), Xi3 = extractionem_0_1();
			Xi5 = 2 * PI * extractionem_0_1();
			Xi4 = sqrt(-log(extractionem_0_1()));
			w_x = calculatå_omega_y(sqrt(-log(Xi2)), PI * Xi3);
			w_y = calculatå_omega_y(Xi4, Xi5);
			w_z = calculatå_omega_z(Xi4, Xi5);
		}
		Velociti_Atom W(w_x, w_y , w_z);
		V = W + u_p;
		U = V_H - V;
		h1 = calculatå_h1(sigma_max, W, U, x);
		Xi6 = extractionem_0_1();
		//if (h1 > Xi6)
		//	cout << "h1 > Xi6\n";
	} while (!(h1 > Xi6 ) ); // && !(check_inf(V.x) && check_inf(V.y) && check_inf(V.z))
	V.x *= Cp, V.y *= Cp, V.z *= Cp;
	this->CopyOnly(V);
}


void Velociti_Atom::extractionem_velociti_reload_P(Velociti_Atom& V_H1)
{
	double P4, Xi1, h1, Xi2, Xi3, Xi4, Xi5, Xi6;
	double w_x, w_y, w_z;
	Velociti_Atom u_p(U_p / Cp, 0, 0);
	Velociti_Atom V, U;
	Velociti_Atom V_H(V_H1.x / Cp, V_H1.y / Cp, V_H1.z / Cp);

	//Velociti_Atom V_H_(V_H1.x *C_H, V_H1.y * C_H, V_H1.z * C_H);

//	double gx = (V_H1.x - U_p), gy = (V_H1.y), gz = (V_H1.z);
	//double x = sqrt(gx * gx + gy * gy + gz * gz)/Cp;

	double gx = (V_H1.x/Cp - U_p), gy = (V_H1.y/Cp), gz = (V_H1.z/Cp);
	double x = sqrt(gx * gx + gy * gy + gz * gz);
	double C1, C2;
	double sigma_max = sigma_(x);
	
	do {
		P4 = calculatå_P4(sqrt(PI) * x / 2);
		Xi1 = extractionem_0_1();

		if (P4 < Xi1)
		{
			Xi2 = extractionem_0_1(), Xi3 = extractionem_0_1();
			Xi4 = extractionem_0_1(), Xi5 = 2 * PI * extractionem_0_1();
			C1 = sqrt(-log(Xi2 * Xi3));
			double omega_x = calculatå_omega_x(Xi4);
			C2 = sqrt(1 - omega_x * omega_x);
			w_x = C1 * omega_x;
			w_y = C1 * calculatå_omega_y(C2, Xi5);
			w_z = C1 * calculatå_omega_z(C2, Xi5);
		}
		else {

			Xi2 = extractionem_0_1(), Xi3 = extractionem_0_1();
			Xi5 = 2 * PI * extractionem_0_1();
			Xi4 = sqrt(-log(extractionem_0_1()));
			w_x = calculatå_omega_y(sqrt(-log(Xi2)), PI * Xi3);
			w_y = calculatå_omega_y(Xi4, Xi5);
			w_z = calculatå_omega_z(Xi4, Xi5);
		}
		Velociti_Atom W(w_x, w_y, w_z);
		V = W + u_p;
		U = V_H - V;
		
		h1 = calculatå_h1_(sigma_max, W, U, x);
		Xi6 = extractionem_0_1();
		//if (h1 > Xi6)
		//	cout << "h1 > Xi6\n";
	} while (!(h1 > Xi6)); // && !(check_inf(V.x) && check_inf(V.y) && check_inf(V.z))
	V.x *= Cp, V.y *= Cp, V.z *= Cp;
	this->CopyOnly(V);
}

void OutPut(vector<double>& V, string& name)
{
	ofstream out;
	cout << "name = " << name << endl;
	out.open(name);
	if (out.is_open())
	{
		for (size_t i = 0; i < V.size(); ++i)
		{
			out << V[i] << "\n";
		}
	}
	out.close();
}

void OutPut(vector<long long int>& V, string& name)
{
	ofstream out;
	cout << "name = " << name << endl;
	out.open(name);
	if (out.is_open())
	{
		for (size_t i = 0; i < V.size(); ++i)
		{
			out << V[i] << "\n";
		}
	}
	out.close();
}

void OutPut(vector<Velociti_Atom>&vec, string& name)
{
	ofstream out;
	cout << "name = " << name << endl;
	out.open(name);
	if (out.is_open())
	{
		for (size_t i = 0; i < vec.size(); ++i)
		{
			out << vec[i] << " ";
		}
	}
	out.close();
}

double  calculatå_h3(double z)
{
	//double U_H_ = U_H / Cp;
	if (z < -U_H)
		return fabs(z + U_H) / (fabs(z) + fabs(U_H));
	return 0;
}

double calculatå_z(double P4, double Xi1)
{
	double Xi2 = extractionem_0_1();
	double Xi3 = extractionem_0_1();
	if (P4 > Xi1)
		return sqrt(-log(Xi2))* cos(PI * Xi3);

	if (Xi2 > 0.5)
		return sqrt(-log(2 * (1 - Xi2)));
	return -sqrt(-log(2 * Xi2));
	

}

double  calculatå_P4(double C)
{
	return (C / (1 + C));
}

double  calculate_P5(double P4)
{
	return 1.0 - P4;
}

double  calculatå_V_x(double z)
{
	return z + U_H;
}

double  calculatå_V_y(double Xi5, double Xi6)
{
	return Xi5 * cos(Xi6);
}

double  calculatå_V_z(double Xi5, double Xi6)
{
	return Xi5 * sin(Xi6);
}

double  calculatå_omega_x(double Xi4)
{
	return 1 - 2*Xi4;
}
double  calculatå_omega_y(double omega_x, double Xi5)
{
	return omega_x * cos(Xi5);
}
double  calculatå_omega_z(double omega_x, double Xi5)
{
	return omega_x * sin(Xi5);
}

double  calculatå_h1(double sigma_max, Velociti_Atom& W, Velociti_Atom& U, double x)
{

	return U.calculete_length() * sigma_max / (sigma_max * (x + W.calculete_length()));
}


double  calculatå_h1_(double sigma_max, Velociti_Atom& W, Velociti_Atom& U, double x)
{
	//double l = U.calculete_length();  U.multiplicationby_scalar(C_H);
	return U.calculete_length() * sigma_(U.calculete_length())  / (sigma_max * (x + W.calculete_length()));
}