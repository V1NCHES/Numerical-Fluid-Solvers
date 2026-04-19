#ifndef velocity_at_border_h
#define velocity_at_border_h
#include"include.h"
class Velociti_Atom;

double  calculatå_h3(double z);
double  calculatå_z(double P4, double Xi1);
double  calculatå_P4(double C);
double  calculate_P5(double P4);
double  calculatå_V_x(double z);
double  calculatå_V_y(double Xi5, double Xi6);
double  calculatå_V_z(double Xi5, double Xi6);


double  calculatå_h1(double sigma_max, Velociti_Atom& W, Velociti_Atom& U, double x);
double  calculatå_h1_(double sigma_max, Velociti_Atom& W, Velociti_Atom& U, double x);


double  calculatå_omega_x(double Xi4);
double  calculatå_omega_y(double omega_x, double Xi5);
double  calculatå_omega_z(double omega_x, double Xi5);
bool departure_check(Velociti_Atom& V, double length, double x);
class Velociti_Atom
{
public:
	double x, y, z;

	Velociti_Atom(){ extractionem_velociti_at_border();} //ðîçûãðûø ñêîðîñòè íà ãðàíèöå

	Velociti_Atom(double x1, double y1, double z1){ x = x1; y = y1; z = z1; }

	double calculete_length(){ return sqrt(x * x + y * y + z * z);}


	double calculete_distances(Velociti_Atom &V)
	{
		double vx = x - V.x, vy = y - V.y, vz = z - V.z;
		return sqrt(vx * vx + vy * vy + vz * vz);
	}
	void extractionem_velociti_at_border();

	void extractionem_velociti_P(Velociti_Atom& V_H1);
	void extractionem_velociti_reload_P(Velociti_Atom& V_H1);
	void set_V(double x1,double y1, double z1){ x = x1; y = y1; z = z1; }
	void normalization() { double l = calculete_length(); x /= l; y /= l; z /= l;}
	Velociti_Atom multiplicationby_scalar(double lamda) { return Velociti_Atom(x * lamda, y * lamda, z * lamda);}
	Velociti_Atom vector_product(Velociti_Atom& V ){ return Velociti_Atom( y*V.z - z*V.y , -x * V.z + z * V.x,  x * V.y - y * V.x ); }
	double scalar_product(Velociti_Atom& V) { return V.x *x + V.y * y+ V.z *z; }
	Velociti_Atom operator+(Velociti_Atom &V){ return Velociti_Atom(V.x + x, V.y + y, V.z + z);}
	void multiplicationby_scalar_I(double lamda) { x *= lamda, y *= lamda, z *= lamda; }
	Velociti_Atom operator-(Velociti_Atom& V) { return Velociti_Atom(x - V.x,y - V.y, z - V.z); }
	//void OutPut();
	void CopyOnly(Velociti_Atom& V) { x = V.x; y = V.y; z = V.z; }
	friend ostream& operator<<(ostream& cout, Velociti_Atom& P);//îïðàòîð âûâîäà
};

//Velociti_Atom operator+(Velociti_Atom& V , Velociti_Atom& W) { return Velociti_Atom(V.x + W.x, V.y + W.y, V.z + W.z); }
void OutPut(vector< Velociti_Atom>& vec, string& name);
void OutPut(vector<double>& V, string &  name);
void OutPut(vector<long long int>& V, string& name);
#endif