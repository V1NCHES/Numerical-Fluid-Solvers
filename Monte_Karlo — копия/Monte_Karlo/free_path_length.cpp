#include"free_path_length.h"


double reload(Velociti_Atom& V, Velociti_Atom& W)
{
	//Velociti_Atom VH(V.x*C_H , V.y * C_H ,  V.z * C_H);

	double sigma = pow(1 - a2 * log(V.calculete_length()/Ve)/sigma_Ve, 2);

	//double gx = (V.x - W.x), gy = (V.y - W.y), gz = (V.z - W.z);
	//	double x = sqrt(gx * gx + gy * gy + gz * gz) / Cp;

	double gx = (V.x / Cp - W.x), gy = (V.y / Cp - W.y), gz = (V.z /Cp - W.z);
	double x = sqrt(gx * gx + gy * gy + gz * gz) ;

	double U_ = exp(-x * x) / sqrt(PI) + (x + 1. / (2 * x)) * erf(x);
	double nu = U_ * Cp * sigma;
	double gamma = extractionem_0_1();
	
	return (-V.calculete_length() * log(gamma) / nu);
}

double hard_balls(Velociti_Atom &V, Velociti_Atom &W)
{
	//double sigma = 1;
	//Velociti_Atom VH(V.x , V.y , V.z );
	//double gx = V.x - W.x , gy = V.y  - W.y, gz = V.z  - W.z;
	//Velociti_Atom VH(V.x/Cp, V.y/Cp, V.z/Cp);
	double gx = (V.x - W.x), gy = (V.y - W.y) , gz = (V.z - W.z) ;

	//double gx = (V.x - W.x ) / Cp, gy = (V.y - W.y ) / Cp, gz = (V.z - W.z) / Cp;
	double x = sqrt(gx * gx + gy * gy + gz * gz)/Cp ;
	//double gx = V.x/Cp, gy = V.y / Cp , gz = V.z / Cp;
	//double x = sqrt(gx * gx + gy * gy + gz * gz);
	double U_ = exp(-x * x)/sqrt(PI)  + (x+ 1./(2*x)) * erf(x);
	double nu = U_ * Cp / CH;
	double gamma = extractionem_0_1();
	// VH(V.x / Cp, V.y / Cp, V.z / Cp);
	return (-V.calculete_length()*log(gamma) /nu);
}