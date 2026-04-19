#include"particle_collision.h"
double extractionem_angle_xi()
{
	return acos(extractionem_1_1()); //gamma in (-1 ; 1)
}

double extractionem_angle_eps()
{
	return 2*PI * extractionem_0_1(); //gamma in (0 ; 1)
}

Velociti_Atom calculete_new_velocity(Velociti_Atom& V, Velociti_Atom& W)
{
	//double Xi = PI, eps_ = extractionem_angle_eps(); // розыгрыш углов //extractionem_angle_xi()
	double Xi = extractionem_angle_xi(), eps_ = extractionem_angle_eps();
	Velociti_Atom g(V.x - (V.x + W.x) / 2, V.y - (V.y + W.y) / 2, V.z - (V.z + W.z) / 2); // высчитываем g - в системе ЦМ
	//cout << "eps_ = " << eps_ << "\n";
	double l = g.calculete_length(); // считаем длину g 
	double C1 = cos(Xi), C2 = sin(Xi) * cos(eps_),  C3 = sin(Xi) * sin(eps_); // считаем константя для векторов e1, e2, e3
	//cout << "C1 = " << C1 << "; C2 = " << C2 << "; C3 = " << C3 << "\n";
	Velociti_Atom e1(g.x/l, g.y / l, g.z / l), e(1,0,0); // высчитываем e1
	if (is_equal(g.y, 0) && is_equal(g.z, 0)) e.set_V(0, 1, 0); // высчитываем e

	//cout << "e1 = " << e1 << "  L = " << e1.calculete_length() << "\n";
	//cout << "e = " <<e  << "  L = " << e.calculete_length() << "\n";

	Velociti_Atom e2 = e.vector_product(e1); // высчитываем e2, e3 спомощью векторного произведения
	//cout << "e2 = " << e2 <<"  L = " <<e2.calculete_length()<< "\n";
	e2.normalization(); // нормируе вектор так как если это не сделать то |g| != |g'| 
	//cout << "e2 = " << e2 << "  L = " << e2.calculete_length() << "\n";
	Velociti_Atom e3 = e1.vector_product(e2); //это вектор не нормирую он и так будет иметь длину 1
	//Velociti_Atom vector_product(Velociti_Atom& V ){ return Velociti_Atom( y*V.z - z*V.y , -x * V.z + z * V.x,  x * V.y - y * V.x ); } // векторное проиведение
	//cout << "e3 = " << e3 << "  L = " << e3.calculete_length() << "\n";	//if (fabs(e1.calculete_length() - 1) > 1e-5)
	//	cout << e1.calculete_length() <<"    fabs(e1.calculete_length() - 1)>1e-5\n";
	//if (fabs(e2.calculete_length() - 1) > 1e-5)
	//	cout << e2.calculete_length() << "    fabs(e2.calculete_length() - 1) > 1e-5\n";
	//if (fabs(e3.calculete_length() - 1) > 1e-5)
	//	cout << e3.calculete_length() << "    fabs(e3.calculete_length() - 1) > 1e-5\n";
	Velociti_Atom g_new( l*(e1.x*C1 + e2.x*C2 + e3.x*C3) , l*(e1.y*C1 + e2.y*C2 + e3.y*C3), l*(e1.z*C1 + e2.z*C2 + e3.z*C3) );

	//Velociti_Atom g_newN( e1.x*g_new.x +e2.x*g_new.y + e3.x*g_new.z, e1.y * g_new.x + e2.y*g_new.y + e3.y*g_new.z,  e1.z*g_new.x + e2.z*g_new.y + e3.z*g_new.z);
	
	return Velociti_Atom(g_new.x + (V.x + W.x)/2, g_new.y + (V.y + W.y) / 2, g_new.z + (V.z + W.z) / 2);
	//return Velociti_Atom(g_new.x + (V.x + W.x) / 2, g_new.y + (V.y + W.y) / 2, g_new.z + (V.z + W.z) / 2);
}
//cout << "e1 = " << e1 << " ; e2 = " << e2 << " ; e3 = " << e3 <<" ; g_new_ ="<< g_new_ << endl;
//cout << "V = " << V << " ; W =" << W <<" ; g = "<<g << " ; g_new = " << g_new << endl;