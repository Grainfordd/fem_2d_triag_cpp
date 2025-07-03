#ifndef ELEMENTO_HPP
#define ELEMENTO_HPP

#include "nodo.hpp"
#include <armadillo>
#include <cmath>
#include <algorithm>


using namespace arma;

class Elemento {
public:
    double E, nu, t;
    double area;
    mat D, K, B, factores;
    /* Nodo nodos[3]; */
	int id;
	vec esfuerzos;
	double esf_von_mises;
	std::vector<int> nodos_id;

    Elemento(int id_, std::vector<int>nodos_id_in, std::vector<Nodo> global_nodos, double E_in = 0, double nu_in = 0, double t_in = 0);

	~Elemento(){};
    
private:
	Nodo ordernar_nodos(Nodo nodos[3]);
    mat obtener_factores(Nodo nodos[3]);
    double calc_area(Nodo nodos[3]);
    mat mat_B(mat fact, double area);

};

#endif
