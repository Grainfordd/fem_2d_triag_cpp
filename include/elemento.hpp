#ifndef ELEMENTO_HPP
#define ELEMENTO_HPP

#include "nodo.hpp"
#include <armadillo>
#include <iostream>

using namespace arma;

class Elemento {
public:
    double E, nu, t;
    double area;
    mat D, K, K2, B, factores;
    Nodo nodos[3];

    Elemento(Nodo nodos_in[3] = 0, double E_in = 0, double nu_in = 0, double t_in = 0);
    
private:
    mat obtener_factores(Nodo nodos[3]);
    double calc_area(Nodo nodos[3]);
    mat mat_B(mat fact, double area);
};

#endif
