#include "../include/elemento.hpp"
#include <armadillo>

mat Elemento::obtener_factores(Nodo nodos[3]) {
    mat factores = zeros(3, 3);

    for (int i = 0; i < 3; i++) {
		int j = (i+1)%3;
		int k = (i+2)%3;

        double ai = nodos[j].x * nodos[k].y - nodos[k].x * nodos[j].y;
        double bi = nodos[j].y - nodos[k].y;
        double ci = nodos[k].x - nodos[j].x;

        factores(0, i) = ai;
        factores(1, i) = bi;
        factores(2, i) = ci;
    }

    return factores;
}

double Elemento::calc_area(Nodo nodos[3]) {
    mat val = {
        {nodos[0].x, nodos[0].y, 1},
        {nodos[1].x, nodos[1].y, 1},
        {nodos[2].x, nodos[2].y, 1},
    };

    double area = 0.5 * det(val);
    return area;
}

mat Elemento::mat_B(mat fact, double area) {
    mat B = {
        {fact(1,0), 0, fact(1,1), 0, fact(1,2), 0},
        {0, fact(2,0), 0, fact(2,1), 0, fact(2,2)},
        {fact(2,0), fact(1,0), fact(2,1), fact(1,1), fact(2,2), fact(1,2)}
    };
    B = B / (2*area);
    return B;
}

Elemento::Elemento(int id_,Nodo nodos_in[3], double E_in, double nu_in, double t_in) {
	id = id_;
    E = E_in;
    nu = nu_in;
    t = t_in;
    area = calc_area(nodos_in);

    factores = obtener_factores(nodos_in);

    for (int i = 0; i < 3; i++) {
        nodos[i] = nodos_in[i];
    }

    D = {
        {1, nu, 0},
        {nu, 1, 0},
        {0, 0, (1-nu)/2}
    };

    D = D * E / (1 - nu*nu);
    B = mat_B(factores, area);
    K = t * area * B.t() * D * B;
}
