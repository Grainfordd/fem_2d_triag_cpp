#include "../include/elemento.hpp"

mat Elemento::obtener_factores(Nodo nodos[3]) {
    mat factores = zeros(3, 3);

    for (int i = 0; i < 3; i++) {
        double ai = nodos[(i+1)%3].x * nodos[(i+2)%3].y - nodos[(i+2)%3].x * nodos[(i+1)%3].y;
        double bi = nodos[(i+1)%3].y - nodos[(i+2)%3].y;
        double ci = nodos[(i+2)%3].x - nodos[(i+1)%3].x;

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
    return std::abs(area);
	/* return area; */
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

Elemento::Elemento(Nodo nodos_in[3], double E_in, double nu_in, double t_in) {
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
    K2 = K;
}
