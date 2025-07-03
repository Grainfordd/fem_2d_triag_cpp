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
	mat factores = obtener_factores(nodos);
    mat val = {
        {nodos[0].x, nodos[0].y, 1},
        {nodos[1].x, nodos[1].y, 1},
        {nodos[2].x, nodos[2].y, 1},
    };


    double area = 0.5 * det(val);
    return area;
}

mat Elemento::mat_B(mat fact, double area) {

	mat B = zeros(3, 6);
	for (int i = 0; i < 3; i++){

		mat bi = {
			{fact(1, i), 0},
			{0, fact(2, i)},
			{fact(2, i), fact(1,i)}
		};

		bi = 1 /(2*area) * bi;

		for (int j = 0; j < 3; j++){
			B(j,i*2 ) = bi(j, 0);
			B(j,i*2 +1) = bi(j, 1);
		}
	}

    return B;
}


Elemento::Elemento(int id_, std::vector<int> nodos_id_in,  std::vector<Nodo> global_nodos, double E_in, double nu_in, double t_in) {
	/* cout << "hola" << endl; */
	id = id_;
    E = E_in;
    nu = nu_in;
    t = t_in;

	for (int i = 0; i < 3; i++){
		nodos_id.push_back(nodos_id_in[i]) ;
	}

	Nodo nodos[3];
	nodos[0] = global_nodos[nodos_id_in[0]-1];
	nodos[1] = global_nodos[nodos_id_in[1]-1];
	nodos[2] = global_nodos[nodos_id_in[2]-1];

    area = calc_area(nodos);

    factores = obtener_factores(nodos);

    D = {
        {1, nu, 0},
        {nu, 1, 0},
        {0, 0, (1-nu)/2}
    };

    D = D * E / (1 - nu*nu);
    B = mat_B(factores, area);
    K = t * area * B.t() * D * B;
}
