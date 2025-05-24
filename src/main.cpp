#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class Nodo {
	public:
		double x, y;
		int num_nodo;

		Nodo(int val_nodo=0, double x_val=0, double y_val=0){
			x = x_val;
			y = y_val;
			num_nodo = val_nodo;

		}
};

class Elemento { 
	public:
		double E, nu, t;
		double area;

		mat A_i;
		mat D = zeros(3,3);
		mat K ;
		mat B;
		mat factores;

		mat obtener_factores(Nodo nodos[3]){
			mat factores = zeros(3, 3);

			for (int i = 0; i < 3 ; i++){
				double ai = nodos[(i+1)%3].x * nodos[(i+2)%3].y - nodos[(i+2)%3].x * nodos[(i+1)%3].y;
				double bi = nodos[(i+1)%3].y - nodos[(i+2)%3].y;
				double ci = nodos[(i+2)%3].x - nodos[(i+1)%3].x;

				factores(0, i) = ai;
				factores(1, i) = bi;
				factores(2, i) = ci;
			}

			return factores;
		}

		double calc_area(Nodo nodos[3]){
			mat val = {
				{nodos[0].x, nodos[0].y, 1},
				{nodos[1].x, nodos[1].y, 1},
				{nodos[2].x, nodos[2].y, 1},
			};

			double area = 0.5 * det(val);
			return area;
		}

		mat mat_B(mat fact, double area){
			mat B =  {
				{fact(1,0), 0, fact(1,1), 0, fact(1,2), 0},
				{0, fact(2,0), 0, fact(2,1), 0, fact(2,2)},
				{fact(2,0), fact(1,0), fact(2,1), fact(1,1), fact(2,2), fact(1,2)}
			};
			B = B / (2*area);

			return B;
		}

		Elemento(Nodo nodos[3]= 0, double E_in=0, double nu_in=0, double t_in=0){
			E = E_in;
			nu = nu_in;
			t = t_in;
			area = calc_area(nodos);
			factores = obtener_factores(nodos);

			D = {
				{1, nu, 0},
				{nu, 1, 0},
				{0, 0, (1-nu)/2,}
			};

			D = D * E / (1 - nu*nu);
			B = mat_B(factores, area);
			K = t * area * B.t() * D * B;
		}
};


int main(void){
	double nodos_info[3][2] = {
		{0, 0},
		{4, 0},
		{0, 3}
	};

	Nodo nodos[3];
	double E = 200e3;
	double nu = 0.3;
	double t = 1;

	for (int i = 0; i < 3; i++){
		nodos[i] = Nodo(i+1, nodos_info[i][0], nodos_info[i][1]);
	}

	Elemento elemento1 = Elemento(nodos, E, nu, t);

	cout << elemento1.K << endl;
}
