#include <armadillo>
#include <vector>

#include "../include/elemento.hpp"
#include "../include/nodo.hpp"
#include "../include/leer_malla.hpp"
#include "../include/utils.hpp"

using namespace std;
using namespace arma;

int main(void){
	double E = 200e3;
	double nu = 0.3;
	double t = 1;

	// ---------------------------  Leer nodos y elementos ----------------------------------
	string malla = "../data/malla.msh";
	vector<Nodo> nodos = leer_nodos(malla);
	vector<vector<int>> elementos_indices = leer_elementos(malla);

	vector<Elemento> elementos;
	for (const auto& indices : elementos_indices) {
		Nodo nodos_elemento[3] = {nodos[indices[0]-1], nodos[indices[1]-1], nodos[indices[2]-1]};
		if (calc_area(nodos_elemento) < 0){
			swap(nodos_elemento[1], nodos_elemento[2]);
		}
		elementos.emplace_back(nodos_elemento, E, nu, t);
	}
	// -------------------------------------------------------------------------------------

	int n = nodos.size() * 2;
	int num_elementos = elementos.size();

	//  ---------------------- Expandir matrices de rigidez para cada elemento ------------
	for (int k = 0; k < num_elementos; k++){
		for (int i = 1; i <= n ; i++){
			bool cond1 = i != elementos[k].nodos[0].id*2 -1 && i != elementos[k].nodos[0].id*2;
			bool cond2 = i != elementos[k].nodos[1].id*2 -1 && i != elementos[k].nodos[1].id*2;
			bool cond3 = i != elementos[k].nodos[2].id*2 -1 && i != elementos[k].nodos[2].id*2;

			if (cond1 && cond2 && cond3){
				elementos[k].K2.insert_rows(i-1, 1);
				elementos[k].K2.insert_cols(i-1, 1);
			}
		}
	}
	// --------------------------------------------------------------------------------------

	// Ensamblar matriz de rigidez global
	mat k_glob = zeros(n, n);
	for (int i = 0; i < num_elementos; i++){
		k_glob += elementos[i].K2;
	}
    // ----------------------------- Leer Fuerzas -------------------------------------------- 
	mat f_read;
	f_read.load(csv_name("../data/fuerzas.csv", {"Nodo", "dir", "valor"}));
	vec f = zeros(n);

	for (int i = 0; i < f_read.n_rows; i++){
		int nodo = f_read(i, 0);
		int dir = f_read(i ,1);
		double indice = (dir == 0) ? nodo*2 - 2 : nodo*2 -1;
		f(indice) = f_read(i, 2);
	}
	vec f_red = f;
	// ---------------------------------------------------------------------------------------

	// ---------------------- Leer desplazamientos y reducir vectore de fuerzas --------------

	vec disp = ones(n);
	mat disp_read;
	disp_read.load(csv_name("../data/disp.csv", {"Nodo", "dir", "valor"}));

	for (int i = disp_read.n_rows-1; i >= 0; i--){
		int nodo = disp_read(i, 0);
		int dir = disp_read(i ,1);
		double indice = (dir == 0) ? nodo*2 - 2 : nodo*2 -1;
		f_red = eliminar_elem_vec(f_red, indice);
		disp(indice) = 0;
	}
	// ---------------------------------------------------------------------------------------

	// Reducir matriz de rigidez
	mat k_red = k_glob;
	for (int i = n-1 ; i >= 0; i--){
		if (disp[i] == 0){
			k_red = eliminar_fila(k_red, i);
			k_red = eliminar_columna(k_red, i);
		}
	}
	vec desplazamientos = solve(k_red, f_red);
	desplazamientos.print();
}
