#include <armadillo>
#include <vector>

#include "../include/elemento.hpp"
#include "../include/nodo.hpp"
#include "../include/leer_malla.hpp"
#include "../include/utils.hpp"

using namespace std;
using namespace arma;

/* mat eliminar_fila(mat A, int fila) { */
/*     int largo = A.n_rows; */
/*     if (fila == 0) return A.rows(fila + 1, largo - 1); */
/*     if (fila + 1 == largo) return A.rows(0, fila - 1); */
/*     return join_vert(A.rows(0, fila - 1), A.rows(fila + 1, largo - 1)); */
/* } */


/* mat eliminar_columna(mat A, int columna) { */
/*     int largo = A.n_cols; */
/*     if (columna == 0) return A.cols(columna + 1, largo - 1); */
/*     if (columna + 1 == largo) return A.cols(0, columna - 1); */
/*     return join_horiz(A.cols(0, columna - 1), A.cols(columna + 1, largo - 1)); */
/* } */


/* vector<Nodo> leer_nodos(const string& filename) { */
/*     vector<Nodo> nodos; */
/*     ifstream file(filename); */
/*     string line; */
    
/*     while (getline(file, line)) { */
/*         if (line == "$Nodes") { */
/*             int num_nodos; */
/*             file >> num_nodos; */
/*             nodos.resize(num_nodos); */
            
/*             for (int i = 0; i < num_nodos; ++i) { */
/*                 int id; */
/*                 double x, y, z; */
/*                 file >> id >> x >> y >> z; */
/*                 nodos[id-1] = Nodo(id, x, y); */
/*             } */
/*             break; */
/*         } */
/*     } */
/*     return nodos; */
/* } */

/* vector<vector<int>> leer_elementos(const string& filename) { */
/*     vector<vector<int>> elementos; */
/*     ifstream file(filename); */
/*     string line; */
    
/*     while (getline(file, line)) { */
/*         if (line == "$Elements") { */
/*             int num_elementos; */
/*             file >> num_elementos; */
            
/*             for (int i = 0; i < num_elementos; ++i) { */
/*                 int id, type, num_tags; */
/*                 file >> id >> type >> num_tags; */
                
/*                 // Saltar los tags */
/*                 for (int j = 0; j < num_tags; ++j) { */
/*                     int tag; */
/*                     file >> tag; */
/*                 } */
                
/*                 if (type == 2) { // Solo nos interesan los triángulos (tipo 2) */
/*                     vector<int> nodos_elemento(3); */
/*                     file >> nodos_elemento[0] >> nodos_elemento[1] >> nodos_elemento[2]; */
/*                     elementos.push_back(nodos_elemento); */
/*                 } else { */
/*                     // Saltar otros tipos de elementos */
/*                     int dummy; */
/*                     for (int j = 0; j < 3; ++j) file >> dummy; */
/*                 } */
/*             } */
/*             break; */
        /* } */
    /* } */

	/* for (auto& vec : elementos){ */
		/* sort(vec.begin(), vec.end()); */
	/* } */

    /* return elementos; */
/* } */
/* double calc_area(Nodo nodos_elem[3]) { */
/*     mat val = { */
/*         {nodos_elem[0].x, nodos_elem[0].y, 1}, */
/*         {nodos_elem[1].x, nodos_elem[1].y, 1}, */
/*         {nodos_elem[2].x, nodos_elem[2].y, 1} */
/*     }; */

/*     double area = 0.5 * det(val); */
/* 	return area; */
/* } */

/* vec eliminar_elem_vec(vec fuerzas, int indice){ */
/* 	vec nuevo_f = join_vert(fuerzas.head(indice), fuerzas.tail(fuerzas.n_elem - indice - 1)); */
/*     return nuevo_f; */
	
/* } */

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
