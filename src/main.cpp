#include <armadillo>
/* #include <iomanip> */
#include <string>
#include <vector>

#include "../include/elemento.hpp"
#include "../include/nodo.hpp"
#include "../include/leer_malla.hpp"
#include "../include/utils.hpp"
#include "../include/vtk_writer.hpp"

using namespace std;
using namespace arma;

void mostrar_elementos(vector<Elemento> elementos, int num_elementos);
void mostrar_desp(vec desp, int num_desp);
void ordenar_nodos(Nodo nodos[], int tamaño);
vec calc_esf(vec desp, Elemento elemento);
double von_mises(vec esfuerzos);

int main(void){
	double E = 200e3;
	double nu = 0.3;
	double t = 1;

	// ---------------------------  Leer nodos y elementos ----------------------------------
	string malla = "../data/malla.msh";
	vector<Nodo> nodos = leer_nodos(malla);
	vector<vector<int>> elementos_indices = leer_elementos(malla);


	vector<Elemento> elementos;
	int i = 1;
	for (const auto& indices : elementos_indices) {
		Nodo nodos_elemento[] = {nodos[indices[0]-1], nodos[indices[1]-1], nodos[indices[2]-1]};

		cout << i << " " << " " << E<< " " << nu << " " << t << endl;

		elementos.emplace_back(i, nodos_elemento, E, nu, t);
		i++;
	}
	// -------------------------------------------------------------------------------------

	int n = nodos.size() * 2;
	int num_elementos = elementos.size();


	// Ensamblar matriz de rigidez global
	mat k_glob = zeros(n, n);

	int n1;
	int n2;
	int n3;

	for (int i = 0; i < num_elementos; i++){
		n1 = elementos[i].nodos[0].id;
		n2 = elementos[i].nodos[1].id;
		n3 = elementos[i].nodos[2].id;

		vector<int>  indices_glob= {n1*2-2, n1*2-1, n2*2-2, n2*2-1, n3*2-2, n3*2-1};
		for (int j = 0; j < 6 ; j++) {
			for (int k = 0; k < 6; k++){
				k_glob(indices_glob[j], indices_glob[k]) += elementos[i].K(j,k);
			}
		}

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


	vec disp_red = solve(k_red, f_red);
	int k = disp_red.n_rows-1;

	for (int i = disp.n_rows-1; i >= 0; i--){
		if (disp[i] != 0){
			disp[i] = disp_red[k];
			k--;
		}
	}

	vec F = k_glob * disp;
	/* F.print(); */

	/* disp.print(); */
	for (int i = 0 ; i <  num_elementos ; i++){
		ordenar_nodos(elementos[i].nodos, 3);
		vec disp_elem = zeros(6);
		for (int j = 0; j < 3; j++) {
			disp_elem(j * 2) = disp(elementos[i].nodos[j].id * 2 - 2);      // componente x
			disp_elem(j * 2 + 1) = disp(elementos[i].nodos[j].id * 2 - 1);  // componente y
    	}

		elementos[i].esfuerzos = calc_esf(disp_elem, elementos[i]);
		elementos[i].esf_von_mises = von_mises(elementos[i].esfuerzos);

		/* cout << elementos[i].esf_von_mises << endl; */
	}

	/* mostrar_desp(desplazamientos, desplazamientos.n_rows); */
	/* mostrar_elementos(elementos, num_elementos); */
	mostrar_desp(disp,  n);
	/* mostrar_desp(F,  n); */
	VTKWriter vtk_writer;
    
    // Crear directorio de resultados si no existe
    system("mkdir -p ../results");
    
    // Opción 1: Exportar resultados básicos
    string archivo = "../results/resultados.vtk";
    vtk_writer.exportar_basico(archivo, nodos, elementos, disp);
	/* k_glob.print(); */
	/* cout << elementos[0].K.is_symmetric() << endl; */
/* cout << k_glob.is_symmetric() << endl; */

	/* F.print(); */

	return 0;
}


vec calc_esf(vec desp, Elemento elemento){
	vec esfuerzos = elemento.D * elemento.B * desp;
	return esfuerzos;
}

double von_mises(vec esfuerzos){
	double sx = esfuerzos[0];
	double sy = esfuerzos[1];
	double tau = esfuerzos[2];
	double esf_von_mises = sqrt(sx*sx + sy*sy - sx*sy + 3*tau*tau);
	return esf_von_mises;

}

void mostrar_desp(vec desp, int num_desp){
	int j = 1;
	for (int i = 0; i < num_desp-4 + 4; i+=2){
		cout << j <<  ": "<< "( " << desp[i] << " , " << desp[i+1] << " )" <<  endl;
		j++;
	}
}


void mostrar_elementos(vector<Elemento> elementos, int num_elementos){
	for (int i = 0; i < num_elementos ; i++){
		cout << endl << "###################### Elemento " << elementos[i].id<<" ##############################################" << endl;
		cout << "Nodos" << endl ;

		for (int j = 0; j < 3; j++){
			Nodo nodo_ = elementos[i].nodos[j];
			
			cout << nodo_.id << "   (" << nodo_.x << " , " << nodo_.y << ")"  << endl;
		}
		cout << endl;

		elementos[i].K.print();
		cout << endl << "################################################################################" << endl;
	}

}

void ordenar_nodos(Nodo nodos[], int tamaño) {
    for (int i = 0; i < tamaño - 1; i++) {
        for (int j = 0; j < tamaño - 1 - i; j++) {
            if (nodos[j].id > nodos[j + 1].id) {
                // Intercambiar nodos
                Nodo temp = nodos[j];
                nodos[j] = nodos[j + 1];
                nodos[j + 1] = temp;
            }
        }
    }
}
