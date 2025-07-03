#include <armadillo>
/* #include <iomanip> */
#include <string>
#include <vector>

#include "../include/elemento.hpp"
#include "../include/nodo.hpp"
#include "../include/leer_malla.hpp"
#include "../include/utils.hpp"

using namespace std;
using namespace arma;

/* void mostrar_elementos(vector<Elemento> elementos, int num_elementos); */
void mostrar_desp(vec desp, int num_desp);
vec calc_esf(vec desp, Elemento elemento);
double von_mises(vec esfuerzos);

void escribir_vtk(const std::string& filename, 
                  const std::vector<Nodo>& nodos, 
                  const std::vector<Elemento>& elementos,
                  const arma::vec& desplazamientos) {
    
    std::ofstream vtk_file(filename);

    if (!vtk_file.is_open()) {
        std::cerr << "Error: No se pudo abrir el archivo para escribir: " << filename << std::endl;
        return;
    }

    // 1. Cabecera del archivo VTK
    vtk_file << "# vtk DataFile Version 3.0" << std::endl;
    vtk_file << "Resultados de Analisis de Elementos Finitos" << std::endl;
    vtk_file << "ASCII" << std::endl;
    vtk_file << "DATASET UNSTRUCTURED_GRID" << std::endl;
    vtk_file << std::endl;

    // 2. Escribir las coordenadas de los nodos (Puntos)
    vtk_file << "POINTS " << nodos.size() << " double" << std::endl;
    for (const auto& nodo : nodos) {
        vtk_file << nodo.x << " " << nodo.y << " " << 0.0 << std::endl;
    }
    vtk_file << std::endl;

    // 3. Escribir la conectividad de los elementos (Celdas)
    vtk_file << "CELLS " << elementos.size() << " " << elementos.size() * 4 << std::endl;
    for (const auto& elem : elementos) {
        vtk_file << "3 " << elem.nodos_id[0] - 1 
                 << " " << elem.nodos_id[1] - 1 
                 << " " << elem.nodos_id[2] - 1 << std::endl;
    }
    vtk_file << std::endl;

    // 4. Escribir los tipos de elementos (Todos son triángulos)
    vtk_file << "CELL_TYPES " << elementos.size() << std::endl;
    for (size_t i = 0; i < elementos.size(); ++i) {
        vtk_file << "5" << std::endl; // 5 es el código para VTK_TRIANGLE
    }
    vtk_file << std::endl;

    // 5. Escribir datos asociados a los puntos (Desplazamientos)
    vtk_file << "POINT_DATA " << nodos.size() << std::endl;
    vtk_file << "VECTORS Desplazamientos double" << std::endl;
    for (const auto& nodo : nodos) {
        int idx_u = nodo.id * 2 - 2;
        int idx_v = nodo.id * 2 - 1;
        vtk_file << desplazamientos(idx_u) << " " << desplazamientos(idx_v) << " 0.0" << std::endl;
    }
    vtk_file << std::endl;
    
    // 6. Escribir datos asociados a las celdas (Esfuerzos)
    vtk_file << "CELL_DATA " << elementos.size() << std::endl;
    // Esfuerzo de Von Mises (escalar)
    vtk_file << "SCALARS Esfuerzo_Von_Mises double 1" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;
    for (const auto& elem : elementos) {
        vtk_file << elem.esf_von_mises << std::endl;
    }
    vtk_file << std::endl;

    // Esfuerzos Sx, Sy, Txy (como campos escalares separados)
    vtk_file << "SCALARS Sigma_x double 1" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;
    for (const auto& elem : elementos) {
        vtk_file << elem.esfuerzos(0) << std::endl;
    }
    vtk_file << std::endl;

    vtk_file << "SCALARS Sigma_y double 1" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;
    for (const auto& elem : elementos) {
        vtk_file << elem.esfuerzos(1) << std::endl;
    }
    vtk_file << std::endl;

    vtk_file << "SCALARS Tau_xy double 1" << std::endl;
    vtk_file << "LOOKUP_TABLE default" << std::endl;
    for (const auto& elem : elementos) {
        vtk_file << elem.esfuerzos(2) << std::endl;
    }
    vtk_file << std::endl;

    vtk_file.close();
    std::cout << "-> Archivo VTK '" << filename << "' generado exitosamente." << std::endl;
}


int main(void){
	double E = 200e3;
	double nu = 0.3;
	double t = 1;

	// ---------------------------  Leer nodos y elementos ----------------------------------
	string malla = "../data/malla.msh";
	vector<Nodo> nodos = leer_nodos(malla);
	vector<vector<int>> elementos_indices = leer_elementos(malla);

	cout << elementos_indices[0][0] << endl;
	cout << elementos_indices[0][1] << endl;
	cout << elementos_indices[0][2] << endl;
	return 0;


	vector<Elemento> elementos;
	int i = 1;
	for (const auto& indices : elementos_indices) {

		elementos.emplace_back(i, indices, nodos, E, nu, t);
		i++;
	}
	// -------------------------------------------------------------------------------------
	int n = nodos.size() * 2;
	int num_elementos = elementos.size();

	// Ensamblar matriz de rigidez global
	mat k_glob = zeros(n, n);

	int n1, n2, n3;
	for (int i = 0; i < num_elementos; i++){
		n1 = nodos[elementos[i].nodos_id[0]-1].id;
		n2 = nodos[elementos[i].nodos_id[1]-1].id;
		n3 = nodos[elementos[i].nodos_id[2]-1].id;

		vector<int>  indices_glob= {
			n1*2-2, n1*2-1,
			n2*2-2, n2*2-1,
			n3*2-2, n3*2-1
		};
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

	for (int i = 0 ; i <  num_elementos ; i++){
		vec disp_elem = zeros(6);
		for (int j = 0; j < 3; j++) {
			disp_elem(j*2) = disp(nodos[elementos[i].nodos_id[j] - 1].id * 2 -2);
			disp_elem(j*2+1) = disp(nodos[elementos[i].nodos_id[j] - 1].id * 2 -1);
    	}
		elementos[i].esfuerzos = calc_esf(disp_elem, elementos[i]);
		elementos[i].esf_von_mises = von_mises(elementos[i].esfuerzos);
	}

	mostrar_desp(disp,  n);
	escribir_vtk("resultados.vtk", nodos, elementos, disp);

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

