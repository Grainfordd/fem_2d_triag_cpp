#ifndef FEM_UTILS_H
#define FEM_UTILS_H

#include <armadillo>
#include "nodo.hpp"

arma::mat eliminar_fila(const arma::mat& A, int fila);
arma::mat eliminar_columna(const arma::mat& A, int columna);
arma::vec eliminar_elem_vec(const arma::vec& fuerzas, int indice);
double calc_area(const Nodo nodos_elem[3]);

#endif 
