#include "../include/utils.hpp"

arma::mat eliminar_fila(const arma::mat& A, int fila) {
    int largo = A.n_rows;
    if (fila == 0) return A.rows(fila + 1, largo - 1);
    if (fila + 1 == largo) return A.rows(0, fila - 1);
    return arma::join_vert(A.rows(0, fila - 1), A.rows(fila + 1, largo - 1));
}

arma::mat eliminar_columna(const arma::mat& A, int columna) {
    int largo = A.n_cols;
    if (columna == 0) return A.cols(columna + 1, largo - 1);
    if (columna + 1 == largo) return A.cols(0, columna - 1);
    return arma::join_horiz(A.cols(0, columna - 1), A.cols(columna + 1, largo - 1));
}

arma::vec eliminar_elem_vec(const arma::vec& fuerzas, int indice) {
    return arma::join_vert(fuerzas.head(indice), fuerzas.tail(fuerzas.n_elem - indice - 1));
}

double calc_area(const Nodo nodos_elem[3]) {
    arma::mat val = {
        {nodos_elem[0].x, nodos_elem[0].y, 1},
        {nodos_elem[1].x, nodos_elem[1].y, 1},
        {nodos_elem[2].x, nodos_elem[2].y, 1}
    };
    return 0.5 * arma::det(val);
}
