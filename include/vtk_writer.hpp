#ifndef VTK_WRITER_HPP
#define VTK_WRITER_HPP

#include <string>
#include <vector>
#include <armadillo>
#include "../include/nodo.hpp"
#include "../include/elemento.hpp"

using namespace std;
using namespace arma;

class VTKWriter {
public:
    // Constructor
    VTKWriter();
    
    // Destructor
    ~VTKWriter();
    
    // Función básica para exportar resultados esenciales
    void exportar_basico(const string& nombre_archivo, 
                        const vector<Nodo>& nodos, 
                        const vector<Elemento>& elementos, 
                        const vec& desplazamientos);
    
    // Función completa para exportar todos los resultados
    void exportar_completo(const string& nombre_archivo, 
                          const vector<Nodo>& nodos, 
                          const vector<Elemento>& elementos, 
                          const vec& desplazamientos,
                          const vec& fuerzas);
    
    // Función para exportar solo la malla (sin resultados)
    void exportar_malla(const string& nombre_archivo,
                       const vector<Nodo>& nodos,
                       const vector<Elemento>& elementos);

private:
    // Funciones auxiliares privadas
    void escribir_header(ofstream& archivo, const string& descripcion);
    void escribir_puntos(ofstream& archivo, const vector<Nodo>& nodos);
    void escribir_celdas(ofstream& archivo, const vector<Elemento>& elementos);
    void escribir_tipos_celda(ofstream& archivo, int num_elementos);
    void escribir_datos_nodos(ofstream& archivo, const vector<Nodo>& nodos, 
                             const vec& desplazamientos);
    void escribir_datos_nodos_completo(ofstream& archivo, const vector<Nodo>& nodos, 
                                      const vec& desplazamientos, const vec& fuerzas);
    void escribir_datos_elementos(ofstream& archivo, const vector<Elemento>& elementos);
    double calcular_magnitud(double x, double y);
};

#endif // VTK_WRITER_HPP
