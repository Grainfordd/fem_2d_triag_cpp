#include "../include/vtk_writer.hpp"
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

// Constructor
VTKWriter::VTKWriter() {
    // Constructor vacío
}

// Destructor
VTKWriter::~VTKWriter() {
    // Destructor vacío
}

// Función básica para exportar resultados esenciales
void VTKWriter::exportar_basico(const string& nombre_archivo, 
                               const vector<Nodo>& nodos, 
                               const vector<Elemento>& elementos, 
                               const vec& desplazamientos) {
    
    ofstream archivo(nombre_archivo);
    
    if (!archivo.is_open()) {
        cerr << "Error: No se pudo crear el archivo " << nombre_archivo << endl;
        return;
    }
    
    // Escribir header
    escribir_header(archivo, "Resultados básicos de Elementos Finitos");
    
    // Escribir geometría
    escribir_puntos(archivo, nodos);
    escribir_celdas(archivo, elementos);
    escribir_tipos_celda(archivo, elementos.size());
    
    // Escribir datos
    escribir_datos_nodos(archivo, nodos, desplazamientos);
    escribir_datos_elementos(archivo, elementos);
    
    archivo.close();
    cout << "Archivo VTK básico creado: " << nombre_archivo << endl;
}

// Función completa para exportar todos los resultados
void VTKWriter::exportar_completo(const string& nombre_archivo, 
                                 const vector<Nodo>& nodos, 
                                 const vector<Elemento>& elementos, 
                                 const vec& desplazamientos,
                                 const vec& fuerzas) {
    
    ofstream archivo(nombre_archivo);
    
    if (!archivo.is_open()) {
        cerr << "Error: No se pudo crear el archivo " << nombre_archivo << endl;
        return;
    }
    
    // Escribir header
    escribir_header(archivo, "Resultados completos de Elementos Finitos");
    
    // Escribir geometría
    escribir_puntos(archivo, nodos);
    escribir_celdas(archivo, elementos);
    escribir_tipos_celda(archivo, elementos.size());
    
    // Escribir datos completos
    escribir_datos_nodos_completo(archivo, nodos, desplazamientos, fuerzas);
    escribir_datos_elementos(archivo, elementos);
    
    archivo.close();
    cout << "Archivo VTK completo creado: " << nombre_archivo << endl;
}

// Función para exportar solo la malla
void VTKWriter::exportar_malla(const string& nombre_archivo,
                              const vector<Nodo>& nodos,
                              const vector<Elemento>& elementos) {
    
    ofstream archivo(nombre_archivo);
    
    if (!archivo.is_open()) {
        cerr << "Error: No se pudo crear el archivo " << nombre_archivo << endl;
        return;
    }
    
    // Escribir header
    escribir_header(archivo, "Malla de Elementos Finitos");
    
    // Escribir geometría
    escribir_puntos(archivo, nodos);
    escribir_celdas(archivo, elementos);
    escribir_tipos_celda(archivo, elementos.size());
    
    // Solo IDs de nodos y elementos
    archivo << "POINT_DATA " << nodos.size() << endl;
    archivo << "SCALARS Node_ID int 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (const auto& nodo : nodos) {
        archivo << nodo.id << endl;
    }
    archivo << endl;
    
    archivo << "CELL_DATA " << elementos.size() << endl;
    archivo << "SCALARS Element_ID int 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (const auto& elemento : elementos) {
        archivo << elemento.id << endl;
    }
    
    archivo.close();
    cout << "Archivo VTK de malla creado: " << nombre_archivo << endl;
}

// Funciones auxiliares privadas
void VTKWriter::escribir_header(ofstream& archivo, const string& descripcion) {
    archivo << "# vtk DataFile Version 3.0" << endl;
    archivo << descripcion << endl;
    archivo << "ASCII" << endl;
    archivo << "DATASET UNSTRUCTURED_GRID" << endl;
    archivo << endl;
}

void VTKWriter::escribir_puntos(ofstream& archivo, const vector<Nodo>& nodos) {
    archivo << "POINTS " << nodos.size() << " float" << endl;
    for (const auto& nodo : nodos) {
        archivo << fixed << setprecision(6) 
                << nodo.x << " " << nodo.y << " 0.0" << endl;
    }
    archivo << endl;
}

void VTKWriter::escribir_celdas(ofstream& archivo, const vector<Elemento>& elementos) {
    int num_elementos = elementos.size();
    archivo << "CELLS " << num_elementos << " " << num_elementos * 4 << endl;
    
    for (const auto& elemento : elementos) {
        archivo << "3 " 
                << elemento.nodos[0].id - 1 << " "  // VTK usa indexación base 0
                << elemento.nodos[1].id - 1 << " "
                << elemento.nodos[2].id - 1 << endl;
    }
    archivo << endl;
}

void VTKWriter::escribir_tipos_celda(ofstream& archivo, int num_elementos) {
    archivo << "CELL_TYPES " << num_elementos << endl;
    for (int i = 0; i < num_elementos; i++) {
        archivo << "5" << endl;  // 5 = triángulo en VTK
    }
    archivo << endl;
}

void VTKWriter::escribir_datos_nodos(ofstream& archivo, const vector<Nodo>& nodos, 
                                    const vec& desplazamientos) {
    int num_nodos = nodos.size();
    archivo << "POINT_DATA " << num_nodos << endl;
    
    // Desplazamientos como vectores
    archivo << "VECTORS Desplazamientos float" << endl;
    for (int i = 0; i < num_nodos; i++) {
        double dx = desplazamientos(i * 2);
        double dy = desplazamientos(i * 2 + 1);
        archivo << fixed << setprecision(6) 
                << dx << " " << dy << " 0.0" << endl;
    }
    archivo << endl;
    
    // Magnitud de desplazamientos
    archivo << "SCALARS Magnitud_Desplazamiento float 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < num_nodos; i++) {
        double dx = desplazamientos(i * 2);
        double dy = desplazamientos(i * 2 + 1);
        double magnitud = calcular_magnitud(dx, dy);
        archivo << fixed << setprecision(6) << magnitud << endl;
    }
    archivo << endl;
}

void VTKWriter::escribir_datos_nodos_completo(ofstream& archivo, const vector<Nodo>& nodos, 
                                             const vec& desplazamientos, const vec& fuerzas) {
    int num_nodos = nodos.size();
    archivo << "POINT_DATA " << num_nodos << endl;
    
    // Desplazamientos como vectores
    archivo << "VECTORS Desplazamientos float" << endl;
    for (int i = 0; i < num_nodos; i++) {
        double dx = desplazamientos(i * 2);
        double dy = desplazamientos(i * 2 + 1);
        archivo << fixed << setprecision(6) 
                << dx << " " << dy << " 0.0" << endl;
    }
    archivo << endl;
    
    // Fuerzas como vectores
    archivo << "VECTORS Fuerzas float" << endl;
    for (int i = 0; i < num_nodos; i++) {
        double fx = fuerzas(i * 2);
        double fy = fuerzas(i * 2 + 1);
        archivo << fixed << setprecision(6) 
                << fx << " " << fy << " 0.0" << endl;
    }
    archivo << endl;
    
    // Magnitud de desplazamientos
    archivo << "SCALARS Magnitud_Desplazamiento float 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < num_nodos; i++) {
        double dx = desplazamientos(i * 2);
        double dy = desplazamientos(i * 2 + 1);
        double magnitud = calcular_magnitud(dx, dy);
        archivo << fixed << setprecision(6) << magnitud << endl;
    }
    archivo << endl;
    
    // Magnitud de fuerzas
    archivo << "SCALARS Magnitud_Fuerza float 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (int i = 0; i < num_nodos; i++) {
        double fx = fuerzas(i * 2);
        double fy = fuerzas(i * 2 + 1);
        double magnitud = calcular_magnitud(fx, fy);
        archivo << fixed << setprecision(6) << magnitud << endl;
    }
    archivo << endl;
    
    // ID del nodo
    archivo << "SCALARS Node_ID int 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (const auto& nodo : nodos) {
        archivo << nodo.id << endl;
    }
    archivo << endl;
}

void VTKWriter::escribir_datos_elementos(ofstream& archivo, const vector<Elemento>& elementos) {
    int num_elementos = elementos.size();
    archivo << "CELL_DATA " << num_elementos << endl;
    
    // Esfuerzo de von Mises
    archivo << "SCALARS Von_Mises_Stress float 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (const auto& elemento : elementos) {
        archivo << fixed << setprecision(6) << elemento.esf_von_mises << endl;
    }
    archivo << endl;
    
    // Esfuerzos individuales
    archivo << "SCALARS Esfuerzo_X float 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (const auto& elemento : elementos) {
        archivo << fixed << setprecision(6) << elemento.esfuerzos(0) << endl;
    }
    archivo << endl;
    
    archivo << "SCALARS Esfuerzo_Y float 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (const auto& elemento : elementos) {
        archivo << fixed << setprecision(6) << elemento.esfuerzos(1) << endl;
    }
    archivo << endl;
    
    archivo << "SCALARS Esfuerzo_Cortante float 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (const auto& elemento : elementos) {
        archivo << fixed << setprecision(6) << elemento.esfuerzos(2) << endl;
    }
    archivo << endl;
    
    // Esfuerzos como vector 3D
    archivo << "VECTORS Esfuerzos float" << endl;
    for (const auto& elemento : elementos) {
        archivo << fixed << setprecision(6) 
                << elemento.esfuerzos(0) << " "  // sx
                << elemento.esfuerzos(1) << " "  // sy
                << elemento.esfuerzos(2) << endl; // tau
    }
    archivo << endl;
    
    // ID del elemento
    archivo << "SCALARS Element_ID int 1" << endl;
    archivo << "LOOKUP_TABLE default" << endl;
    for (const auto& elemento : elementos) {
        archivo << elemento.id << endl;
    }
    archivo << endl;
}

double VTKWriter::calcular_magnitud(double x, double y) {
    return sqrt(x*x + y*y);
}
