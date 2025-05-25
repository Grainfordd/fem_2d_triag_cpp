#ifndef MESH_READER_H
#define MESH_READER_H

#include <vector>
#include <string>
#include "nodo.hpp"

std::vector<Nodo> leer_nodos(const std::string& filename);
std::vector<std::vector<int>> leer_elementos(const std::string& filename);

#endif // MESH_READER_H
