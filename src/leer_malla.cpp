#include "../include/leer_malla.hpp"
#include <fstream>
#include <algorithm>

std::vector<Nodo> leer_nodos(const std::string& filename) {
    std::vector<Nodo> nodos;
    std::ifstream file(filename);
    std::string line;
    
    while (getline(file, line)) {
        if (line == "$Nodes") {
            int num_nodos;
            file >> num_nodos;
            nodos.resize(num_nodos);
            
            for (int i = 0; i < num_nodos; ++i) {
                int id;
                double x, y, z;
                file >> id >> x >> y >> z;
                nodos[id-1] = Nodo(id, x, y);
            }
            break;
        }
    }
    return nodos;
}

std::vector<std::vector<int>> leer_elementos(const std::string& filename) {
    std::vector<std::vector<int>> elementos;
    std::ifstream file(filename);
    std::string line;
    
    while (getline(file, line)) {
        if (line == "$Elements") {
            int num_elementos;
            file >> num_elementos;
            
            for (int i = 0; i < num_elementos; ++i) {
                int id, type, num_tags;
                file >> id >> type >> num_tags;
                
                // Saltar los tags
                for (int j = 0; j < num_tags; ++j) {
                    int tag;
                    file >> tag;
                }
                
                if (type == 2) { // Solo triángulos (tipo 2)
                    std::vector<int> nodos_elemento(3);
                    file >> nodos_elemento[0] >> nodos_elemento[1] >> nodos_elemento[2];
                    elementos.push_back(nodos_elemento);
                } else {
                    // Saltar otros elementos
                    int dummy;
                    for (int j = 0; j < 3; ++j) file >> dummy;
                }
            }
            break;
        }
    }

    for (auto& vec : elementos){
        std::sort(vec.begin(), vec.end());
    }

    return elementos;
}
