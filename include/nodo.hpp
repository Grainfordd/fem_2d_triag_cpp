#ifndef NODO_HPP
#define NODO_HPP

class Nodo {
public:
    double x, y;
    int id;

    Nodo(int val_nodo = 0, double x_= 0, double y_= 0);

	~Nodo(){};
};


#endif
