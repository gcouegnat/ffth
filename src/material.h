#ifndef MATERIAL_H_Y80J6HTZ
#define MATERIAL_H_Y80J6HTZ

#include <iostream>

struct Material {
    double props[42];

    void print() const {
        for(int i = 0; i < 42; ++i) {
            std::cout << props[i] << " ";
        }
        std::cout << std::endl;
    }
};

struct Phase {
    int matid;
    double angles[3];
};


#endif /* end of include guard: MATERIAL_H_Y80J6HTZ */

