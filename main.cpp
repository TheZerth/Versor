#include <iostream>
#include "libs/Versor.h"

int main() {

    VRSR::Versor v1 = {1.0, 1.0, 1.0, 1.0};
    VRSR::Versor v2 = {2.0, 2.0, 2.0, 2.0};
    VRSR::Versor v3 = {3.0, 3.0, 3.0, 3.0};
    VRSR::Versor v4 = {4.0, 4.0, 4.0, 4.0};
    VRSR::Versor v5 = {5.0, 5.0, 5.0, 5.0};
    VRSR::Versor v6 = {6.0, 6.0, 6.0, 6.0};
    VRSR::Versor v7 = {7.0, 7.0, 7.0, 7.0};
    VRSR::Versor v8 = {8.0, 8.0, 8.0, 8.0};
    VRSR::Versor v9 = {9.0, 9.0, 9.0, 9.0};
    VRSR::Versor v10 = {10.0, 10.0, 10.0, 10.0};
    VRSR::Versor versors[10] = {v1, v2, v3, v4, v5, v6, v7, v8, v9, v10};

    for (auto &versor : versors) {
        std::cout << (versor+versor) << std::endl;
        std::cout << (versor+3.14159) << std::endl;
        std::cout << (versor+3.14159f) << std::endl;
        std::cout << (versor+3l) << std::endl;
    }

    return 0;
}