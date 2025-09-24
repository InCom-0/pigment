#include <spectral.hpp>

#include <iostream>

int main() {

    spectral_km::Color col1("#002185");
    spectral_km::Color col2("#FCD200");


    auto resCol = spectral_km::mix({{col1, 0.5}, {col2, 0.5}});

    std::cout << resCol.toString();
    return 0;
}