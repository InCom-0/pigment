#include <incom/pigment/pigment.hpp>

#include <iostream>

int main() {

    incpig::Color col1("#002185");
    incpig::Color col2("#FCD200");


    auto resCol = incpig::Color::mix({{col1, 0.5}, {col2, 0.5}});

    std::cout << resCol.toString();
    return 0;
}