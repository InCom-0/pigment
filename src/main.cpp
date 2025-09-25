#include <incom/pigment/pigment.hpp>

#include <iostream>

int main() {

    incpigm::Color col1("#ff0000");
    incpigm::Color col2("#00FF00");


    std::vector<std::pair<std::tuple<int, int, int>, double>> test{{{255, 0, 0}, 2}, {{0, 255, 0}, 2}};


    auto resCol = incpigm::Color::mix({{col1, 1}, {col2, 1}});
    auto resCol_2 = incpigm::mix_sRBG_f(test);


    std::cout << resCol.toString();
    return 0;
}