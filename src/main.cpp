#include <iostream>

#include "DoubleSlitExperiment.h"

int main() {
    size_t N = 256;
    double box_size = 1;
    double wave_speed = 1;
    double T_end = 2;

    DoubleSlitExperiment experiment(N, box_size, T_end, wave_speed);
    experiment.run();

    std::cout << "Results:\n" << experiment.results_as_string();

}
