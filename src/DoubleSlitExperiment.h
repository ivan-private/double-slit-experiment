#pragma once

#include <vector>
#include <string>

class DoubleSlitExperiment {

// Public methods
public:
    DoubleSlitExperiment(size_t domain_num_vertices,
                         double domain_size,
                         double T_end,
                         double wave_speed
                         ) :

                         N(domain_num_vertices),
                         box_size(domain_size),
                         T_end(T_end),
                         wave_speed(wave_speed)
    {}

    void run();

    std::string results_as_string() const;


// Private methods
private:
    void initialize_mask(std::vector<bool> &mask) const;

    void calculate_laplacian(std::vector<double> &laplacian, double fac) const;

    void calculate_next_U(const std::vector<double> &laplacian, const std::vector<bool> &mask,
            std::vector<double> &U_prev);

// Private members
private:
    std::vector<double> U; // lazy initialization, gets initialized when run() is called
    const size_t N;
    const double box_size, T_end, wave_speed;

};

