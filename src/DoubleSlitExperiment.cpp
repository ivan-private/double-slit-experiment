#include "DoubleSlitExperiment.h"

#include <cmath>
#include <numbers>
#include <algorithm>
#include <sstream>
#include <iomanip>

using std::vector;

void
DoubleSlitExperiment::run()
{
    if (U.empty())
        U.resize(N * N);

    std::fill(U.begin(), U.end(), 0);

    vector<double> U_prev(N * N, 0.0);
    vector<double> laplacian(N * N, 0.0);
    vector<double> x_lin(N, 0.0);
    vector<bool> mask(N * N, false);

    const double dx = box_size / static_cast<double>(N);

    // making dt dependent on dx gives better accuracy
    const double dt = (std::sqrt(2)/2) * dx / wave_speed;

    // all constants the Laplacian multiplies with together in one factor
    const double fac = std::pow(dt, 2) * pow(wave_speed, 2) / pow(dx, 2);


    for (size_t i = 0; i < N; i++)
        x_lin[i] = 0.5*dx + static_cast<double>(i) * dx;

    initialize_mask(mask);

    double t = 0;
    while (t < T_end)
    {
        calculate_laplacian(laplacian, fac);

        calculate_next_U(laplacian, mask, U_prev);

        // apply boundary conditions (Dirichlet/inflow)
        // inflow from the top row
        for (size_t j = 1; j < N-1; j++)
            U[j] = std::sin(20 * std::numbers::pi * t) * std::pow(std::sin(std::numbers::pi * x_lin[j]) , 2);


        t += dt;
    }
}

void
DoubleSlitExperiment::initialize_mask(vector<bool> &mask) const
{
    for (size_t i = 0; i < N; i++)
    {
        mask[N*i] = true;
        mask[N*i + N-1] = true;
        mask[i] = true;
        mask[N*(N-1) + i] = true;
    }

    for (size_t i = N/4; i < (N*9) / 32; i++)
        for (size_t j = 0; j < N-1; j++)
            mask[i*N + j] = true;

    for (size_t i = 1; i < N-1; i++)
        for (size_t j = (N*5) / 16; j < (N*3) / 8; j++)
            mask[i*N + j] = false;

    for (size_t i = 1; i < N-1; i++)
        for (size_t j = (N*5) / 8; j < (N*11) / 16; j++)
            mask[i*N + j] = false;
}

void
DoubleSlitExperiment::calculate_laplacian(vector<double> &laplacian, const double fac) const
{
    // the edges are handled by wrapping around ()
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            laplacian[i*N + j] = fac * (U[i*N + (N+j-1) % N]
                                 + U[i*N + (j+1) % N]
                                 + U[((N+i-1) % N) * N + j]
                                 + U[((i+1) % N) * N + j]
                                 - 4 * U[i*N + j]);
        }
    }
}

void
DoubleSlitExperiment::calculate_next_U(const vector<double> &laplacian, const vector<bool> &mask,
                                       vector<double> &U_prev)
{
    double U_new;
    // take a time step and calculate new U based on the formula:
    // U_new = 2*U - U_prev + laplacian
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = 0; j < N; j++)
        {
            if (!mask[i * N + j])
            {
                U_new = 2 * U[i*N + j] - U_prev[i*N + j] + laplacian[i*N + j];
                U_prev[i*N + j] = U[i*N + j];
                U[i*N + j] = U_new;
            }
            else
                U[i*N + j] = 0.0;
        }
    }
}

std::string
DoubleSlitExperiment::results_as_string() const
{
    std::ostringstream oss;

    // Set formatting options for the output
    oss << std::fixed << std::setprecision(2); // Fixed-point notation with 2 decimal places

    for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < N; ++j) {
            const size_t index = i * N + j;
            if (index < U.size()) {
                oss << std::setw(10) << U[index];
            }
        }
        oss << "\n";
    }

    return oss.str();
}
