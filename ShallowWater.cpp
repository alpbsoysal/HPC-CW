/**
 * @file ShallowWater.cpp
 * @author Alp Soysal
 * @brief A program that solves the shallow water equations using 6th order central differencing and 4th order runge-kutta, for the HPC module coursework.
 * @version 0.1
 * @date 2023-03-19
 */

#include <iostream>

using namespace std;

#include <boost/program_options.hpp>

namespace po = boost::program_options;

/**
 * @class ShallowWater
 * @brief Implements a numerical solution to the shallow water equations.
 */
class ShallowWater {

    private:
        int nx, ny;
        double dx, dy;
        double* u, v, h;

    public:

        void SetInitialConditions();

        void TimeIntegrate();

        void CalculateFluxLoop();

        void CalculateFluxBLAS();

        void deri_x(const double* var, double* der);

        void deri_y(const double* var, double* der);

};

void ShallowWater::SetInitialConditions() {
    return;
}

// Implement runge-kutta 4
void ShallowWater::TimeIntegrate() {
// calls CalculateFlux_loop() or CalculateFlux_matrix()
    return;
}

// Calculate the fluxes for all variables
void ShallowWater::CalculateFluxLoop() {
    // calls deri_x and deri_y
    return;
}

void ShallowWater::CalculateFluxBLAS() {
    return;
}

/**
 * @brief Function that calculates the x-derivative of a field VAR and stores it in DER using a loop-based method.
 * 
 * @param var Pointer to the input field matrix
 * @param der Pointer to the output matrix
 */
void ShallowWater::deri_x(const double* var, double* der) {

    for (int row = 0; row < ny; row++)
    {
        // first 3 elements in row
        der[0*ny + row] = (-1/60*var[(nx-3)*ny+row] + 3/20*var[(nx-2)*ny+row] - 3/4*var[(nx-1)*ny+row] + 3/4*var[ny+row] - 3/20*var[2*ny+row] + 1/60*var[3*ny+row])/dx;
        der[1*ny + row] = (-1/60*var[(nx-2)*ny+row] + 3/20*var[(nx-1)*ny+row] - 3/4*var[row] + 3/4*var[2*ny+row] - 3/20*var[3*ny+row] + 1/60*var[4*ny+row])/dx;
        der[2*ny + row] = (-1/60*var[(nx-1)*ny+row] + 3/20*var[row] - 3/4*var[ny+row] + 3/4*var[3*ny+row] - 3/20*var[4*ny+row] + 1/60*var[5*ny+row])/dx;

        // last 3 elements in row
        der[(nx-3)*ny + row] = (-1/60*var[(nx-6)*ny+row] + 3/20*var[(nx-5)*ny+row] - 3/4*var[(nx-4)*ny+row] + 3/4*var[(nx-2)*ny+row] - 3/20*var[(nx-1)*ny+row] + 1/60*var[row])/dx;
        der[(nx-2)*ny + row] = (-1/60*var[(nx-5)*ny+row] + 3/20*var[(nx-4)*ny+row] - 3/4*var[(nx-3)*ny+row] + 3/4*var[(nx-1)*ny+row] - 3/20*var[row] + 1/60*var[ny+row])/dx;
        der[(nx-1)*ny + row] = (-1/60*var[(nx-4)*ny+row] + 3/20*var[(nx-3)*ny+row] - 3/4*var[(nx-2)*ny+row] + 3/4*var[row] - 3/20*var[ny+row] + 1/60*var[2*ny+row])/dx;

        // elements between first 3 and last 3 elements
        for (int col = 3; col < nx-3; col++)
        {
            der[col*ny + row] = (-1/60*var[(col-3)*ny+row] + 3/20*var[(col-2)*ny+row] - 3/4*var[(col-1)*ny+row] + 3/4*var[(col+1)*ny+row] - 3/20*var[(col+2)*ny+row] + 1/60*var[(col+3)*ny+row])/dx;
        }
    }
}

/**
 * @brief Funciton that calculates the y-dericative of a field VAR and stores it in DER using a loop-based method.
 * 
 * @param var Pointer to the input field matrix
 * @param der Pointer to the output matrix
 */
void ShallowWater::deri_y(const double* var, double* der) {

    for (int col = 0; col < nx; col++)
    {
        // first 3 elements in row
        der[col*ny + 0] = (-1/60*var[col*ny+ny-3] + 3/20*var[col*ny+ny-2] - 3/4*var[col*ny+ny-1] + 3/4*var[col*ny+1] - 3/20*var[col*ny+2] + 1/60*var[col*ny+3])/dy;
        der[col*ny + 1] = (-1/60*var[col*ny+ny-2] + 3/20*var[col*ny+ny-1] - 3/4*var[col*ny] + 3/4*var[col*ny+2] - 3/20*var[col*ny+3] + 1/60*var[col*ny+4])/dy;
        der[col*ny + 2] = (-1/60*var[col*ny+ny-1] + 3/20*var[col*ny] - 3/4*var[col*ny+1] + 3/4*var[col*ny+3] - 3/20*var[col*ny+4] + 1/60*var[col*ny+5])/dy;

        // last 3 elements in row
        der[col*ny + ny-3] = (-1/60*var[col*ny+ny-6] + 3/20*var[col*ny+ny-5] - 3/4*var[col*ny+ny-4] + 3/4*var[col*ny+ny-2] - 3/20*var[col*ny+ny-1] + 1/60*var[col*ny])/dy;
        der[col*ny + ny-2] = (-1/60*var[col*ny+ny-5] + 3/20*var[col*ny+ny-4] - 3/4*var[col*ny+ny-3] + 3/4*var[col*ny+ny-1] - 3/20*var[col*ny] + 1/60*var[col*ny+1])/dy;
        der[col*ny + ny-1] = (-1/60*var[col*ny+ny-4] + 3/20*var[col*ny+ny-3] - 3/4*var[col*ny+ny-2] + 3/4*var[col*ny] - 3/20*var[col*ny+1] + 1/60*var[col*ny+2])/dy;

        // elements between first 3 and last 3 elements
        for (int row = 3; row < ny-3; row++)
        {
            der[col*ny + row] = (-1/60*var[col*ny+row-3] + 3/20*var[col*ny+row-2] - 3/4*var[col*ny+row-1] + 3/4*var[col*ny+row+1] - 3/20*var[col*ny+row+2] + 1/60*var[col*ny+row+3])/dy;
        }
    }
}

int main(int argc, char* argv[]) {

    po::options_description opts("Solves the shallow water equations");

    opts.add_options()
        ("help,h", "Print help message")
        ("dT", po::value<double>()->default_value(0.1), "Timestep to use")
        ("T", po::value<double>()->default_value(80), "Time to integrate for")
        ("Nx", po::value<int>()->default_value(100), "Number of points in X direction")
        ("Ny", po::value<int>()->default_value(100), "NUmber of points in Y direction")
        ("ic", po::value<int>(), "Initial condition to use, 1-4");

    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    if (vm.count("help"))
    {
        cout << opts << endl;
        return 0;
    }

    const double dT = vm["dT"].as<double>();
    const double T  = vm["T"].as<double>();
    const int Nx = vm["Nx"].as<int>();
    const int Ny = vm["Ny"].as<int>();

    if (vm.count("ic")){

        const int ic = vm["ic"].as<int>();

        if (ic < 1 || 4 < ic)
        {
            cout << "Initial conditions not valid. Valid range is 1-4. Exiting program." << endl;
            return -1;
        }

    } else {
        cout << "Initial conditions not specified. Exiting program." << endl;
        return -1;
    }
}