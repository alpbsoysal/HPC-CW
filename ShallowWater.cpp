/**
 * @file ShallowWater.cpp
 * @author Alp Soysal
 * @brief A program that solves the shallow water equations using 6th order central differencing and 4th order runge-kutta, for the HPC module coursework.
 * @version 0.1
 * @date 2023-03-19
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

using namespace std;

#include <boost/program_options.hpp>
#include <cblas.h>

namespace po = boost::program_options;

constexpr double G = 9.81;

void PrintMatrix(int n, int m, double* H) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < m; ++j) {
            cout << setw(10) << H[j*n+i] << " ";
        }
        cout << endl;
    }
}

void PrintVector(int n, double* u) {
    for (int i = 0; i < n; ++i) {
        cout << u[i] << endl;
    }
}

/**
 * @class ShallowWater
 * @brief Implements a numerical solution to the shallow water equations.
 */
class ShallowWater {

    private:
        int nx, ny;
        int method;
        double dx, dy, dt;
        double x0, x1, y0, y1, T;
        double *u, *v, *h, *A;

    public:

        ShallowWater(const double& pDT, const double& pT, const int& pNx, const int& pNy, const int& pIc, const int& pMethod);
        ~ShallowWater();

        void SetInitialConditions(const int& ic);
        void TimeIntegrate();

        void CalculateFluxLoop(double* pU, double* pV, double* pH, double* pKU, double* pKV, double* pKH);
        void DeriXLoop(const double* var, double* der);
        void DeriYLoop(const double* var, double* der);

        void CalculateFluxBLAS(double* pU, double* pV, double* pH, double* pKU, double* pKV, double* pKH);
        void PopulateMatrix();
        void DeriXBLAS(const double* var, double* der);
        void DeriYBLAS(const double* var, double* der);

        void FileOutput();

};

/**
 * @brief Construct a new Shallow Water object
 * 
 * @param pDT Timestep to be used
 * @param pT  Time to integrate for
 * @param pNx Number of X discretisations
 * @param pNy Number of Y discretisations
 * @param pIc Initial condition to use (1-4)
 */
ShallowWater::ShallowWater(const double& pDT, const double& pT, const int& pNx, const int& pNy, const int& pIc, const int& pMethod) {

    // Initialise default values
    x0 = 0.0;
    y0 = 0.0;

    dx = 1;
    dy = 1;

    // Initialise user inputs
    dt = pDT;
    T = pT;
    nx = pNx;
    ny = pNy;
    method = pMethod;

    // Initialise calculated values
    x1 = x0 + (nx - 1)*dx;
    y1 = y0 + (ny - 1)*dy;

    u = new double[nx*ny];
    v = new double[nx*ny];
    h = new double[nx*ny];

    // Set initial conditions
    SetInitialConditions(pIc);

    // Populate A matrix if using BLAS
    if (method)
    {
        PopulateMatrix();
    }
}

/**
 * @brief Destroy the Shallow Water object
 */
ShallowWater::~ShallowWater() {
    // Deallocate pointers
    delete[] u;
    delete[] v;
    delete[] h;
    delete[] A;
}

/**
 * @brief Set the initial conditions for the solution.
 * 
 * @param ic The index of the initial condition to use
 */
void ShallowWater::SetInitialConditions(const int& ic) {

    // Mean surface height
    double Hm = 10.0;

    double x, y;

    for (int row = 0; row < ny; row++)
    {
        for (int col = 0; col < nx; col++)
        {
            // Calculate local x and y values
            x = x0 + col*dx;
            y = y0 + row*dy;

            // Set both velocity fields to zero
            u[col*ny + row] = 0;
            v[col*ny + row] = 0;

            // Set the height field based on the initial condition chosen
            switch (ic)
            {
            case 1:
                h[col*ny + row] = Hm + exp(-pow((x - 50), 2) / 25);
                break;

            case 2:
                h[col*ny + row] = Hm + exp(-pow((y - 50), 2) / 25);
                break;

            case 3:
                h[col*ny + row] = Hm + exp(-(pow((x - 50), 2) + pow((y - 50), 2)) / 25);
                break;

            case 4:
                h[col*ny + row] = Hm + exp(-(pow((x - 25), 2) + pow((y - 25), 2)) / 25) + exp(-(pow((x - 75), 2) + pow((y - 75), 2)) / 25);
                break;
            }
        }
    }
}

/**
 * @brief Solves the initialised shallow water equation systems.
 */
void ShallowWater::TimeIntegrate() {

    double* k1u = new double[nx*ny];
    double* k1v = new double[nx*ny];
    double* k1h = new double[nx*ny];

    double* k2u = new double[nx*ny];
    double* k2v = new double[nx*ny];
    double* k2h = new double[nx*ny];

    double* k3u = new double[nx*ny];
    double* k3v = new double[nx*ny];
    double* k3h = new double[nx*ny];

    double* k4u = new double[nx*ny];
    double* k4v = new double[nx*ny];
    double* k4h = new double[nx*ny];

    double* tempU = new double[nx*ny];
    double* tempV = new double[nx*ny];
    double* tempH = new double[nx*ny];

    // Total number of timesteps
    int nt = T/dt + 1;

    // Iterate over timesteps and solve the equation
    for (int i = 0; i < nt; i++)
    {
        // Calculate all k1 matrices
        switch (method)
        {
        case 0:
            CalculateFluxLoop(u, v, h, k1u, k1v, k1h);
            break;
        case 1:
            CalculateFluxBLAS(u, v, h, k1u, k1v, k1h);
            break;
        }

        // Calculate y_n + dt*k1/2
        for (int col = 0; col < nx; col++)
        {
            for (int row = 0; row < ny; row++)
            {
                tempU[col*ny + row] = u[col*ny + row] + dt*k1u[col*ny + row]/2;
                tempV[col*ny + row] = v[col*ny + row] + dt*k1v[col*ny + row]/2;
                tempH[col*ny + row] = h[col*ny + row] + dt*k1h[col*ny + row]/2;
            }
        }

        // Calculate k2 matrices
        switch (method)
        {
        case 0:
            CalculateFluxLoop(tempU, tempV, tempH, k2u, k2v, k2h);
            break;
        case 1:
            CalculateFluxBLAS(tempU, tempV, tempH, k2u, k2v, k2h);
            break;
        }

        // Calculate y_n + dt*k2/2
        for (int col = 0; col < nx; col++)
        {
            for (int row = 0; row < ny; row++)
            {
                tempU[col*ny + row] = u[col*ny + row] + dt*k2u[col*ny + row]/2;
                tempV[col*ny + row] = v[col*ny + row] + dt*k2v[col*ny + row]/2;
                tempH[col*ny + row] = h[col*ny + row] + dt*k2h[col*ny + row]/2;
            }
        }

        // Calculate k3 matrices
        switch (method)
        {
        case 0:
            CalculateFluxLoop(tempU, tempV, tempH, k3u, k3v, k3h);
            break;
        case 1:
            CalculateFluxBLAS(tempU, tempV, tempH, k3u, k3v, k3h);
            break;
        }

        // Calculate y_n + dt*k3
        for (int col = 0; col < nx; col++)
        {
            for (int row = 0; row < ny; row++)
            {
                tempU[col*ny + row] = u[col*ny + row] + dt*k3u[col*ny + row];
                tempV[col*ny + row] = v[col*ny + row] + dt*k3v[col*ny + row];
                tempH[col*ny + row] = h[col*ny + row] + dt*k3h[col*ny + row];
            }
        }

        // Calculate k4 matrices
        switch (method)
        {
        case 0:
            CalculateFluxLoop(tempU, tempV, tempH, k4u, k4v, k4h);
            break;
        case 1:
            CalculateFluxBLAS(tempU, tempV, tempH, k4u, k4v, k4h);
            break;
        }

        // Calculate next iteration
        for (int col = 0; col < nx; col++)
        {
            for (int row = 0; row < ny; row++)
            {
                u[col*ny + row] = u[col*ny + row] + dt/6*(k1u[col*ny + row] + 2*k2u[col*ny + row] + 2*k3u[col*ny + row] + k4u[col*ny + row]);
                v[col*ny + row] = v[col*ny + row] + dt/6*(k1v[col*ny + row] + 2*k2v[col*ny + row] + 2*k3v[col*ny + row] + k4v[col*ny + row]);
                h[col*ny + row] = h[col*ny + row] + dt/6*(k1h[col*ny + row] + 2*k2h[col*ny + row] + 2*k3h[col*ny + row] + k4h[col*ny + row]);
            }
        }
    }

    // Output final state to file
    FileOutput();

    delete[] k1u;
    delete[] k1v;
    delete[] k1h;

    delete[] k2u;
    delete[] k2v;
    delete[] k2h;

    delete[] k3u;
    delete[] k3v;
    delete[] k3h;

    delete[] k4u;
    delete[] k4v;
    delete[] k4h;

    delete[] tempU;
    delete[] tempV;
    delete[] tempH;
}

/**
 * @brief Calculates the flux (right hand side) of the shallow water equations using variables pU, pV, pH; puts calculated fluxes in pKU, pKV, pKH.
 * 
 * @param pU  Input matrix for U
 * @param pV  Input matrix for V
 * @param pH  Input matrix for H
 * @param pKU Output matrix for the flux of U
 * @param pKV Output matrix for the flux of V
 * @param pKH Output matrix for the flux of H
 */
void ShallowWater::CalculateFluxLoop(double* pU, double* pV, double* pH, double* pKU, double* pKV, double* pKH) {

    double* hu = new double[nx*ny];
    double* hv = new double[nx*ny];

    double* du_dx = new double[nx*ny];
    double* du_dy = new double[nx*ny];
    double* dh_dx = new double[nx*ny];

    double* dv_dx = new double[nx*ny];
    double* dv_dy = new double[nx*ny];
    double* dh_dy = new double[nx*ny];

    double* dhu_dx = new double[nx*ny];
    double* dhv_dy = new double[nx*ny];

    // Calculate hu and hv
    for (int col = 0; col < nx; col++)
    {
        for (int row = 0; row < ny; row++)
        {
            hu[col*ny + row] = pH[col*ny + row] * pU[col*ny + row];
            hv[col*ny + row] = pH[col*ny + row] * pV[col*ny + row];
        }
    }

    // Calculate all derivatives
    DeriXLoop(pU, du_dx);
    DeriYLoop(pU, du_dy);
    DeriXLoop(pH, dh_dx);

    DeriXLoop(pV, dv_dx);
    DeriYLoop(pV, dv_dy);
    DeriYLoop(pH, dh_dy);

    DeriXLoop(hu, dhu_dx);
    DeriYLoop(hv, dhv_dy);

    // Calculate fluxes
    for (int col = 0; col < nx; col++)
    {
        for (int row = 0; row < ny; row++)
        {
            // Calculate the flux of U
            pKU[col*ny + row] = pU[col*ny + row] * du_dx[col*ny + row] + pV[col*ny + row] * du_dy[col*ny + row] + G*dh_dx[col*ny + row];
            // Calculate the flux of V
            pKV[col*ny + row] = pU[col*ny + row] * dv_dx[col*ny + row] + pV[col*ny + row] * dv_dy[col*ny + row] + G*dh_dy[col*ny + row];
            // Calculate the flux of H
            pKH[col*ny + row] = dhu_dx[col*ny + row] + dhv_dy[col*ny + row];
        }
    }

    delete[] hu;
    delete[] hv;

    delete[] du_dx;
    delete[] du_dy;
    delete[] dh_dx;

    delete[] dv_dx;
    delete[] dv_dy;
    delete[] dh_dy;

    delete[] dhu_dx;
    delete[] dhv_dy;
}

/**
 * @brief Function that calculates the x-derivative of a field VAR and stores it in DER using a loop-based method.
 * 
 * @param var Pointer to the input field matrix
 * @param der Pointer to the output matrix
 */
void ShallowWater::DeriXLoop(const double* var, double* der) {

    for (int row = 0; row < ny; row++)
    {
        // first 3 elements in row
        der[0*ny + row] = (-1.0/60.0*var[(nx-3)*ny+row] + 3.0/20.0*var[(nx-2)*ny+row] - 3.0/4.0*var[(nx-1)*ny+row] + 3.0/4.0*var[ny+row] - 3.0/20.0*var[2*ny+row] + 1.0/60.0*var[3*ny+row])/dx;
        der[1*ny + row] = (-1.0/60.0*var[(nx-2)*ny+row] + 3.0/20.0*var[(nx-1)*ny+row] - 3.0/4.0*var[row] + 3.0/4.0*var[2*ny+row] - 3.0/20.0*var[3*ny+row] + 1.0/60.0*var[4*ny+row])/dx;
        der[2*ny + row] = (-1.0/60.0*var[(nx-1)*ny+row] + 3.0/20.0*var[row] - 3.0/4.0*var[ny+row] + 3.0/4.0*var[3*ny+row] - 3.0/20.0*var[4*ny+row] + 1.0/60.0*var[5*ny+row])/dx;

        // last 3 elements in row
        der[(nx-3)*ny + row] = (-1.0/60.0*var[(nx-6)*ny+row] + 3.0/20.0*var[(nx-5)*ny+row] - 3.0/4.0*var[(nx-4)*ny+row] + 3.0/4.0*var[(nx-2)*ny+row] - 3.0/20.0*var[(nx-1)*ny+row] + 1.0/60.0*var[row])/dx;
        der[(nx-2)*ny + row] = (-1.0/60.0*var[(nx-5)*ny+row] + 3.0/20.0*var[(nx-4)*ny+row] - 3.0/4.0*var[(nx-3)*ny+row] + 3.0/4.0*var[(nx-1)*ny+row] - 3.0/20.0*var[row] + 1.0/60.0*var[ny+row])/dx;
        der[(nx-1)*ny + row] = (-1.0/60.0*var[(nx-4)*ny+row] + 3.0/20.0*var[(nx-3)*ny+row] - 3.0/4.0*var[(nx-2)*ny+row] + 3.0/4.0*var[row] - 3.0/20.0*var[ny+row] + 1.0/60.0*var[2*ny+row])/dx;

        // elements between first 3 and last 3 elements
        for (int col = 3; col < nx-3; col++)
        {
            der[col*ny + row] = (-1.0/60.0*var[(col-3)*ny+row] + 3.0/20.0*var[(col-2)*ny+row] - 3.0/4.0*var[(col-1)*ny+row] + 3.0/4.0*var[(col+1)*ny+row] - 3.0/20.0*var[(col+2)*ny+row] + 1.0/60.0*var[(col+3)*ny+row])/dx;
        }
    }
}

/**
 * @brief Funciton that calculates the y-dericative of a field VAR and stores it in DER using a loop-based method.
 * 
 * @param var Pointer to the input field matrix
 * @param der Pointer to the output matrix
 */
void ShallowWater::DeriYLoop(const double* var, double* der) {

    for (int col = 0; col < nx; col++)
    {
        // first 3 elements in row
        der[col*ny + 0] = (-1.0/60.0*var[col*ny+ny-3] + 3.0/20.0*var[col*ny+ny-2] - 3.0/4.0*var[col*ny+ny-1] + 3.0/4.0*var[col*ny+1] - 3.0/20.0*var[col*ny+2] + 1.0/60.0*var[col*ny+3])/dy;
        der[col*ny + 1] = (-1.0/60.0*var[col*ny+ny-2] + 3.0/20.0*var[col*ny+ny-1] - 3.0/4.0*var[col*ny] + 3.0/4.0*var[col*ny+2] - 3.0/20.0*var[col*ny+3] + 1.0/60.0*var[col*ny+4])/dy;
        der[col*ny + 2] = (-1.0/60.0*var[col*ny+ny-1] + 3.0/20.0*var[col*ny] - 3.0/4.0*var[col*ny+1] + 3.0/4.0*var[col*ny+3] - 3.0/20.0*var[col*ny+4] + 1.0/60.0*var[col*ny+5])/dy;

        // last 3 elements in row
        der[col*ny + ny-3] = (-1.0/60*var[col*ny+ny-6] + 3.0/20.0*var[col*ny+ny-5] - 3.0/4.0*var[col*ny+ny-4] + 3.0/4.0*var[col*ny+ny-2] - 3.0/20.0*var[col*ny+ny-1] + 1.0/60.0*var[col*ny])/dy;
        der[col*ny + ny-2] = (-1.0/60*var[col*ny+ny-5] + 3.0/20.0*var[col*ny+ny-4] - 3.0/4.0*var[col*ny+ny-3] + 3.0/4.0*var[col*ny+ny-1] - 3.0/20.0*var[col*ny] + 1.0/60.0*var[col*ny+1])/dy;
        der[col*ny + ny-1] = (-1.0/60*var[col*ny+ny-4] + 3.0/20.0*var[col*ny+ny-3] - 3.0/4.0*var[col*ny+ny-2] + 3.0/4.0*var[col*ny] - 3.0/20.0*var[col*ny+1] + 1.0/60.0*var[col*ny+2])/dy;

        // elements between first 3 and last 3 elements
        for (int row = 3; row < ny-3; row++)
        {
            der[col*ny + row] = (-1.0/60.0*var[col*ny+row-3] + 3.0/20.0*var[col*ny+row-2] - 3.0/4.0*var[col*ny+row-1] + 3.0/4.0*var[col*ny+row+1] - 3.0/20.0*var[col*ny+row+2] + 1.0/60.0*var[col*ny+row+3])/dy;
        }
    }
}

/**
 * @brief Calculates the flux (right hand side) of the shallow water equations using variables pU, pV, pH; puts calculated fluxes in pKU, pKV, pKH.
 * 
 * @param pU  Input matrix for U
 * @param pV  Input matrix for V
 * @param pH  Input matrix for H
 * @param pKU Output matrix for U
 * @param pKV Output matrix for V
 * @param pKH Output matrix for h
 */
void ShallowWater::CalculateFluxBLAS(double* pU, double* pV, double* pH, double* pKU, double* pKV, double* pKH) {

    double* hu = new double[nx*ny];
    double* hv = new double[nx*ny];

    double* du_dx = new double[nx*ny];
    double* du_dy = new double[nx*ny];
    double* dh_dx = new double[nx*ny];

    double* dv_dx = new double[nx*ny];
    double* dv_dy = new double[nx*ny];
    double* dh_dy = new double[nx*ny];

    double* dhu_dx = new double[nx*ny];
    double* dhv_dy = new double[nx*ny];

    // Calculate hu and hv
    cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, 1.0, pH, 1, pU, 1, 0.0, hu, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, 1.0, pH, 1, pV, 1, 0.0, hv, 1);

    // Calculate all derivatives

    DeriXLoop(pU, du_dx);
    DeriYLoop(pU, du_dy);
    DeriXLoop(pH, dh_dx);

    DeriXLoop(pV, dv_dx);
    DeriYLoop(pV, dv_dy);
    DeriYLoop(pH, dh_dy);

    DeriXLoop(hu, dhu_dx);
    DeriYLoop(hv, dhv_dy);

    // Calculate fluxes
    // u
    cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, 1.0, pU, 1, du_dx, 1, 0.0, pKU, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, 1.0, pV, 1, du_dy, 1, 1.0, pKU, 1);
    cblas_daxpy(nx*ny, G, dh_dx, 1, pKU, 1);

    // v
    cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, 1.0, pU, 1, dv_dx, 1, 0.0, pKV, 1);
    cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, 1.0, pV, 1, dv_dy, 1, 1.0, pKV, 1);
    cblas_daxpy(nx*ny, G, dh_dy, 1, pKV, 1);

    // h
    cblas_daxpy(nx*ny, 1.0, dhu_dx, 1, dhv_dy, 1);
    cblas_dcopy(nx*ny, dhv_dy, 1, pKH, 1);

    delete[] hu;
    delete[] hv;

    delete[] du_dx;
    delete[] du_dy;
    delete[] dh_dx;

    delete[] dv_dx;
    delete[] dv_dy;
    delete[] dh_dy;

    delete[] dhu_dx;
    delete[] dhv_dy;
}

/**
 * @brief Populate the A matrix for BLAS based derivatives
 * 
 */
void ShallowWater::PopulateMatrix() {

    int n;

    if (nx > ny)
    {
        n = nx;
    } else {
        n = ny;
    }

    A = new double[n*n];

    for (int col = 0; col < n; col++)
    {
        A[col*7 + 0] = -1.0/60.0;
        A[col*7 + 1] = 3.0/20.0;
        A[col*7 + 2] = -3.0/4.0;
        A[col*7 + 3] = 0.0;
        A[col*7 + 4] = 3.0/4.0;
        A[col*7 + 5] = -3.0/20.0;
        A[col*7 + 6] = -1.0/60.0; 
    }
}

/**
 * @brief Outputs the current state of the solution to an output file
 */
void ShallowWater::FileOutput() {

    ofstream out("output.txt", ios::out | ios::trunc);

    out.precision(5);

    double x, y;

    for (int row = 0; row < ny; row++)
    {
        for (int col = 0; col < nx; col++)
        {
            x = x0 + col*dx;
            y = y0 + row*dy;
            // Output x and y location
            out << x << " " << y << " ";

            // Output variables
            out << u[col*ny + row] << " " << v[col*ny + row] << " " << h[col*ny + row] << endl;
        }

        // Output empty line after column
        out << endl;
    }

    out.close();
}

int main(int argc, char* argv[]) {

    po::options_description opts("Solves the shallow water equations");

    opts.add_options()
        ("help,h", "Print help message")
        ("dT", po::value<double>()->default_value(0.1), "Timestep to use")
        ("T", po::value<double>()->default_value(80), "Time to integrate for")
        ("Nx", po::value<int>()->default_value(100), "Number of points in X direction")
        ("Ny", po::value<int>()->default_value(100), "NUmber of points in Y direction")
        ("method", po::value<int>()->default_value(0), "Method to use. 0 for LOOP-based and 1 for BLAS-based.")
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
    int ic;
    int method;

    if (vm.count("ic")){

        ic = vm["ic"].as<int>();

        if (ic < 1 || 4 < ic)
        {
            cout << "Initial conditions not valid. Valid range is 1-4. Exiting program." << endl;
            return -1;
        }

    } else {
        cout << "Initial conditions not specified. Exiting program." << endl;
        return -1;
    }

    if (vm.count("method")){

        method = vm["method"].as<int>();

        if (method < 0 || 1 < method)
        {
            cout << "Method not valid. Valid range is 0-1. Exiting program." << endl;
            return -1;
        }

    } else {
        cout << "Method not specified. Exiting program." << endl;
        return -1;
    }

    // Create solver object
    ShallowWater solver = ShallowWater(dT, T, Nx, Ny, ic, method);

    // Solve system
    solver.TimeIntegrate();
}