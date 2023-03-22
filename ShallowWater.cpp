/**
 * @file ShallowWater.cpp
 * @author Alp Soysal
 * @brief A program that solves the shallow water equations using 6th order central differencing and 4th order runge-kutta, for the HPC module coursework.
 * @version 0.2
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
        double *u, *v, *h, *Ax, *Ay;

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
    delete[] Ax;
    delete[] Ay;
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
        for (int i = 0; i < nx*ny; i++)
        {

            tempU[i] = u[i] + dt*k1u[i]*0.5;
            tempV[i] = v[i] + dt*k1v[i]*0.5;
            tempH[i] = h[i] + dt*k1h[i]*0.5;
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
        for (int i = 0; i < nx*ny; i++)
        {
            tempU[i] = u[i] + dt*k2u[i]*0.5;
            tempV[i] = v[i] + dt*k2v[i]*0.5;
            tempH[i] = h[i] + dt*k2h[i]*0.5;
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
        for (int i = 0; i < nx*ny; i++)
        {
            tempU[i] = u[i] + dt*k3u[i];
            tempV[i] = v[i] + dt*k3v[i];
            tempH[i] = h[i] + dt*k3h[i];
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
        for (int i = 0; i < nx*ny; i++)
        {
                u[i] = u[i] + dt/6*(k1u[i] + 2*k2u[i] + 2*k3u[i] + k4u[i]);
                v[i] = v[i] + dt/6*(k1v[i] + 2*k2v[i] + 2*k3v[i] + k4v[i]);
                h[i] = h[i] + dt/6*(k1h[i] + 2*k2h[i] + 2*k3h[i] + k4h[i]);
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
    for (int i = 0; i < nx*ny; i++)
    {
        hu[i] = pH[i] * pU[i];
        hv[i] = pH[i] * pV[i];
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
    for (int i = 0; i < nx*ny; i++)
    {
        // Calculate the flux of U
        pKU[i] = pU[i] * du_dx[i] + pV[i] * du_dy[i] + G*dh_dx[i];
        // Calculate the flux of V
        pKV[i] = pU[i] * dv_dx[i] + pV[i] * dv_dy[i] + G*dh_dy[i];
        // Calculate the flux of H
        pKH[i] = dhu_dx[i] + dhv_dy[i];
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

    double c1 = -1.0/60.0;
    double c2 =  3.0/20.0;
    double c3 = -3.0/4.0;
    double c4 =  3.0/4.0;
    double c5 = -3.0/20.0;
    double c6 =  1.0/60.0;
    double ddx = 1.0/dx;
    int index;

    for (int row = 0; row < ny; row++)
    {
        // first 3 elements in row
        der[row] = (c1*var[(nx-3)*ny+row] + c2*var[(nx-2)*ny+row] + c3*var[(nx-1)*ny+row] + c4*var[ny+row] + c5*var[2*ny+row] + c6*var[3*ny+row])*ddx;
        der[1*ny + row] = (c1*var[(nx-2)*ny+row] + c2*var[(nx-1)*ny+row] + c3*var[row] + c4*var[2*ny+row] + c5*var[3*ny+row] + c6*var[4*ny+row])*ddx;
        der[2*ny + row] = (c1*var[(nx-1)*ny+row] + c2*var[row] + c3*var[ny+row] + c4*var[3*ny+row] + c5*var[4*ny+row] + c6*var[5*ny+row])*ddx;

        // elements between first 3 and last 3 elements
        for (int col = 3; col < nx-3; col++)
        {
            index = col*ny+row;
            der[index] = (c1*var[index-3*ny] + c2*var[index-2*ny] + c3*var[index-ny] + c4*var[index+ny] + c5*var[index+2*ny] + c6*var[index+3*ny])*ddx;
        }

        // last 3 elements in row
        der[(nx-3)*ny + row] = (c1*var[(nx-6)*ny+row] + c2*var[(nx-5)*ny+row] + c3*var[(nx-4)*ny+row] + c4*var[(nx-2)*ny+row] + c5*var[(nx-1)*ny+row] + c6*var[row])*ddx;
        der[(nx-2)*ny + row] = (c1*var[(nx-5)*ny+row] + c2*var[(nx-4)*ny+row] + c3*var[(nx-3)*ny+row] + c4*var[(nx-1)*ny+row] + c5*var[row] + c6*var[ny+row])*ddx;
        der[(nx-1)*ny + row] = (c1*var[(nx-4)*ny+row] + c2*var[(nx-3)*ny+row] + c3*var[(nx-2)*ny+row] + c4*var[row] + c5*var[ny+row] + c6*var[2*ny+row])*ddx;
    }
}

/**
 * @brief Funciton that calculates the y-dericative of a field VAR and stores it in DER using a loop-based method.
 * 
 * @param var Pointer to the input field matrix
 * @param der Pointer to the output matrix
 */
void ShallowWater::DeriYLoop(const double* var, double* der) {

    double c1 = -1.0/60.0;
    double c2 =  3.0/20.0;
    double c3 = -3.0/4.0;
    double c4 =  3.0/4.0;
    double c5 = -3.0/20.0;
    double c6 =  1.0/60.0;
    double ddy = 1.0/dy;
    int index, colindex;

    for (int col = 0; col < nx; col++)
    {
        colindex = col*ny;
        // first 3 elements in row
        der[colindex + 0] = (c1*var[colindex+ny-3] + c2*var[colindex+ny-2] + c3*var[colindex+ny-1] + c4*var[colindex+1] + c5*var[colindex+2] + c6*var[colindex+3])*ddy;
        der[colindex + 1] = (c1*var[colindex+ny-2] + c2*var[colindex+ny-1] + c3*var[colindex] + c4*var[colindex+2] + c5*var[colindex+3] + c6*var[colindex+4])*ddy;
        der[colindex + 2] = (c1*var[colindex+ny-1] + c2*var[colindex] + c3*var[colindex+1] + c4*var[colindex+3] + c5*var[colindex+4] + c6*var[colindex+5])*ddy;

        // elements between first 3 and last 3 elements
        for (int row = 3; row < ny-3; row++)
        {
            index = colindex+row;
            der[index] = (c1*var[index-3] + c2*var[index-2] + c3*var[index-1] + c4*var[index+1] + c5*var[index+2] + c6*var[index+3])*ddy;
        }

        // last 3 elements in row
        der[colindex + ny-3] = (c1*var[colindex+ny-6] + c2*var[colindex+ny-5] + c3*var[colindex+ny-4] + c4*var[colindex+ny-2] + c5*var[colindex+ny-1] + c6*var[colindex])*ddy;
        der[colindex + ny-2] = (c1*var[colindex+ny-5] + c2*var[colindex+ny-4] + c3*var[colindex+ny-3] + c4*var[colindex+ny-1] + c5*var[colindex] + c6*var[colindex+1])*ddy;
        der[colindex + ny-1] = (c1*var[colindex+ny-4] + c2*var[colindex+ny-3] + c3*var[colindex+ny-2] + c4*var[colindex] + c5*var[colindex+1] + c6*var[colindex+2])*ddy;
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

    DeriXBLAS(pU, du_dx);
    DeriYBLAS(pU, du_dy);
    DeriXBLAS(pH, dh_dx);

    DeriXBLAS(pV, dv_dx);
    DeriYBLAS(pV, dv_dy);
    DeriYBLAS(pH, dh_dy);

    DeriXBLAS(hu, dhu_dx);
    DeriYBLAS(hv, dhv_dy);

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
 */
void ShallowWater::PopulateMatrix() {

    Ax = new double[nx*nx]();
    Ay = new double[ny*ny]();

    double* Arow = new double[ny]();
    double* Acol = new double[nx]();

    Arow[1] = Acol[1] = -3.0/4.0;
    Arow[2] = Acol[2] = 3.0/20.0;
    Arow[3] = Acol[3] = -1.0/60.0;
    Arow[ny-3] = Acol[nx-3] = 1.0/60.0;
    Arow[ny-2] = Acol[nx-2] = -3.0/20.0;
    Arow[ny-1] = Acol[nx-1] = 3.0/4.0;

    for (int row = 0; row < ny; row++)
    {
        for (int col = 0; col < ny; col++)
        {
            Ay[col*ny+row] = Arow[(ny-col+row)%ny];
        }
    }

    for (int row = 0; row < nx; row++)
    {
        for (int col = 0; col < nx; col++)
        {
            Ax[col*nx+row] = Acol[(nx-col+row)%nx];
        }
    }
}

/**
 * @brief Function that calculates the x-derivative of a field VAR and stores it in DER using a matrix-based method.
 * 
 * @param var Pointer to the input field matrix
 * @param der Pointer to the output matrix
 */
void ShallowWater::DeriXBLAS(const double* var, double* der) {

    // Calculate derivative using matrix-matrix multiplication
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ny, nx, nx, 1.0, var, ny, Ax, nx, 0.0, der, ny);

}

/**
 * @brief Function that calculates the y-derivative of a field VAR and stores it in DER using a matrix-based method.
 * 
 * @param var Pointer to the input field matrix
 * @param der Pointer to the output matrix
 */
void ShallowWater::DeriYBLAS(const double* var, double* der) {

    // Calculate derivative using matrix-matrix multiplication
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, ny, nx, ny, 1.0, Ay, ny, var, ny, 0.0, der, ny);
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