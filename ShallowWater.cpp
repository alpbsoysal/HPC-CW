/**
 * @file ShallowWater.cpp
 * @author Alp Soysal
 * @brief A program that solves the shallow water equations using 6th order central differencing and 4th order runge-kutta, for the HPC module coursework.
 * @version 1.0
 * @date 2023-03-19
 */

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cmath>

#include <boost/program_options.hpp>
#include <cblas.h>
#include <omp.h>

using namespace std;
namespace po = boost::program_options;

// Set gravity compile constant
constexpr double G = 9.81;

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
        double c1, c2, c3, c4, c5, c6, ddx, ddy;
        double *u, *v, *h, *A, *Aut, *Alt;

    public:

        ShallowWater(const double& pDT, const double& pT, const int& pNx, const int& pNy, const int& pIc, const int& pMethod);
        ~ShallowWater();

        void SetInitialConditions(const int& ic);
        void TimeIntegrate();

        void CalculateFluxLoop(const double* pU, const double* pV, const double* pH, double* pKU, double* pKV, double* pKH);
        void DeriXLoop(const double* var, double* der);
        void DeriYLoop(const double* var, double* der);

        void CalculateFluxBLAS(const double* pU, const double* pV, const double* pH, double* pKU, double* pKV, double* pKH);
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

    c1 = -1.0/60.0;
    c2 =  3.0/20.0;
    c3 = -3.0/4.0;
    c4 =  3.0/4.0;
    c5 = -3.0/20.0;
    c6 =  1.0/60.0;
    ddx = 1.0/dx;
    ddy = 1.0/dy;

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
    if(method)
    {
        delete[] A;
        delete[] Aut;
        delete[] Alt;
    }
}

/**
 * @brief Set the initial conditions for the solution.
 * 
 * @param ic The index of the initial condition to use
 */
void ShallowWater::SetInitialConditions(const int& ic) {

    // Mean surface height
    double Hm = 10.0;

    // Declare private variables
    double x, y, val1, val2;

    // Loop inside switch for efficiency
    switch (ic)
    {
    case 1:
        #pragma omp parallel for default(shared) private(x, y, val1, val2) schedule(static)
        for (int col = 0; col < nx; col++)
        {
            // Precompute possible values
            x = x0 + col*dx;
            val1 = Hm + exp(-(x - 50)*(x - 50)*0.04);

            for (int row = 0; row < ny; row++)
            {
                h[col*ny + row] = val1;
            }
        }
        break;

    case 2:
        #pragma omp parallel for default(shared) private(x, y) schedule(static)
        for (int row = 0; row < ny; row++)
        {
            y = y0 + row*dy;
            val1 = Hm + exp(-(y - 50)*(y - 50)*0.04);

            for (int col = 0; col < nx; col++)
            {
                h[col*ny + row] = val1;
            }
        }
        break;

    case 3:
        #pragma omp parallel for default(shared) private(x, y) schedule(static)
        for (int col = 0; col < nx; col++)
        {
            x = x0 + col*dx;
            val1 = ((x - 50)*(x - 50));

            for (int row = 0; row < ny; row++)
            {
                y = y0 + row*dy;
                h[col*ny + row] = Hm + exp(-(val1 + (y - 50)*(y - 50))*0.04);
            }
        }
        break;

    case 4:
        #pragma omp parallel for default(shared) private(x, y) schedule(static)
        for (int col = 0; col < nx; col++)
        {
            x = x0 + col*dx;
            val1 = ((x - 25)*(x - 25));
            val2 = ((x - 75)*(x - 75));

            for (int row = 0; row < ny; row++)
            {
                y = y0 + row*dy;
                h[col*ny + row] = Hm + exp(-(val1 + (y - 25)*(y - 25))*0.04) + exp(-(val2 + (y - 75)*(y - 75))*0.04);
            }
        }
        break;
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
    int nt = T/dt;

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
        #pragma omp for schedule(static)
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
        #pragma omp for schedule(static)
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
        #pragma omp for schedule(static)
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
        #pragma omp for schedule(static)
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
void ShallowWater::CalculateFluxLoop(const double* pU, const double* pV, const double* pH, double* pKU, double* pKV, double* pKH) {

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

    // Calculate all derivatives
    DeriXLoop(pU, du_dx);
    DeriYLoop(pU, du_dy);
    DeriXLoop(pH, dh_dx);

    DeriXLoop(pV, dv_dx);
    DeriYLoop(pV, dv_dy);
    DeriYLoop(pH, dh_dy);

    // Calculate fluxes
    #pragma omp parallel for schedule(static)
    for (int i = 0; i < nx*ny; i++)
    {
        // Calculate hu and hv using product rule
        dhu_dx[i] = pH[i]*du_dx[i] + pU[i]*dh_dx[i];
        dhv_dy[i] = pH[i]*dv_dy[i] + pV[i]*dh_dy[i];

        // Calculate the flux of U
        pKU[i] = - pU[i] * du_dx[i] - pV[i] * du_dy[i] - G*dh_dx[i];
        // Calculate the flux of V
        pKV[i] = - pU[i] * dv_dx[i] - pV[i] * dv_dy[i] - G*dh_dy[i];
        // Calculate the flux of H
        pKH[i] = - dhu_dx[i] - dhv_dy[i];
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

    int index;

    #pragma omp parallel for private(index) default(shared) schedule(static)
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

    int index, colindex;

    #pragma omp parallel for private(index, colindex) default(shared) schedule(static)
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
void ShallowWater::CalculateFluxBLAS(const double* pU, const double* pV, const double* pH, double* pKU, double* pKV, double* pKH) {

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

    // Calculate all derivatives
    DeriXBLAS(pU, du_dx);
    DeriYBLAS(pU, du_dy);
    DeriXBLAS(pH, dh_dx);
 
    DeriXBLAS(pV, dv_dx);
    DeriYBLAS(pV, dv_dy);
    DeriYBLAS(pH, dh_dy);

    // Calculate fluxes
    #pragma omp parallel sections default(shared)
    {
        // u
        #pragma omp section
        {
            cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, -1.0, pU, 1, du_dx, 1, 0.0, pKU, 1);
            cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, -1.0, pV, 1, du_dy, 1, 1.0, pKU, 1);
            cblas_daxpy(nx*ny, -G, dh_dx, 1, pKU, 1);
        }
        // v
        #pragma omp section
        {
            cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, -1.0, pU, 1, dv_dx, 1, 0.0, pKV, 1);
            cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, -1.0, pV, 1, dv_dy, 1, 1.0, pKV, 1);
            cblas_daxpy(nx*ny, -G, dh_dy, 1, pKV, 1);
        }
        // h using product rule
        #pragma omp section
        {
            cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, -1.0, pH, 1, du_dx, 1, 0.0, pKH, 1);
            cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, -1.0, pU, 1, dh_dx, 1, 1.0, pKH, 1);
            cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, -1.0, pH, 1, dv_dy, 1, 1.0, pKH, 1);
            cblas_dgbmv(CblasColMajor, CblasNoTrans, nx*ny, nx*ny, 0, 0, -1.0, pV, 1, dh_dy, 1, 1.0, pKH, 1);
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
 * @brief Populate the A matrix for BLAS based derivatives
 */
void ShallowWater::PopulateMatrix() {

    int n;

    // Choose largest width, creating one stencil Matrix is sufficient
    if (nx > ny)
    {
        n = nx;
    } else {
        n = ny;
    }

    A = new double[7*n];

    // Populate A using banded storage
    #pragma omp parallel for default(shared) schedule(static)
    for (int col = 0; col < n; col++)
    {
        A[col*7 + 0] = c6;
        A[col*7 + 1] = c5;
        A[col*7 + 2] = c4;
        A[col*7 + 3] = 0.0;
        A[col*7 + 4] = c3;
        A[col*7 + 5] = c2;
        A[col*7 + 6] = c1; 
    }

    // For boundary conditions, populate packed triangular 3x3 matrices
    // Upper right triangle
    Aut = new double[6];
    Aut[0] = c1;
    Aut[1] = c2;
    Aut[2] = c1;
    Aut[3] = c3;
    Aut[4] = c2;
    Aut[5] = c1;

    // Lower left triangle
    Alt = new double[6];
    Alt[0] = c6;
    Alt[1] = c5;
    Alt[2] = c4;
    Alt[3] = c6;
    Alt[4] = c5;
    Alt[5] = c6;
}

/**
 * @brief Function that calculates the x-derivative of a field VAR and stores it in DER using a matrix-based method.
 * 
 * @param var Pointer to the input field matrix
 * @param der Pointer to the output matrix
 */
void ShallowWater::DeriXBLAS(const double* var, double* der) {

    // Allocate temporary storage arrays
    double* tvar = new double[nx*ny];
    double* bc;

    // Copy var as BLAS modifies it
    cblas_dcopy(nx*ny, var, 1, tvar, 1);

    #pragma omp parallel private(bc) default(shared)
    {
        bc = new double[6];

        #pragma omp  for schedule(static)
        for (int row = 0; row < ny; row++)
        {
            // Copy first 3 and last 3 elements of each row into boundary condition array
            cblas_dcopy(3, tvar+row, ny, bc, 1);
            cblas_dcopy(3, tvar+ny*(nx-3)+row, ny, bc+3, 1);

            // Multiply stencil with row of input field. Row accessed using pointer arithmethic
            cblas_dgbmv(CblasColMajor, CblasNoTrans, nx, nx, 3, 3, 1.0, A, 7, tvar+row, ny, 0, der+row, ny);

            // Calculate contributions due to periodic boundary conditions and store back in BC
            cblas_dtpmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, 3, Alt, bc, 1);
            cblas_dtpmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 3, Aut, bc+3, 1);

            // Add boundary condition contributions to derivative matrixS
            cblas_daxpy(3, 1.0, bc, 1, der+ny*(nx-3)+row, ny);
            cblas_daxpy(3, 1.0, bc+3, 1, der+row, ny);
        }
    }
}

/**
 * @brief Function that calculates the y-derivative of a field VAR and stores it in DER using a matrix-based method.
 * 
 * @param var Pointer to the input field matrix
 * @param der Pointer to the output matrix
 */
void ShallowWater::DeriYBLAS(const double* var, double* der) {

    // Allocate temporary storage arrays
    double* tvar = new double[nx*ny];
    double* bc;

    // Copy var as BLAS modifies it
    cblas_dcopy(nx*ny, var, 1, tvar, 1);

    #pragma omp parallel private(bc) default(shared)
    {
        bc = new double[6];

        #pragma omp for schedule(static)
        for (int col = 0; col < nx; col++)
        {
            // Copy first 3 and last 3 elements of each column into boundary condition array
            cblas_dcopy(3, tvar+col*ny, 1, bc, 1);
            cblas_dcopy(3, tvar+col*ny+ny-3, 1, bc+3, 1);

            // Multiply stencil with column of input field. Column accessed using pointer arithmethic
            cblas_dgbmv(CblasColMajor, CblasNoTrans, ny, ny, 3, 3, 1.0, A, 7, tvar+col*ny, 1, 0, der+col*ny, 1);

            // Calculate contributions due to periodic boundary conditions and store back in BC
            cblas_dtpmv(CblasColMajor, CblasLower, CblasNoTrans, CblasNonUnit, 3, Alt, bc, 1);
            cblas_dtpmv(CblasColMajor, CblasUpper, CblasNoTrans, CblasNonUnit, 3, Aut, bc+3, 1);

            // Add boundary condition contributions to derivative matrixS
            cblas_daxpy(3, 1.0, bc, 1, der+col*ny+ny-3, 1);
            cblas_daxpy(3, 1.0, bc+3, 1, der+col*ny, 1);
        }
    }
}

/**
 * @brief Outputs the current state of the solution to an output file
 */
void ShallowWater::FileOutput() {

    // Open file
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
            out << u[col*ny + row] << " " << v[col*ny + row] << " " << h[col*ny + row] << "\n";
        }

        // Output empty line after column
        out << "\n";
    }

    // Close file
    out.close();
}

int main(int argc, char* argv[]) {

    po::options_description opts("Solves the shallow water equations");

    // Define command line options
    opts.add_options()
        ("help,h", "Print help message")
        ("dt", po::value<double>()->default_value(0.1), "Timestep to use")
        ("T", po::value<double>()->default_value(80), "Time to integrate for")
        ("Nx", po::value<int>()->default_value(100), "Number of points in X direction")
        ("Ny", po::value<int>()->default_value(100), "NUmber of points in Y direction")
        ("method", po::value<int>()->default_value(0), "Method to use. 0 for LOOP-based and 1 for BLAS-based.")
        ("ic", po::value<int>(), "Initial condition to use, 1-4");

    // Parse command line options
    po::variables_map vm;
    po::store(po::parse_command_line(argc, argv, opts), vm);
    po::notify(vm);

    // Produce output text
    if (vm.count("help"))
    {
        cout << opts << endl;
        return 0;
    }

    // Store inputs
    const double dT = vm["dt"].as<double>();
    const double T  = vm["T"].as<double>();
    const int Nx = vm["Nx"].as<int>();
    const int Ny = vm["Ny"].as<int>();
    int ic;
    int method;

    // Checking if --ic is valid
    if (vm.count("ic"))
    {
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

    // Checking if --method is valid
    if (vm.count("method"))
    {
        method = vm["method"].as<int>();

        if (method < 0 || 1 < method)
        {
            cout << "Method not valid. Valid range is 0-1. Exiting program." << endl;
            return -1;
        }
    }

    // Create solver object
    ShallowWater solver = ShallowWater(dT, T, Nx, Ny, ic, method);

    // Solve system
    solver.TimeIntegrate();
}
