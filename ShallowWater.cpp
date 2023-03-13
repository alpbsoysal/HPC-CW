
class ShallowWater {

    private:
        double* u, v, h;

    public:
        void SetInitialConditions() {
            return;
        }

        // Implement runge-kutta 4
        void TimeIntegrate() {
            // calls CalculateFlux_loop() or CalculateFlux_matrix()
            return;
        }

        // Calculate the fluxes for all variables
        void CalculateFlux_loop() {
            // calls deri_x and deri_y
            return;
        }

        void CalculateFlux_blas() {
            return;
        }

        // Calculate an x derivative using a loop-based method
        void deri_x(double* var, int nx, int ny, double dx, double* der) {

        }

        void deri_y(double* var, int nx, int ny, double dy, double* der) {

        }
};