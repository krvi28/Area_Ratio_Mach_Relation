// Newton Method to solve the area-mach number relation iteratively

#include <iostream>
#include <iomanip>
#include <cmath>

double gamma {1.196};          // specific heat ratio
double area_mach_relation(double &, double &);

int main()
{   
    std::cout << std::fixed << std::setprecision(4);
    double exit_mach, mach_est, area_ratio;
    area_ratio = 77.5;          // A/A*
    // Assume Mach number for initialization
    // Do not start at 0 or 1
    // If the expected mach is < 1, initialize at M slightly > 0
    // If the expected mach is > 1, initialize at M > 3
    mach_est = 3.0;            

    exit_mach = area_mach_relation(area_ratio, mach_est);       // function for solving non-linear equation

    std::cout << "Exit Mach Number is " << exit_mach << std::endl;
}

// Area expansion ratio vs exit Mach number relation (non-linear)
double area_mach_relation(double& area_exp_ratio, double& mach_est) {
    double F_M {}, Deriv_M {}, error {100}, allow_err {0.001};
    double M {mach_est}, M_sq;
    double A, B, C, D, E, F;
    
    while (error >= allow_err) {
        M_sq = pow(M, 2);
        A = 2.0 / (gamma + 1.0);
        B = 1.0 + ((gamma - 1.0) * M_sq * 0.5);
        C = (gamma + 1.0) / (2.0 * (gamma - 1.0));

        D = pow(2, (1.0 - (3.0 * gamma)) / (2.0 - (2.0 * gamma)));
        E = M_sq * (2.0 + (M_sq * (gamma - 1)));
        F = M_sq - 1;

        F_M = ((1.0 / M) * pow(A * B, C)) - area_exp_ratio;
        Deriv_M = D * (F / E) * pow((B / (gamma + 1.0)), C);
        M = M - (F_M / Deriv_M);

        error = abs((1.0 / M) * pow(A * B, C) - area_exp_ratio) / area_exp_ratio;

    }
    return M;
}
