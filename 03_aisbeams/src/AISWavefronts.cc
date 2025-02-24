#include "AISWavefronts.hh"

// flat wavefront and gradient wavefront
double flatWavefront(const std::array<double, 3>& pos, const double& param1, const double& param2)
{
    return 0.0;
}

std::array<double, 3> flatGradientWavefront(const std::array<double, 3>& pos, const double& paaram1, const double& param2)
{
    return {0.0, 0.0, 0.0};
}

// gaussian wavefront and gradient
double gaussianWavefront(const std::array<double, 3>& pos, const double& kz,const double& w0)
{
    double r   = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    double z   = pos[2];
    double zR  = pi * w0 * w0 / lambdaSr87;
    double RInv = z / (z*z + zR*zR);

    // if the wavevector is negative, invert the sign phase
    double sign;
    if(kz < 0)
    {
        sign = -1.0;
    }
    else
    {
        sign = 1.0;
    }

    double phase = abs(kz) * r * r * RInv / 2.0 - atan(z/zR);

    return sign*phase *0; // wavefront curvature accounted for in the zernike coeffs in this branch
}

std::array<double,3> gaussianGradientWavefront(const std::array<double, 3>& pos, const double& kz, const double& w0)
{
    double r   = sqrt(pos[0]*pos[0] + pos[1]*pos[1]);
    double x   = pos[0];
    double y   = pos[1];
    double z   = pos[2];
    double zR  = pi * w0 * w0 / lambdaSr87;
    double RInv = z / (z*z + zR*zR);

    double sign;
    if(kz < 0)
    {
        sign = -1.0; // invert the sign of the gouy phase if the wavevector is negative (see gaussianWavefront)
    }
    else
    {
        sign = 1.0;
    }

    double dWdX = abs(kz) * x * RInv;
    double dWdY = abs(kz) * y * RInv;
    double dWdZ = - 1.0 / (zR*(1 + (z/zR)*(z/zR))) - abs(kz) * r * r * RInv * RInv + abs(kz) * r * r * RInv / (2.0*z);

    std::array<double,3> gradient = {sign * dWdX, sign * dWdY, sign * dWdZ};

    return {0.0, 0.0, 0.0}; // wavefront curvature accounted for in the zernike coeffs in this branch   
}

double Rmn(int m, int n, double rho)
{
    // note: rho is rho/R
    double res;

    if(m==n)
    {
        res = pow(rho, n);
    }
    else if(m == n-2)
    {
        res = n*Rmn(n,n,rho) + (n-1)*Rmn(n-2,n-2,rho);
    }
    else
    {
        double K1 = 4*n*(n-1)/((n+m)*(n-m));
        double K2 = 2*((m*m*(n-1)+n*(n-1)*(n-2)))/((n+m)*(m-n)*(n-2));
        double K3 = -n*(m+n-2)*(m-n+2)/((m+n)*(m-n)*(n-2));

        res = (K1 * rho * rho + K2) * Rmn(m, n-2, rho) + K3 * Rmn(m, n-4, rho);
    }

    return res;
}
double Zmn(int m, int n, double rho, double theta)
{
    // note: rho is rho/R
    if(m == 0)
    {
        return Rmn(m, n, rho);
    }
    else if(m > 0)
    {
        return Rmn(m, n, rho) * cos(m*theta);
    }
    else
    {
        return Rmn(-m, n, rho) * sin(-m*theta);
    }
}