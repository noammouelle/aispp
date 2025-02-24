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
    else if(n - m == 2)
    {
        res = ((m+2)*pow(rho,2) - (m+1))*pow(rho,m);
    }
    else
    {
        double K1 = 2*(n-2)*((n-2+m)/2 + 1)*((n-2-m)/2 + 1);
        double K2 = 2*(n-2)*(n-1)*n;
        double K3 = -pow(m,2)*(n-1) - n*(n-1)*(n-2);
        double K4 = -2*n*((n+m-2)/2)*((n-m-2)/2);

        res = ((K2 * pow(rho,2) + K3) * Rmn(m, n-2, rho) + K4 * Rmn(m, n-4, rho))/K1;
    }

    return res;
}
double Zmn(int m, int n, double rho, double theta)
{
    // note: rho is rho/R

    double prefactor = sqrt((2.0 * n + 2.0) / (1.0 + (m == 0)));
    double res;
    if(m == 0)
    {
        res = Rmn(m, n, rho);
    }
    else if(m > 0)
    {
        res = Rmn(m, n, rho) * cos(m*theta);
    }
    else
    {
        res = Rmn(-m, n, rho) * sin(-m*theta);
    }

    return prefactor * res;
}

double a(int n, int m, int nprime, int mprime)
{
    double res;

    res = sqrt((2-(m==0))/(2-(mprime==0))*(n+1)*(nprime+1));

    return res;
}

std::array<double,3> GradZmn(int m, int n, double rho, double theta)
{
    double alpham;
    if(m >= 0)
    {
        alpham = 1;
    }
    else
    {
        alpham = -1;
    }

    // compute the x-component of the gradient
    double gradX = 0.0;
    for(int nprime = 0; nprime <= n-1; nprime++)
    {
        if((n >=abs(m)) && (nprime >= abs(m-1)) && ((n-abs(m))%2 == 0) && ((nprime-abs(m-1))%2 == 0))
        {
            gradX += a(n, m, nprime, m-1) * Zmn(alpham*abs(m-1), nprime, rho, theta);
        }
        if((n >= abs(m)) && (nprime >= abs(m+1)) && ((n-abs(m))%2 == 0) && ((nprime-abs(m+1))%2 == 0))
        {
            gradX += alpham*sgn(m+1)*a(n,m,nprime,m+1)*Zmn(alpham*abs(m+1),nprime,rho,theta);
        }
    }

    // compute the y-component of the gradient
    double gradY = 0.0;
    for(int nprime = 0; nprime <= n-1; nprime++)
    {
        if((n >= abs(m)) && (nprime >= abs(m-1)) && ((n-abs(m))%2 == 0) && ((nprime-abs(m-1))%2 == 0))
        {
            gradY += -alpham*sgn(m-1)*a(n,m,nprime,m-1)*Zmn(-alpham*abs(m-1),nprime,rho,theta);
        }
        if((n >= abs(m)) && (nprime >= abs(m+1)) && ((n-abs(m))%2 == 0) && ((nprime-abs(m+1))%2 == 0))
        {
            gradY += a(n,m,nprime,m+1)*Zmn(-alpham*abs(m+1),nprime,rho,theta);
        }
    }

    // compute the z-component of the gradient
    double gradZ = 0.0;

    return {gradX, gradY, gradZ};
}

double sgn(double val)
{
    if(val > 0)
    {
        return 1.0;
    }
    else if(val < 0)
    {
        return -1.0;
    }
    else
    {
        return 0.0;
    }
}
