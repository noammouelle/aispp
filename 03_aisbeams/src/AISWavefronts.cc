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

    double phase = abs(kz) * r * r * RInv / 2.0 - atan(z/zR);

    return phase; 
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


std::complex<double> gaussianCarrierField(const std::array<double, 3>& shiftedPos,
                                          const double& kz, const double& w0)
{
    // gaussianEnvelope's own formula, inlined to avoid depending on the pulse
    // time arguments it takes but does not use for the spatial profile.
    double r  = std::sqrt(shiftedPos[0]*shiftedPos[0] + shiftedPos[1]*shiftedPos[1]);
    double z  = shiftedPos[2];
    double zR = pi * w0 * w0 / lambdaSr87;
    double w  = w0 * std::sqrt(1 + (z/zR)*(z/zR));
    double amplitude = std::exp(-r*r/(w*w)) * w0 / w;

    double phase = gaussianWavefront(shiftedPos, kz, w0);
    return std::polar(amplitude, phase);
}

std::complex<double> fourierModeField(const std::array<double, 3>& shiftedPos,
                                      const double& transverseX,
                                      const double& transverseY,
                                      const double& dzFromMirror,
                                      const double& kz, const double& w0,
                                      const FourierMode& mode)
{
    std::complex<double> carrier = gaussianCarrierField(shiftedPos, kz, w0);
    if(!mode.active())
    {
        return carrier;
    }

    const double k = std::abs(kz);
    const double q = mode.q();

    // Each sideband walks off by q dz / k, along the mode's own direction.
    const double walk = q * dzFromMirror / k;
    const double ux = (q > 0.0) ? mode.qx / q : 0.0;
    const double uy = (q > 0.0) ? mode.qy / q : 0.0;

    std::array<double,3> plusPos  = {shiftedPos[0] - walk*ux, shiftedPos[1] - walk*uy, shiftedPos[2]};
    std::array<double,3> minusPos = {shiftedPos[0] + walk*ux, shiftedPos[1] + walk*uy, shiftedPos[2]};

    const double carrierPhase = mode.qx*transverseX + mode.qy*transverseY + mode.psi;
    const std::complex<double> talbot =
        std::polar(1.0, -q*q*dzFromMirror/(2.0*k));

    std::complex<double> sidebands =
          std::polar(1.0,  carrierPhase) * gaussianCarrierField(plusPos,  kz, w0)
        + std::polar(1.0, -carrierPhase) * gaussianCarrierField(minusPos, kz, w0);

    return carrier + std::complex<double>(0.0, 0.5*mode.epsilon) * talbot * sidebands;
}


void fourierModeFirstOrder(const std::array<double, 3>& shiftedPos,
                           const double& transverseX,
                           const double& transverseY,
                           const double& dzFromMirror,
                           const double& kz, const double& w0,
                           const FourierMode& mode,
                           bool conjugateCarrier,
                           double& phaseOut,
                           double& amplitudeFactorOut)
{
    phaseOut = 0.0;
    amplitudeFactorOut = 1.0;
    if(!mode.active())
    {
        return;
    }

    const double k = std::abs(kz);
    const double q = mode.q();
    const double walk = q * dzFromMirror / k;
    const double ux = (q > 0.0) ? mode.qx / q : 0.0;
    const double uy = (q > 0.0) ? mode.qy / q : 0.0;

    std::array<double,3> plusPos  = {shiftedPos[0] - walk*ux, shiftedPos[1] - walk*uy, shiftedPos[2]};
    std::array<double,3> minusPos = {shiftedPos[0] + walk*ux, shiftedPos[1] + walk*uy, shiftedPos[2]};

    std::complex<double> carrier = gaussianCarrierField(shiftedPos, kz, w0);
    std::complex<double> plus    = gaussianCarrierField(plusPos,  kz, w0);
    std::complex<double> minus   = gaussianCarrierField(minusPos, kz, w0);
    if(conjugateCarrier)
    {
        carrier = std::conj(carrier);
        plus    = std::conj(plus);
        minus   = std::conj(minus);
    }
    if(std::abs(carrier) <= 0.0)
    {
        return;
    }

    const double carrierPhase = mode.qx*transverseX + mode.qy*transverseY + mode.psi;
    const std::complex<double> talbot = std::polar(1.0, -q*q*dzFromMirror/(2.0*k));
    const std::complex<double> ratio =
        std::complex<double>(0.0, 0.5*mode.epsilon) * talbot
        * (std::polar(1.0, carrierPhase) * plus + std::polar(1.0, -carrierPhase) * minus)
        / carrier;

    // E = E0 (1 + ratio)  =>  to first order  phase = Im(ratio), |E|/|E0| = 1 + Re(ratio).
    phaseOut = ratio.imag();
    amplitudeFactorOut = 1.0 + ratio.real();
}
