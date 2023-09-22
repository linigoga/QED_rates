#include "QEDProcesses.h"
using namespace std;
using namespace boost::math::constants;
using namespace boost::math;
using BremmstrahlungFunc = double (*)(double, double, double);

const float h_bar = 1.0545718e-34; // J s
const float c = 299792458; // m/s
const float e_charge = 1.60217662e-19; // C
const float m_e = 9.10938356e-31; // kg
const float epsilon_0 = 8.854187817e-12; // F/m
const float PI = constants::pi<double>();
const float alpha = pow(e_charge,2) / (4 * PI * epsilon_0 * h_bar * c); // fine structure constant
const float r_e = pow(e_charge,2) / (4 * PI * epsilon_0 * m_e * pow(c,2)); // classical electron radius


QED::QED()
{
}

    double QED::thomasFermiWavelength(int atomic_number)
    {
        /*
        Thomas-Fermi wavelength
        */
        double lambda_c = h_bar / (m_e * c);
        double lambda_tf = 4.0 * PI * epsilon_0 * pow(h_bar,2) * pow(atomic_number,-1.0/3.0)/(m_e*pow(c,2));
        return 0.855*lambda_tf / lambda_c;
    }

    double QED::i1Function(double delta, double screening_length, double q)
    {
        /*
        I1 function
        */        
        double part1 = screening_length * delta * (atan(delta * screening_length) - atan(screening_length));
        double part2 = - pow(screening_length,2) * pow((1 - delta),2) / (1 + pow(screening_length,2));
        double part3 = 0.5 * log( (1.0 + pow(screening_length,2)) / (1.0 + pow(delta * screening_length,2)) );

        double i1 = part1 + part2 + part3;

        return i1;
    }

    double QED::i2Function(double delta, double screening_length, double q)
    {
        /*
        I2 function
        */
        double part1 = 4.0 * pow(screening_length * delta,3) * (atan(delta * screening_length) - atan(screening_length));
        double part2 = (1.0 + 3.0 * pow(screening_length*delta,2)) * log(( 1.0 + pow(screening_length,2)) / (1.0 + pow(delta * screening_length,2)));
        double part3 = (6.0 * pow(screening_length,4) * pow(delta,2)) / (1.0 + pow(screening_length,2)) * log(delta);
        double part4 = (pow(screening_length,2) * (delta - 1.0 ) * (delta + 1.0 - 4.0 * pow(screening_length*delta,2) ) ) / (1.0 + pow(screening_length,2));

        double i2 = (part1 + part2 + part3 + part4) / 2.0;
        return i2;
    }

    double QED::nonRelativisticBremsstrahlungCrossSection(double gamma1, double k, double atomic_number)
    {
        /*
            Calculate the non-relativistic bremmstrahlung differential cross section.

    
            Parameters
            ----------
            Input:
                atomic_number : int
                    The atomic number of the target material.
                k : float
                    The energy of the photon in units of m_e * c**2.
                gamma1 : float
                    The energy of the electron in units of m_e * c**2.

            Returns
            -------
            the Bremmstrahlung differential cross section
        */

       if (k > 0.0 && k < gamma1 - 1.0)
       {
            double gamma2 = gamma1 - k;

            double p1 = sqrt(pow(gamma1,2) - 1.0);
            double p2 = sqrt(pow(gamma2,2) - 1.0);
        
            double delta_p_plus = p1 + p2;
            double delta_p_minus = p1 - p2;

            double beta1 = sqrt(1.0 - 1.0 / pow(gamma1,2));
            double beta2 = sqrt(1.0 - 1.0 / pow(gamma2,2));

            double lambda_tf = thomasFermiWavelength(atomic_number);

            double part1 = 8.0 * pow(r_e,2) * atomic_number *  alpha / (3.0 * k * pow(p1 , 2));
            double part2 = log( (pow(delta_p_plus,2) * pow(lambda_tf,2) + 1.0) / (pow(delta_p_minus,2) * pow(lambda_tf,2) + 1.0));
            double part3 = 1.0 / (pow(delta_p_plus* lambda_tf , 2)  + 1.0) - 1.0 / (pow(delta_p_minus,2) * pow(lambda_tf,2) + 1.0);

            if (atomic_number * alpha* (1.0/beta2 - 1.0/beta1) < 1.0)
            {
                double elwert = beta1/beta2 * (1.0 - exp(- 2*PI * atomic_number * alpha / beta1) )/( 1.0 - exp(- 2*PI * atomic_number * alpha / beta2));
                return part1 * (part2 + part3) * elwert;    
            }
            else
            {
                return  part1 * (part2 + part3);    
            } 
        }
        else
        {
            return 0.0;
        }
    }

    double QED::relBremmstrahlungDifferentialCrossSection(double atomic_number, double k, double gamma1)
    {
        /*      
            Calculate the relativistic bremmstrahlung differential cross section.

    
            Parameters
            ----------
            Input:
                atomic_number : int
                    The atomic number of the target material.
                k : float
                    The energy of the photon in units of m_e * c**2.
                gamma1 : float
                    The energy of the electron in units of m_e * c**2.

            Returns
            -------
            the Bremmstrahlung differential cross section
        */

        double lambda_tf = thomasFermiWavelength(atomic_number);
    
        double delta = k / (2.0 * gamma1 * (gamma1 - k));

        double part1 = 4.0 * pow(atomic_number,2) * pow(r_e,2.0) * alpha / k * (1.0 + pow((gamma1 - k) / gamma1,2.0) * i1Function(delta, lambda_tf));
        double part2 = - 2.0/3.0 * (gamma1 - k)/gamma1 *( i2Function(delta, lambda_tf) + 5.0 / 6.0);

        double diff_sigma_brem = part1 + part2;
        return diff_sigma_brem;
    }

    double QED::coulombCorrectionTerm(double atomic_number, int n)
    {
        double part1 = pow(atomic_number * alpha,2) / (1.0 + pow(atomic_number * alpha,2));
        double sum = 0.0;
        for (int i = 1; i < n; i++)
        {
            sum += pow(alpha * atomic_number,i) * (boost::math::zeta(2*i + 1) - 1.0);
        }
        return part1 * sum;
    }

    double QED::ultraRelBremmstrahlungDifferentialCrossSection(double atomic_number, double k,double gamma1)
    {
    /*
        Calculate the ultra relativistic bremmstrahlung differential cross section.

        
        Parameters
        ----------
        Input:
            atomic_number : int
                The atomic number of the target material.
            k : float
                The energy of the photon in units of m_e * c**2.
            gamma1 : float
                The energy of the electron in units of m_e * c**2.

        Returns
        -------
        the Bremmstrahlung differential cross section
    
    */
    if (k > 0.0 && k < gamma1)
        {
            double lambda_c = h_bar / (m_e * c);
            double lambda_tf = thomasFermiWavelength(atomic_number)/lambda_c;

            double delta = k / (2.0 * gamma1 * (gamma1 - k));

            double part1 = 4.0 * pow(atomic_number,2) * pow(r_e,2) * alpha / k;
            double part2 = (1.0 + pow((gamma1 - k) / gamma1,2 ) * i1Function(delta, lambda_tf));
            double part3 = -2.0/3.0 * (gamma1 - k)/gamma1 * (i2Function(delta, lambda_tf) + 5.0 / 6.0 - coulombCorrectionTerm(atomic_number));

            return part1 * (part2 + part3);
        }
        else
        {
            return 0.0;
        }
    }

   
    double QED::bremmstrahlungCrossSection(double atomic_number, double gamma, double lower_bound, double upper_bound) 
    {
        /*

        Calculate the bremmstrahlung cross section.

        Parameters
        ----------
        Input:
            atomic_number : int
                The atomic number of the target material.
            gamma : float
                The energy of the electron in units of m_e * c**2.
            lower_bound : float
                The lower bound of the integration.
            upper_bound : float
                The upper bound of the integration.

        Returns
        -------
            the Bremmstrahlung cross section
        
        */


        // Get the appropriate function based on gamma

        std::function<double(double, double, double)> bremFunc;


        // Determine the appropriate function based on gamma
        if (gamma >= 2.0)
        {
            bremFunc = [this](double atomic_number, double k, double gamma1) -> double {
                return this->nonRelativisticBremsstrahlungCrossSection(atomic_number, k, gamma1);
            };

        }
        else if (gamma > 1.0 && gamma <= 100.0)
        {
            bremFunc = [this](double atomic_number, double k, double gamma1) -> double {
                QED qed;
                return this->relBremmstrahlungDifferentialCrossSection(atomic_number, k, gamma1);
            };
        }
        else
        {
            bremFunc = [this](double atomic_number, double k, double gamma1) -> double {
                return this->ultraRelBremmstrahlungDifferentialCrossSection(atomic_number, k, gamma1);
            };
        }

        // Use Boost's trapezoidal integration method
        auto integrand = [&bremFunc, &atomic_number, &gamma](double k) {
            return bremFunc(atomic_number, k, gamma);
        };

        double error_estimate = 0.0;
        double L1_norm = 1;
        double result = boost::math::quadrature::trapezoidal(integrand, lower_bound, upper_bound, error_estimate, L1_norm);

        return result;
    }

    double QED::diffBetheHeitlerCrossSection(double gamma_p, double k, double atomic_number)
    {
        /*
        Calculate the differential Bethe-Heitler cross section. 


        Parameters
        ----------
        Input:
        gamma_p : float
            The energy of the positron in units of m_e * c**2.
        
        k : float
            The energy of the photon in units of m_e * c**2.
        
        atomic_number : int
            The atomic number of the target material.

        Returns
        -------
        the differential Bethe-Heitler cross section
        */

        if (k >= 2.0 &&  gamma_p >= 1.0 && gamma_p >= k - 1.0)
        {
            double coulomb_correction;
            double gamma_e = k - gamma_p;
            double delta = k / (2.0 * gamma_p * gamma_e);
            double lambda_tf = thomasFermiWavelength(atomic_number);

            if (k > 200.0)
            {
                coulomb_correction = coulombCorrectionTerm(atomic_number);
            }
            else
            {
                coulomb_correction = 0.0;
            }

            double part1 = 4.0 * pow(atomic_number * r_e,2) * alpha / pow(k,3);
            double part2 = pow(gamma_p*gamma_e,2) *  (i1Function(delta, lambda_tf) + 1.0 - coulomb_correction);
            double part3 = 2.0/3.0 * gamma_e * gamma_p * (i2Function(delta, lambda_tf) + 5.0 / 6.0 - coulomb_correction);

            double diff_sigma_bethe_heitler = part1 * (part2 + part3);
            return diff_sigma_bethe_heitler;
        }
        else
        {
            return 0.0;
        }
    }

    double QED::betheHeitlerCrossSection(double k, double atomic_number, double lower_bound, double upper_bound)
    {
        /*
        Calculate the Bethe-Heitler cross section.
        */
        if (k > 2.0)
        {          
            auto integrand = [this, &k, &atomic_number](double gamma_p) {
                return diffBetheHeitlerCrossSection(gamma_p, k, atomic_number);
            };

            double error_estimate = 0.0;
            double L1_norm = 1.0;
            double result = boost::math::quadrature::trapezoidal(integrand, lower_bound, upper_bound, error_estimate, L1_norm);
            
            return result;
        }
        else
        {
            return 0.0;
        }
    }    
