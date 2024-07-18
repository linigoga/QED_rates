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
const float conversion_factor = 1.0e-6 * m_e * pow(c,2) / e_charge;  // Conversion factor from MeV to SI units
const double m =  m_e*pow(c,2) / e_charge / 1e9;

QED::QED()
{
}

    double QED::thomasFermiWavelength(int atomic_number)
    {
        /*
        Thomas-Fermi wavelength
        */
        double lambda_c = h_bar / (m_e * c);
        double lambda_tf = 4.0 * PI * epsilon_0 * pow(h_bar,2) * pow(atomic_number,-1.0/3.0)/(m_e*pow(e_charge,2));
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

   
    double QED::bremmstrahlungCrossSection(double atomic_number, double gamma) 
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
            the Bremmstrahlung cross section in cm^2
        
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
        double lower_bound = 1.0;
        double upper_bound = gamma - 1.0;
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

        if (k >= 2.0 &&  gamma_p >= 1.0 && gamma_p <= k - 1.0)
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
            double part2 = (pow(gamma_p,2) + pow(gamma_e,2)) *  (i1Function(delta, lambda_tf) + 1.0 - coulomb_correction);
            double part3 = 2.0/3.0 * gamma_e * gamma_p * (i2Function(delta, lambda_tf) + 5.0 / 6.0 - coulomb_correction);

            double diff_sigma_bethe_heitler = part1 * (part2 + part3);
            return diff_sigma_bethe_heitler;
        }
        else
        {
            return 0.0;
        }
    }

    double QED::betheHeitlerCrossSection(double k, double atomic_number)
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
            double lower_bound = 1.0;
            double upper_bound = k - 1.0;
            

            double result = boost::math::quadrature::trapezoidal(integrand, lower_bound, upper_bound, error_estimate, L1_norm);
            
            return result;
        }
        else
        {
            return 0.0;
        }
    }    

    std::tuple<double,double,double,double,double> QED::calculateCCoefficient(double gamma)
    {
        double x = 1.0 / gamma;

        double c1 = 1.0;
        double c2 = 1.0;

        double c_ = 4.0 * pow(x,2) / (1.0 - pow(x,2)) * log(1.0 / pow(x,2)) - 4.0/3.0*pow(x,2) + 1.0/6.0 * pow(x,4);
        
        double c_z = 3.0 * pow(x,2) / (1.0 - pow(x,2)) * ( 1.0 - pow(x,2)/(1.0 - pow(x,2)) * log(1.0/pow(x,2)) ) \
            - 13.0/5.0 * pow(x,2) + 7.0/4.0 * pow(x,4) - 9.0/10.0 * pow(x,6) + 1.0/5.0 * pow(x,8);
        
        double c_r = -3.0 / 2.0 * pow(x,2) / (1.0 - pow(x,2)) * ( 1.0 - pow(x,2)/(1.0 - pow(x,2)) * log(1.0/pow(x,2)) ) \
            + 4.0/5.0 * pow(x,2) - 1.0/8.0 * pow(x,4) - 1.0/20.0 * pow(x,6) + 1.0/40.0 * pow(x,8);

        return std::make_tuple(c1,c2,c_,c_z,c_r);
    }


    double QED::nonRelativisticCoulombTridentDiffCrossSection(double atomic_number, double gamma, double positron_energy)
    {
        /*
        Calculate the non-relativistic Coulomb trident differential cross section.


        Parameters
        ----------
        Input:
            atomic_number : int
            The atomic number of the target material.

            gamma : float
            The energy of the electron in units of m_e * c**2.

            positron_energy : float
            The kinetic energy of the positron in units of m_e * c**2.

        Returns
        -------
        the non-relativistic Coulomb trident differential cross section


        */
        double e_p = (positron_energy - 2.0) * m_e * pow(c,2);

        std::tuple<double,double,double,double,double> c_coefficients = calculateCCoefficient(gamma);        

        double c_ = std::get<2>(c_coefficients); // get functions obtains the element of the matrix
        double c_z = std::get<3>(c_coefficients);
        double c_r = std::get<4>(c_coefficients);


        double part1 = pow(atomic_number * alpha * r_e,2) / 32.0;
        double part2 = (log(pow(gamma,2)) - 161.0/60.0 + c_ + c_r + c_z)*pow(e_p,3) / pow(m_e*c,3);

        return part1 * part2;
    }

    double QED::relativisticCoulombTridentDiffCrossSection(double atomic_number, double gamma, double positron_energy)
    {
        /*
        Calculate the relativistic Coulomb trident differential cross section.


        Parameters
        ----------

        Input:
            atomic_number : int
            The atomic number of the target material.

            gamma : float
            The energy of the electron in units of m_e * c**2.

            positron_energy : float
            The kinetic energy of the positron in units of m_e * c**2.

        Returns
        -------
        the relativistic Coulomb trident differential cross section
        */

        std::tuple<double,double,double,double,double> c_coefficients = calculateCCoefficient(gamma);        

        double c1 = std::get<0>(c_coefficients); // get functions obtains the element of the matrix
        double c2 = std::get<1>(c_coefficients);

        double e_p = (positron_energy - 2.0) * m_e * pow(c,2); //Energy of the positron

        double part1 = 56.0 / (9.0 * PI) * pow(atomic_number * r_e * alpha,2);
        double part2 = log(c1 * e_p/( m_e*pow(c,2)));
        double part3 = log(c2 * m_e*pow(c,2) * gamma / e_p )* m_e*pow(c,2) / (e_p);

        return part1 * part2 * part3;
    }

    double QED::totalCoulombTridentCrossSectionLowerLimit(double atomic_number, double gamma)
    {
        /*
        Calculate the lower limit of the Coulomb trident cross section.


        Parameters
        ----------
        Input:
            atomic_number : int
            The atomic number of the target material.

            gamma : float
            The energy of the electron in units of m_e * c**2.

        Returns
        -------
        the lower limit of the Coulomb trident cross section
        */

       return (7.0 / 2304.0) * (atomic_number * r_e * pow(alpha,2)) * pow(gamma - 3.0,3);
    }


    double QED::totalCoulombTridentCrossSectionUpperLimit(double atomic_number, double gamma)
    {
        /*
        Calculate the upper limit of the Coulomb trident cross section.


        Parameters
        ----------
        Input:
            atomic_number : int
            The atomic number of the target material.

            gamma : float
            The energy of the electron in units of m_e * c**2.

        Returns
        -------
        the upper limit of the Coulomb trident cross section
        */

       return 28.0 / (27.0 * PI) * pow( atomic_number * r_e * alpha,2) * pow(log(gamma),3);
    }

    double QED::totalCoulombTridentCrossSection(double atomic_number, double gamma)
    {
        /*
        Calculate the Coulomb trident cross section.

        
        Parameters
        ----------
        Input:
            atomic_number : int
            The atomic number of the target material.

            gamma : float
            The energy of the electron in units of m_e * c**2.

        Returns
        -------
        the Coulomb trident cross section
        */

       if (gamma > 3.0 && gamma < 200)
        {       
            double energy_si = (gamma - 1.0) * conversion_factor;
            return  1.0e-34 * 5.22 * pow(atomic_number,2) * pow(log((2.3 + energy_si) / 3.52) , 3);
        }
        else
        {
            return totalCoulombTridentCrossSectionUpperLimit(atomic_number, gamma);
        }        
    }

    double QED::differentialCoulombTridentCrossSectionLowerLimit(double atomic_number, double gamma, double positron_energy)
    {
        /*
        Calculate the lower limit of the differential Coulomb trident cross section.


        Parameters
        ----------
        Input:
            atomic_number : int
            The atomic number of the target material.

            gamma : float
            The energy of the electron in units of m_e * c**2.

            positron_energy : float
            The kinetic energy of the positron in units of m_e * c**2.

        Returns
        -------
        the lower limit of the differential Coulomb trident cross section
        */


        double total_energy = (positron_energy - 2.0) * m_e * pow(c,2);

        std::tuple<double,double,double,double,double> c_coefficients = calculateCCoefficient(gamma);

        double c_ = std::get<2>(c_coefficients); // get functions obtains the element of the matrix
        double c_z = std::get<3>(c_coefficients);
        double c_r = std::get<4>(c_coefficients);

        double diff_cross_section = pow( atomic_number * r_e * alpha,2) / 32.0 \
            * (log(pow(gamma,2)) - 161.0/60.0 + c_ + c_r + c_z) * pow(total_energy,3) \
            / pow(m_e * pow(c,2),4) * m_e * pow(c,2);

        return diff_cross_section;

    }


    double QED::differentialCoulombTridentCrossSectionUpperLimit(double atomic_number, double gamma, double gamma_positron)
    {
        /*
        Calculate the upper limit of the differential Coulomb trident cross section.


        Parameters
        ----------
        Input:
            atomic_number : int
            The atomic number of the target material.

            gamma : float
            The energy of the electron in units of m_e * c**2.

            positron_gamma : float
            The Lorentz factor of the positron.

        Returns
        -------
        the upper limit of the differential Coulomb trident cross section
        */


        double energy_positron = gamma_positron * m_e * pow(c,2);

        std::tuple<double,double,double,double,double> c_coefficients = calculateCCoefficient(gamma);

        double c1 = std::get<0>(c_coefficients); // get functions obtains the element of the matrix
        double c2 = std::get<1>(c_coefficients);

        double diff_cross_section = 56.0 / (9.0 * PI) * pow(atomic_number * r_e * alpha,2) \
            * log( c1 * energy_positron / (m_e * pow(c,2))) * log(c2 * m_e * pow(c,2) * gamma / energy_positron) \
            * m_e * pow(c,2) / (energy_positron);

        return diff_cross_section;

    }


    double QED::differentialCoulombTridentCrossSection(double atomic_number, double gamma, double gamma_positron)
    {
        /*
        Calculate the differential Coulomb trident cross section.


        Parameters
        ----------
        Input:
            atomic_number : int
            The atomic number of the target material.

            gamma : float
            The energy of the electron in units of m_e * c**2.

            positron_energy : float
            The kinetic energy of the positron in units of m_e * c**2.

        Returns
        -------
        the differential Coulomb trident cross section
        */

       if( gamma_positron > 2.0 && gamma_positron < gamma - 1.0 ) 
       {
            double lower_limit_cross_section = differentialCoulombTridentCrossSectionLowerLimit(atomic_number, gamma, gamma_positron);
            double upper_limit_cross_section = differentialCoulombTridentCrossSectionUpperLimit(atomic_number, gamma, gamma_positron);

            return lower_limit_cross_section * upper_limit_cross_section / (upper_limit_cross_section + lower_limit_cross_section);
       }

         else
         {
            return 0.0;
         }
    }

    double QED::differentialBreitWheelerPairProduction(double gamma_photon, double chi_photon, double chi_positron) {
        /*
        Calculate the differential Breit-Wheeler pair production rate.

        Parameters
        ----------
        Input:
            gamma : double
            Incident photon energy normalized in mc^2 units.

            chi_photon : float
            Quantum nonlinearity parameter for the incident photon

            chi_positron : float
            quantum nonlinearity parameter for the produced positron

        Returns
        -------
        the differential Coulomb trident cross section
        */
        bool conditions = !(gamma_photon > 2.0 && chi_positron < chi_photon && chi_positron > chi_photon / gamma_photon);

        if (conditions) {
            return 0.0;
        }

        const double chi_electron = chi_photon - chi_positron;
        const double x = pow(chi_photon / (chi_electron * chi_positron), 2.0 / 3.0);

        double part1 = m_e * c * c * chi_photon * alpha / (h_bar * gamma_photon);
        double part2 = 1.0 / (PI * sqrt(3.0) * pow(chi_photon, 2));

        auto T1 = [](double s) {
            return sqrt(s) * cyl_bessel_k(1.0 / 3.0, (2.0 / 3.0) * pow(s, 3.0 / 2.0));
        };

        auto T2 = [x, chi_photon](double s) {
            return (2.0 - chi_photon * pow(x, 3.0 / 2.0)) * cyl_bessel_k(2.0 / 3.0, (2.0 / 3.0) * pow(x, 3.0 / 2.0));
        };

        double lower_limit = x;
        double upper_limit = 1e6; // or a large number if infinity is not supported

        // Integrate T1 from x to a large upper limit
        double result_T1 = boost::math::quadrature::trapezoidal(T1, lower_limit, upper_limit);

        double result_T2 = T2(x);

        bool test = (chi_electron > chi_photon / gamma_photon) && (chi_electron < chi_photon * (1.0 - 1.0 / gamma_photon));

        return test * part1 * part2 * (result_T1 - result_T2);
    }



QEDBlackburn::QEDBlackburn()
{
}

    double QEDBlackburn::auxilaryFunctionT(double chi)
    {
        /*
        Calculate the auxilary function T.


        Parameters
        ----------
        Input:
            chi : float
            The quantum nonlinearity parameter chi.  

        Returns
        -------
        the auxilary function T
        */

       double T = 0.16 * pow(boost::math::cyl_bessel_k(1.0/3.0, 4.0 /(3.0 * chi) ),2)/ chi;
       return T;
    }

    double QEDBlackburn::auxilaryFunctionR(double x)
    {
        /*
        Calculate the auxilary function R. It contains a modified bessel funcition
        of the second kind.


        Parameters
        ----------
        Input:
            x : float
            The quantum nonlinearity parameter chi.

        Returns
        -------
        aux_r the auxilary function R
        */

        double numerator = 0.453 * pow(boost::math::cyl_bessel_k(1.0/3.0, 4.0 /(3.0 * x)),2);
        double denominator = 1.0 +  0.145 * pow(x,1.0/4.0) * log(1.0 + 2.26*x) + 0.330*x;
        double aux_r = numerator / denominator;

        return aux_r;
    }


    double QEDBlackburn::auxilaryFunctionG(double x)
    {
        /*
        Calculate the auxilary function S. It contains a modified bessel funcition
        of the second kind.


        Parameters
        ----------
        Input:
            x : float
            The quantum nonlinearity parameter chi.

        Returns
        -------
        aux_g the auxilary function G
        */

       double  g_func = pow(1.0 + 4.8 * (1.0 + x) * log(1.0 + 1.7*x) + 2.44*pow(x,2) , -2.0/3.0);
       return g_func;
    }

    double QEDBlackburn::calculateChi(double gamma_e, double a0, double omega_0, double phi, double n)
    {
        /*
        Calculate the quantum nonlinearity parameter chi for an electron, equation 7 from
        Blackburn.


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            The laser frequency.

            phi : float
            The phase of the laser.

            n : float
            The refractive index of the medium.

        Returns
        -------
        chi the quantum nonlinearity parameter
        */

        double part1 = 2.0 * gamma_e * a0 * omega_0 * m;
        double part2 = exp(log(2) * pow(phi,2) / (2.0 / pow(PI*n,2) ));

        double chi = part1 * part2;
        return chi;
    }

    double QEDBlackburn::calculatePairProductionProbability(double a0, double gamma_e, double omega_0, double omega, double n)
    {
        /*
        Calculate the pair production probability for an electron, equation 7 from
        Blackburn.


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            The laser frequency.

            phi : float
            The phase of the laser.

            n : float
            The refractive index of the medium.


        Returns
        -------
        pair_production_probability the pair production probability
        */

        double x = 2.0 * a0 * omega_0 * omega / pow(m,2);
        double aux_r = auxilaryFunctionR(x);

        double pair_production_probability = alpha * a0 * n * aux_r;

        return pair_production_probability;
    }


    double QEDBlackburn::calculateRadiatedEnergy(double gamma_e, double a0, double omega_0, double chi_c)
    {
        /*
        Calculate the radiated energy for an electron.


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            The laser frequency.

            chi_c : float
            The quantum nonlinearity parameter chi.

        Returns
        -------
        radiated_energy the radiated energy
        */
        double part1 = sqrt(2.0 * PI) * gamma_e * m;
        double numerator = 2.0 * std::log(2.0 * gamma_e * a0 * omega_0 /(m * chi_c));
        double denominator = 1.0 + 2.0 * std::log(2.0 * gamma_e * a0 * omega_0 /(m * chi_c));

                            

        double radiated_energy = part1 * sqrt(numerator / denominator);

        return radiated_energy;
    }

    double QEDBlackburn::calculateCriticalChiRR(double gamma_e,double a0, double omega_0, double chi)
    {
        /*
        Calculate the critical quantum nonlinearity parameter chi for radiation reaction.


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            The laser frequency.

            chi : float
            The quantum nonlinearity parameter chi.

        Returns
        -------
        chi_c the critical quantum nonlinearity parameter chi as a double
        */

        double radiated_energy = calculateRadiatedEnergy(gamma_e, a0, omega_0, chi);
        double chi_critica_rr = chi / (1.0 + radiated_energy / (2.0 * gamma_e * m));

        return chi_critica_rr;
    }




    double QEDBlackburn::criticalChi(double gamma_e, double a0, double omega_0, double n)
        {   
            /*
            Calculate the critical quantum nonlinearity parameter chi.


            Parameters
            ----------
            Input:
                gamma_e : float
                The Lorentz factor of the electron.

                a0 : float
                The normalized vector potential.

                omega_0 : float
                The laser frequency.

                n : float
                The refractive index of the medium.

            Returns
            -------
            chi_c the critical quantum nonlinearity parameter chi as a double
            */

        // Precompute constants
        const double CONSTANT_PART = 72.0 * log(2.0) / pow(PI * alpha, 2) * pow(gamma_e * omega_0 / (n * m), 2);

        auto chi_func = [this, &a0, &gamma_e, &omega_0, &CONSTANT_PART](double chi) -> double {
            double g_chi = this->auxilaryFunctionG(chi);
            return pow(chi, 4) * g_chi * g_chi - CONSTANT_PART * log(2.0 * (gamma_e * a0 * omega_0) / (m * chi));
        };

        double min_guess = 1e-10;
        double max_guess = 1e-5;
        const int MAX_ADJUSTMENTS = 100;
        int adjustments = 0;

        // Ensure that chi_func(min_guess) and chi_func(max_guess) have opposite signs
        while (chi_func(min_guess) * chi_func(max_guess) > 0 && adjustments < MAX_ADJUSTMENTS) {
            min_guess /= 10.0;
            max_guess *= 10.0;
            adjustments++;
        }

        if (adjustments == MAX_ADJUSTMENTS) {
            throw std::runtime_error("Failed to bracket the root after maximum adjustments.");
        }

        auto termination_condition = [](double min, double max) {
            return abs(min - max) <= boost::math::tools::eps_tolerance<double>()(std::max(abs(min), abs(max)), 32);
        };

        std::pair<double, double> result = boost::math::tools::bisect(chi_func, min_guess, max_guess, termination_condition);

        return (result.first + result.second) / 2.0;
    }




    double QEDBlackburn::criticalChiModulus(double gamma_e, double a0, double omega_0, double n)
{
    // Precompute constants
    const double CONSTANT_PART = 72.0 * log(2.0) / pow(PI * alpha, 2) * pow(gamma_e * omega_0 / (n * m), 2);

    auto chi_func = [&a0, &gamma_e, &omega_0, &CONSTANT_PART](double chi) -> double {
        return pow(chi, 4) - CONSTANT_PART * log(2.0 * (gamma_e * a0 * omega_0) / (m * chi));
    };

    double min_guess = 1e-10;
    double max_guess = 1e2;

    // Ensure that chi_func(min_guess) and chi_func(max_guess) have opposite signs
    while (chi_func(min_guess) * chi_func(max_guess) > 0) {
        if (chi_func(min_guess) > 0) {
            min_guess *= 10.0;
        } else {
            max_guess /= 10.0;
        }

        // Check if we've exhausted our search interval
        if (min_guess >= max_guess) {
            throw std::runtime_error("The provided interval does not bracket a root.");
        }
    }

    auto termination_condition = [](double min, double max) {
        return abs(min - max) <= boost::math::tools::eps_tolerance<double>()(std::max(abs(min), abs(max)), 32);
    };

    std::pair<double, double> result = boost::math::tools::bisect(chi_func, min_guess, max_guess, termination_condition);

    return (result.first + result.second) / 2.0;
}



    double QEDBlackburn::criticalFrequency(double gamma_e, double a0, double omega_0, double chi_c)
    {
        /*
        Calculate the critical frequency.


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            reference frequency.

            chi_c : float
            The critical quantum nonlinearity parameter chi.

        Returns
        -------
            omega_c the critical frequency as a double
        */

        double critical_chi_rr = calculateCriticalChiRR(gamma_e, a0, omega_0, chi_c);
        double sqrt_argument = 2.0*critical_chi_rr * m / (a0 * gamma_e * omega_0);

        double numerator = gamma_e * m * sqrt(sqrt_argument);
        double denominator = 1.0 + sqrt(sqrt_argument);

        double omega_c = numerator / denominator;
        return omega_c;

    }

    double QEDBlackburn::calculateCriticalPhase(double gamma_e, double a0, double omega_0, double n)
    {
        /*
        Calculate the critical phase.


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            reference frequency.

            n : float
            The refractive index of the medium.


        Returns
        -------
            phi_c the critical phase as a double
        */




        double chi_c = criticalChi(gamma_e, a0, omega_0, n);
        long double part1 = 2.0 * pow(PI*n,2) / log(2.0);
        long double part2 = log(2.0 * gamma_e * a0 * omega_0 / (chi_c * m));
                        
        double phi_critical = sqrt(part1 * part2);

        return phi_critical;
    }


     double QEDBlackburn::calculateCriticalPhaseModulus(double gamma_e, double a0, double omega_0, double n)
    {
        /*
        Calculate the critical phase insie the modulus.


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            reference frequency.

            n : float
            number of cycles.


        Returns
        -------
            phi_c the critical phase as a double
        */


        double chi_c = criticalChiModulus(gamma_e, a0, omega_0, n);
        double part1 = 2.0 * pow(PI*n,2) / std::log(2.0)*std::log((2.0 * gamma_e * a0 * omega_0 )/ (chi_c * m));
        
        double phi_critical = std::sqrt(part1);

    

        return phi_critical;
    }




    double QEDBlackburn::correctionFactorFhe(double n, double phi_c)
    {
        /*
        Calculate the correction factor Fhe.


        Parameters
        ----------
        Input:
            n : float
            The refractive index of the medium.

            phi_c : float
            The critical phase.

        Returns
        -------
            f_he the correction factor Fhe as a double
        */

       double erf_function_argument = sqrt(2.0* log(2)) * phi_c/(2 * PI * n);

       double f_he = 0.5*(1.0 - erf(erf_function_argument));



       return f_he;
    }

double QEDBlackburn::calculatePhotonEnergySpectrum(double gamma_e, double a0, double omega_0, double omega, double n)
{
    /*
    Calculate the photon energy spectrum.


    Parameters
    ----------
    Input:
        gamma_e : float
        The Lorentz factor of the electron.

        a0 : float
        The normalized vector potential.

        omega_0 : float
        The laser frequency.

        omega : float
        The photon frequency.

        n : float
        The refractive index of the medium.

    Returns
    -------
        photon_energy_spectrum the photon energy spectrum as a double
    */

    double phi_c = criticalFrequency(gamma_e, a0, omega_0, n);
    double chi_c = criticalChi(gamma_e, a0, omega_0, n);
    double corr_factor = correctionFactorFhe(n, phi_c);
    double critical_chi_rr = calculateCriticalChiRR(gamma_e, a0, omega_0, chi_c);

    double energy_0 = gamma_e * m;
    
    double chi_0 = 2.0 * gamma_e * a0 * omega_0 / m;

    double part1 = sqrt(3.0) * PI * alpha * corr_factor / sqrt( 2.0 * std::log(2.0));
                   
    double part2 = a0 * n * critical_chi_rr / chi_0 / (sqrt(energy_0)* sqrt(1.0 + 2.0 * log(chi_0/chi_c)));
    double part3 = exp(- 2.0 * omega / (3.0 * critical_chi_rr * (energy_0 - omega))) / sqrt( 3.0 * critical_chi_rr * (energy_0 - omega) + 4.0 * omega);


    double dn_domega = part1 * part2 * part3;

    return dn_domega;
}

    double QEDBlackburn::positronYield(double gamma_e, double a0, double omega_0,double n)
    {
        /*
        Calculate the positron yield.


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            The laser frequency.

            n : float
            The refractive index of the medium.

        Returns
        -------
            positron_yield the positron yield as a double
        */

        
        double critical_chi = criticalChi(gamma_e, a0, omega_0, n);
        double omega_critical = criticalFrequency(gamma_e, a0, omega_0, n);
        double pair_prod_probability = calculatePairProductionProbability(a0, gamma_e, omega_0, omega_critical, n);
        double chi_critical_rr = calculateCriticalChiRR(gamma_e, a0, omega_0, critical_chi);
        double dn_domega = calculatePhotonEnergySpectrum(gamma_e, a0, omega_0, omega_critical, n);


        double part1 = 3.0 * PI * pair_prod_probability * chi_critical_rr / sqrt(2);
        double part2 = pow(gamma_e * m - omega_critical,2) / (gamma_e * m) * dn_domega;

        double positron_yield = part1 * part2;

        return positron_yield;
    }
    double QEDBlackburn::calculateVariance(double gamma_e, double a0, double omega_0, double n, double chi_c)
    {
        /*
        Calculate the variance of the quantum nonlinearity parameter chi.


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            The laser frequency.

            n : float
            The refractive index of the medium.

            chi : float
            The quantum nonlinearity parameter chi. 

        Returns
        -------
            variance the variance of the quantum nonlinearity parameter chi as a double
        */

       double variance = pow(PI*n,2)/std::log(2.0) / ( 1.0 + 2.0 * std::log(2.0 * gamma_e * a0 * omega_0 / (m * chi_c)));

       return variance;
    }

    double QEDBlackburn::calculateChiAsFunctionOfPhi(double gamma_e, double a0, double omega_0, double phi, double n)
    {
        /*
        Calculate the quantum nonlinearity parameter chi as a function of the phase phi.
        equation 17 of the manuscript


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            The laser frequency.

            phi : float
            The phase of the laser.

            n : float
            number of cycles.

        Returns
        -------
            chi the quantum nonlinearity parameter as a double
        */

        double critical_chi = criticalChi(gamma_e, a0, omega_0, n);
        double critical_phi = calculateCriticalPhase(gamma_e, a0, omega_0, n);
        double radiated_energy = calculateRadiatedEnergy(gamma_e, a0, omega_0, critical_chi);

        double variance = calculateVariance(gamma_e, a0, omega_0, n,critical_chi);

        double chi = critical_chi / (1.0 + radiated_energy / (2.0 * gamma_e * m)) * exp(-pow(phi - critical_phi,2) / (2.0 * variance));
                     
        return chi;
    }

    double QEDBlackburn::calculateGammaAsFunctionOfPhi(double gamma_e, double a0, double omega_0, double phi, double n)
    {
        /*
        Calculate the Lorentz factor gamma as a function of the phase phi.
        equation 16 of the manuscript


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            The laser frequency.

            phi : float
            The phase of the laser.

            n : float
            number of cycles.

        Returns
        -------
            gamma the Lorentz factor as a double
        */

        double critical_chi = criticalChi(gamma_e, a0, omega_0, n);
        double critical_phi = calculateCriticalPhase(gamma_e, a0, omega_0, n);
        
        double radiated_energy = calculateRadiatedEnergy(gamma_e, a0, omega_0, critical_chi); 
        
        double standard_deviation = pow(PI*n,2)/std::log(2.0) * 1.0 / ( 1.0 + 2.0 * std::log(2.0 * gamma_e * a0 * omega_0 / (m * critical_chi)));

        double gamma_f = (2.0 * gamma_e * m - radiated_energy) / (2.0 * gamma_e * m + radiated_energy) * gamma_e;

        double gamma = gamma_f + gamma_e * radiated_energy / (2.0 * gamma_e * m + radiated_energy) * ( 1.0 + erf( phi - critical_phi) / 2.0 * sqrt(standard_deviation));

        return gamma;

    }

QEDReconstructionMethods::QEDReconstructionMethods()
{
}


    double QEDReconstructionMethods::particleBinningAmaro(std::map<std::string,float> cache, std::string beam_type)
    {
        /*
        Calculate the particle binning for the Amaro method.


        Parameters
        ----------
        Input:
            cache : map
            The cache of the simulation.

            beam_type : string
            The type of beam.

        Returns
        -------
            particle_binning : double
                The differential rate of photon production by nonlinear inverse
                compton scattering.
        */

        if (beam_type == "widebeam")
            {
            double n0 = cache["n0"];
            double w0 = cache["w0"];
            double a0 = cache["a0"];
            
            double zmin = cache["zmin"];
            double zmax = cache["zmax"];
            double z_r = cache["z_r"];
            
            
            double part1 = 2 * PI* n0 * pow(w0 ,2) / a0;

            auto integrand = [&z_r](double z) -> double {
                return 1.0 + pow(z/z_r,2);
            };
            double error_estimate = 0.0;
            double L1_norm = 1.0;
            double integrated_value = boost::math::quadrature::trapezoidal(integrand, zmin,zmax, error_estimate, L1_norm);

            double dnda = part1*integrated_value;

            return dnda;
            }
        else if (beam_type == "shortbeam")
            {   

            double a0_eff = cache["a0_eff"];
            double n0 = cache["n0"];
            double W0 = cache["W0"];
            double a0 = cache["a0"];
            double R = cache["R"];
            double delta = cache["delta"];
            
            double part1 = n0 * pow(W0,2) /(pow(R,2) * a0_eff) * pow(a0_eff/a0, pow(W0/R,2));
            
            double exp1 = exp(-pow(delta,2) / pow(R,2));

            double bess_argument = 2.0 * delta * W0 / pow(R,2) * sqrt(log(a0/a0_eff));
            double part2 = boost::math::cyl_bessel_i(0,bess_argument);
            
            
            double dnda = part1 * exp1 * part2;

            return dnda;
            }
        else if (beam_type == "thinbeam")
            {

            double N = cache["N"];
            double z_r = cache["z_r"];
            double a0 = cache["a0"];
            double a0_eff = cache["a0_eff"];
            double L = cache["L"];

            double part1 = 4*N*z_r/L;
            double part2 = pow(a0,2) * pow(a0_eff,2);
            double part3 = 1/sqrt(pow(a0,2) - pow(a0_eff,2));

            return part1*part2*part3;
            }

        else{
            cout << "Beam type not recognized" << endl;
            return 0.0;
        }
    }

    double QEDReconstructionMethods::make3dGaussDistribution(double x, double y, double z, double a0, double waist, double wavelength)
    {
        /*
        Calculate the 3D Gaussian distribution.


        Parameters
        ----------
        Input:
            x : float
            The x coordinate.

            y : float
            The y coordinate.

            z : float
            The z coordinate.

            a0 : float
            The normalized vector potential.

            waist : float
            The waist of the laser.

            wavelength : float
            The wavelength of the laser.

        Returns
        -------
            gauss_distribution : double
                The 3D Gaussian distribution.
        */


        double z_r = PI * pow(waist,2) / wavelength;

        double r2 = pow(x,2) + pow(y,2);

        double gaussian_dist = a0/sqrt(1.0 + pow(z/z_r,2)) * exp(-r2/pow(waist,2) * 1.0/(1.0 + pow(z/z_r,2)));
        return gaussian_dist;
    }

    double QEDReconstructionMethods::calculatePositronsProducedFromBeam(double gamma_e, double a0, double omega_0,double delta, double r_, double waist,double ne,double n)
    {
        /*
        Calculate the positrons produced from the beam.


        Parameters
        ----------
        Input:
            gamma_e : float
            The Lorentz factor of the electron.

            a0 : float
            The normalized vector potential.

            omega_0 : float
            The laser frequency.

            delta : float
            The distance from the center of the beam.

            r_ : float
            The radius of the beam.

            waist : float
            The waist of the laser.

            ne : float
            The electron number density.

            n : float
            The refractive index of the medium.

        Returns
        -------
            positrons_produced : double
                The positrons produced from the beam.
        */

        QEDBlackburn blackburn;                 //One has to initialize the class and then call the function
        double chi_c = blackburn.criticalChi(gamma_e, a0, omega_0, n);
        double omega_c = blackburn.criticalFrequency(gamma_e, a0, omega_0, chi_c);
        double number_of_positrons = blackburn.positronYield(gamma_e, a0, omega_0, n);

        double part1 = 0.727 * a0 * omega_0 * omega_c /pow(m,2);
        double part2 = pow(waist,2) *exp( - pow(delta,2) / pow(r_,2)) / pow(r_,2);

        double positrons_from_beam = part1 * part2 * number_of_positrons * ne;
        return positrons_from_beam;
    }









