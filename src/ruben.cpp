#include <R.h>
#include "Rmath.h"

using std::exp;
using std::log;
using std::fabs;
using std::sqrt;

extern "C" {

    void ruben(double *lambda, int *mult, double *delta, int *n, double *c, double *mode, int *maxit, double *eps, double *dnsty, int *ifault, double *res) {
      
        // Algorithm AS 204 Appl. Statist. (1984) Vol. 33, No.3
        // ruben evaluates the probability that a positive definite quadratic form in Normal variates is less than a given value
      
        if ((n[0] < 1) || (c[0] <= 0) || (maxit[0] < 1) || (eps[0] <= 0.0)) {
            res[0] = -2.0;
            ifault[0] = 2;
            return;
        } else {
            int i,k,m,j;
            double pnorm(double q, double mean, double sd, int lower_tail, int log_p);
            double ao, aoinv, z, bbeta, eps2, hold, hold2, sum, sum1, dans, lans, pans, prbty, tol;
            double *gamma, *theta, *a, *b;
            gamma = new double[n[0]];
            theta = new double[n[0]];
            a = new double[maxit[0]];
            b = new double[maxit[0]];
            tol = -200.0;

            // preliminaries
            
            sum = lambda[0];
            bbeta = sum;
            for (i = 1; i <= n[0]; i++) {
                hold = lambda[i - 1];
              	if ((hold <= 0.0) || (mult[i - 1] < 1) || (delta[i - 1] < 0.0)) {
                    res[0] = -7.0;
                    ifault[0] = -i;
                    delete[] gamma;
                    delete[] theta;
                    delete[] a;
                    delete[] b;
                    return;
                }
              	if (bbeta > hold) bbeta = hold;
              	if (sum < hold) sum = hold;
            }
            if (mode[0] > 0.0) {
                bbeta = mode[0] * bbeta;                   // AS204A
            } else {
                bbeta = 2.0 / (1.0 / bbeta + 1.0 / sum);   // AS204B
            }
            k = 0;
            sum = 1.0;
            sum1 = 0.0;
            for (i = 1; i <= n[0]; i++) {
                hold = bbeta / lambda[i - 1];
                gamma[i - 1] = 1.0 - hold;
                sum = sum * R_pow(hold, mult[i - 1]); 
                sum1 = sum1 + delta[i - 1];
                k = k + mult[i - 1];
                theta[i - 1] = 1.0;
            }
            ao = std::exp(0.5 * (std::log(sum) - sum1));
            if (ao <= 0.0) {
                res[0] = 0.0;
                dnsty[0] = 0.0;
                ifault[0] = 1;
            } else { 
                z = c[0] / bbeta;

                // evaluate probability and density of chi-squared 
                // on k degrees of freedom. The constant 0.22579135264473 
                // is ln(sqrt(pi / 2))

                if ((k % 2) == 0) {
                    i = 2;
                    lans = -0.5 * z;
                    dans = std::exp(lans);
                    pans = 1.0 - dans;
                } else {
                    i = 1;
                    lans = -0.5 * (z + std::log(z)) - 0.22579135264473;
                    dans = std::exp(lans);
                    pans = pnorm(std::sqrt(z), 0.0, 1.0, 1, 0) - pnorm(-std::sqrt(z), 0.0, 1.0, 1,0); 
                }
                k = k - 2;
                for (j = i; j <= k; j = j + 2) {
                    if (lans < tol) {
                        lans = lans + std::log(z / (double)j);
                        dans = std::exp(lans);
                    } else {
                        dans = dans * z / (double)j;
                    }
                    pans = pans - dans;
                }
        
                // evaluate successive terms of expansion
        
                prbty = pans;
                dnsty[0] = dans;
                eps2 = eps[0] / ao;
                aoinv = 1.0 / ao;
                sum = aoinv - 1.0;
                for (m = 1; m <= maxit[0]; m++) {
                    sum1 = 0.0;
                    for (i = 1; i <= n[0]; i++) {
                        hold = theta[i - 1];
                        hold2 = hold * gamma[i - 1];
                        theta[i - 1] = hold2;
                        sum1 = sum1 + hold2 * mult[i - 1] + m * delta[i - 1] * (hold - hold2);
                    }
                    sum1 = 0.5 * sum1;
                    b[m - 1] = sum1;
                    for (i = m - 1; i >= 1; i--) {
                        sum1 = sum1 + b[i - 1] * a[m - i - 1]; 
                    }
                    sum1 = sum1 / (double)m;
                    a[m - 1] = sum1;
                    k = k + 2;
                    if (lans < tol) {
                      	lans = lans + std::log(z / (double)k);
                      	dans = std::exp(lans);
                    } else {
                        dans = dans * z / (double)k;
                    }
                    pans = pans - dans;
                    sum = sum - sum1;
                    dnsty[0] = dnsty[0] + dans * sum1;
                    sum1 = pans * sum1;
                    prbty = prbty + sum1;
                    if (prbty < (-aoinv)) {
                        res[0] = -3.0;
                        ifault[0] = 3;
                        delete[] gamma;
                        delete[] theta;
                        delete[] a;
                        delete[] b;
                        return;
                    }
                    if (std::fabs(pans * sum) < eps2) {
                        if (std::fabs(sum1) < eps2) {
                            ifault[0] = 0;
                            break;
                        }
                    }
                }
                if (m >= maxit[0]) ifault[0] = 4;
                dnsty[0] = ao * dnsty[0] / (bbeta + bbeta);
                prbty = ao * prbty;
                if (prbty < 0.0 || prbty > 1.0) {
                    ifault[0] = ifault[0] + 5;
                } else {
                    if (dnsty[0] < 0.0) ifault[0] = ifault[0] + 6;
                }
                res[0] = prbty;
            }
            delete[] gamma;
            delete[] theta;
            delete[] a;
            delete[] b;
            return;
        }
    }
}

