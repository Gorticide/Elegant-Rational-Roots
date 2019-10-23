/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
 *
 *    I N T R O D U C T O R Y   A N A L Y S I S
 *
 *    October 2019:  Chapter 5, Section 8 Computer Exercises
 *
 *   3. Write a computer program to find and print all possible
 *      rational roots of
 *
 *    P(x) =  a_0*x^n + a_1*x^(n-1) + a_2*x^(n-2) + ... + a_(n-1)*x + a_n = 0
 *
 *      where the user supplies the coefficients of P(x).
 *
 *   5.  Modify program of Exercise 3 to find and test all possible
 *       rational roots of P(x) = 0 and report if any are actual
 *       roots of the equation.
 *
 *       Use the program of Exercise 1 that tests for an upper bound
 *       of real roots to eliminate unnecessary tests:
 *
 *    Also find the greatest negative integer that is a lower bound
 *    for the real roots of P(x).
 *
 *    To do this, create P(-x) or -P(x), depending on which one has positive a_0,
 *    then find the least positive integer that is an upper bound for the real roots of
 *    P(-x) or -P(x).
 *
 *       Note:  Check for multiple roots!
 *       This is what will distinguish ratroots2 from original
 *       ratroots.cpp from algebra2/polynomial_factoring.
 *
 *   Michael William Hentrich
 *   original ratroots.cpp:  2 November 2018
 *   Improved on 27 February 2019 using Fraction class
 *   ratroots_synthDiv.cpp: October 2019
 *          - still going to use Fraction class for Posterity
 *          - and show multiplicity of roots
 *          - not evaluating P(x) to test P(x) == 0
 *          - but using synthetic division
 *
 *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
 */


 #include "fraction.h"
 #include <iostream>
 #include <cmath>
 #include <vector>
 #include <string>
 #include <sstream>      // for print_poly
 #include <utility>      // for std::pair
 #include <algorithm>    // std::for_each

 typedef std::vector<long> Integers;
 typedef std::vector<Fraction> Rationals;

 Integers build_Polynomial(long n);
 std::string print_poly(const Integers& P, long n);
 long get_Upper_Bound(const Integers& T, int degree);
 long get_Lower_Bound(const Integers& T, long degree);
 std::pair<long, long> get_Bounds(const Integers& T, long degree);
 Integers get_factors (long num);

 std::pair<Fraction, Integers>
 synthetic_division(Fraction M, Integers A, long NP);

 Integers Rationals_as_Integers(Rationals Q);
 Rationals Integers_as_Rationals(Integers Z);
 void Test_Rational_Candidates(Integers p, Integers q,
                              std::pair<long, long> BOUNDS,
                              Integers poly,
                              Rationals &solutions,
                              long &slns);

 int main()  {

   long itr_Roots(0);

   // May replace with
   // factors_num_n.size() and factors_den_0.size()
   long N_num_n(0);  // for factors of constant a-n
   long N_den_0(0);  // for factors of leading coefficient a_0

   std::cout << "\nFirst, enter the degree of the polynomial : ";
   long degree;
   std::cin >> degree;

   Rationals roots(degree);  // just the rational roots, at most, degree
   Integers Poly( build_Polynomial(degree) );

   // Since there can be a maximum of a_n positive factors in factors_num_n:
   long a_n = (long)abs(Poly.at(0));
   Integers factors_num_n(a_n);  // a_n is the constant at Poly[0]
   for (long i = 0; i < a_n; ++i) factors_num_n.at(i) = (long)0;

   // Since there can be a maximum of a_0 positive factors in factors_den_0:
   long a_0 = (long)abs(Poly.at(degree));
   Integers factors_den_0(a_0);
   // a_0 is leading coefficient of highest degree term
   for (long i = 0; i < a_0; ++i) factors_den_0.at(i) = (long)0;

   //Display polynomial P(x)
   std::cout << "\nP(x) = " << print_poly(Poly, degree) << "\n\n";

   std::pair<long, long> B = get_Bounds(Poly, degree);
   std::cout << "P(x) has upper bound of " << B.second;
   std::cout << "\n\nP(x) has lower bound of " << B.first;
   std::cout << "\n\n";


   // Store all possible [positive] numerators and [positive] denominators
   for (long j = 1; j <= a_n; j++) {
      if ( (a_n % j) != 0 ) continue;
      factors_num_n.at(N_num_n) = j;
      N_num_n++;
   }
   factors_num_n.erase(remove(factors_num_n.begin(), factors_num_n.end(),
                               (long)0 ),
                        factors_num_n.end() );
   factors_num_n.shrink_to_fit();

   std::cout << "\nThe factors of constant " << a_n << ": { ";
   for (auto f : factors_num_n)  std::cout << f << ' ';
   std::cout << "}\n\n";

   for (long j = 1; j <= a_0; j++) {
      if ( (a_0 % j) != 0) continue;
      factors_den_0.at(N_den_0) = j;
      N_den_0++;
   }
   factors_den_0.erase(remove(factors_den_0.begin(), factors_den_0.end(),
                               (long)0 ),
                        factors_den_0.end() );
   factors_den_0.shrink_to_fit();

   std::cout << "The factors of leading coefficient " << a_0 << ": { ";
   for (auto f : factors_den_0)  std::cout << f << ' ';
   std::cout << "}\n";

   Test_Rational_Candidates(factors_num_n, factors_den_0,
                            B, Poly, roots, itr_Roots);

   if (itr_Roots != 0) {
     roots.erase(remove(roots.begin(), roots.end(), (long)0 ),
                        roots.end() );
     roots.shrink_to_fit();
     std::cout << "\nThe Rational Roots of P(x) = 0 are ";
     std::cout << "{ ";
     for (auto r : roots)  {
       std::cout << r << ' ';
     }
     std::cout << "}\n\n";
   }
   else std::cout << "\nThere are no rational roots.\n\n";
}

long get_Upper_Bound(const Integers& T, long degree)   {
   bool stop = false;
   long Last_M(1), N(degree);

   for (long M = 1; M < abs(T.at(0)); M++)  {
       if (stop) break;
       long Q = T.at(N);

       for (long i = N; i > 0; --i)  {

          long P = Q*M;
          Q = T.at(i - 1) + P;

          if (Q < 0)  {
              Last_M = M + 1;
              break;
           }
          if (i == 1) stop = true;
       }
   }
   return Last_M;
}

long negate_Indeterminate_X(long C, long exp)  {
    if (exp%2 == 1) return -C;
    else return C;
}

Integers negate_Polynomial(const Integers& T, long degree)  {

   Integers polynomial(degree+1);

   for (long i = degree; i > 0; --i)  polynomial.at(i) = -T.at(i);
   polynomial.at(0) = -T.at(0);
   return polynomial;
}

long get_Lower_Bound(const Integers& T, long degree)  {

  // First, create P(-x)
  Integers polynomial(degree+1);

  for (long i = degree; i > 0; --i)  {
     polynomial.at(i) = negate_Indeterminate_X(T.at(i), i);
   }
   polynomial.at(0) = negate_Indeterminate_X(T.at(0), 0);

    // Next, check that a_0 > 0 for P(-x)

    if (polynomial.at(degree) < 0) {
        polynomial = negate_Polynomial(polynomial, degree);
        //Display polynomial -P(x)
        std::cout << "-P(-x) = ";
    }
    else {  // else Use P(-x)
        //Display polynomial P(x)
        std::cout << "P(-x) = ";
    }
    std::cout << print_poly(polynomial, degree) << "\n\n";
    return -(long) get_Upper_Bound(polynomial, degree);
}


std::string print_poly(const Integers& T, long n) {
   std::ostringstream out;
   long idx(n), degree(n);  // idx and degree are intitialized
                                 // with n (the "degree" of T)
//  if (zero_factor) out << "(x)*(";
                                                 // Using Lambda (C++11)
   std::for_each(T.rbegin(), T.rend(),   // range
            [&out, &idx, &degree](long t) {             // operation
              char sign;
              std::string val, expon;
              if (t < 0)  {
                sign = '-';
                val = std::to_string(-t);
              }
              else  {
                sign = '+';
                val = std::to_string(t);
              }
              if (t != 0) {
                if (idx != degree) out << sign;  // nasty bug?  NO!
                // We know sign is positive for leading coefficient
                // where T[degree] is this number, by definition positive
                // so we do not print sign

                if (idx != 0)  {
                  if (idx != 1)  expon = " x^" + std::to_string(idx) + ' ';
                  else expon= " x " ;
                }
                // when coefficients are 0: Indeterminate 'x' not displayed
                else expon= ' ';

                if  ( ( (abs(t) != 1) && (idx != 0) )
                     ||
                      (idx == 0)
                    )
                {
                    out << ' ';
                    out << val;
                }
                // else // t = 1 or t = -1 so leave out the space
                out << expon;

              }   // else t == 0 and no term is displayed
              // nor do we print Indeterminate 'x'
              // when the coefficient, t,  is 0.

              idx -= 1;                // decrement index, idx
           });  // end for_each loop

  // if (zero_factor) out << ") ";
   out << "== 0";
   return out.str();
}


Integers build_Polynomial(long n)   {
  long degree = n;
  long coefficient;

  Integers polynomial(degree+1);

  for (long i = degree; i > 0; --i)  {
     std::cout << "Enter (INTEGRAL) coefficient for x^" << i << " : ";
     std::cin >> coefficient;
     polynomial.at(i) = coefficient;
   }
   std::cout << "\nconstant term: ";
   std::cin >> coefficient;
   polynomial.at(0) = coefficient;
   return polynomial;
}

 std::pair<long, long> get_Bounds(const Integers& T, long degree)  {
     std::pair<long, long> B;
     long UB = get_Upper_Bound(T, degree);
     long LB = get_Lower_Bound(T, degree);
     B.first = LB;
     B.second = UB;
     return B;
 }


 Integers  get_factors (long num)  {
   num = abs(num);
   Integers list;
   long divisor = 0;

    do {
     divisor += 1;
     if (num%divisor == 0) list.push_back(divisor);
   }  while (divisor != num);

   return list;
 }


// Transform Integers into Fractions for Computations
Rationals Integers_as_Rationals(Integers Z)  {
    Rationals Q;
    for (auto INTEGER : Z)  {
      Q.push_back(Fraction(INTEGER, 1));
    }
    return Q;
  }



// For printing Polynomials with Rational Coeffs as Integers
Integers Rationals_as_Integers(Rationals Q)  {
   Integers Z;
   for (auto FRACTION : Q)  {
     Z.push_back(FRACTION.getNum());
   }
   return Z;
 }


std::pair<Fraction, Integers>
synthetic_division(Fraction M, Integers A, long NP)  {  // NP is degree of A

    std::pair<Fraction, Integers> depressed;
    long NQ = NP - 1;  // NQ is degree of Q
    Rationals Q(NQ+1);
    Rationals F( Integers_as_Rationals(A) );
    Q.at(NQ) = F.at(NP);
    if (NQ == 0)  {
      depressed.first = F.at(0) + Q.at(0)*M;
      depressed.second = A;
      return depressed;
    }
    for (long i = NQ; i > 0; --i)  {
      Fraction P = Q.at(i)*M;
      Q.at(i - 1) = F.at(i) + P;
    }
    depressed.first = F.at(0) + Q.at(0)*M;
    Integers Q_Integers( Rationals_as_Integers(Q) );
    //Q_Integers = Rationals_as_Integers(Q);
    depressed.second = Q_Integers;
    return depressed;
}


void Test_Rational_Candidates(Integers p, Integers q,
                             std::pair<long, long> BOUNDS,
                             Integers A,
                             Rationals &solutions,
                             long &slns)  {

    Integers poly(A);
    for (auto h : p)  {
      for (auto k : q)  {
         for (long S = -1; S < 2; S += 2)  {
           Fraction T = Fraction(S, 1) * Fraction(h,k);
           if (
             (T > Fraction(BOUNDS.second, 1) )
             ||
             (T < Fraction(BOUNDS.first, 1))
           )
           {
             continue;
             // advance to next S (-T)
             // OR next possible rational root h/k
           }
           else  {  // T is in interval [LB, UB]

// * Use Synthetic Division to see if T is a rational root *

            long NP = poly.size()- 1;
            long NQ = NP - 1;
            Integers Q(NQ+1);
            std::pair<Fraction, Integers>
            Results = synthetic_division(T, poly, NP);
            Q = Results.second;  // created in above

        //std::cout << "\nTesting candidate " << T << '\n';
        if (Results.first == Fraction(0, 1))  {
          if (slns != 0) {
            for (long i = 0; i < slns; ++i)  {
              if (solutions.at(i) != T)   {
                 solutions.at(slns) = T;
                 slns++;
                 std::cout << "\nRoot " << slns << " = "
                           << T << '\n';
              }
              else continue;
           }
        }
        else  // solution set is { }
          {
           solutions.at(slns) = T;  // roots.at(0) = T
           slns++;
           std::cout << "\nRoot " << slns << " = "
                     << T << '\n';
          }

  // *  Depress P(x) and continue search:

      if (NQ != 0) {
           NP = NQ;
           for (long i = 0; i < NQ + 1; i++)  {
             poly.at(i) = Q.at(i);  // depress equation
           }
           poly.pop_back();
           Results = synthetic_division(T, poly, NP);
           Q = Results.second;

           if (Results.first == Fraction(0, 1)) {
            solutions.at(slns) = T;
            slns++;
            std::cout << "\nRoot " << slns << " = "
                      << T << '\n';
           }
          else {
              std::cout << "\nRational Number " << T
                        << " is not a multiple root ...\n";
              if (NQ == 0)  continue;
          }
        }
        else continue;
        }
        else continue;  // T is not a rational root
       }  // end else (T is in interval)
      }  // end S
    }  // end k
   } // end h
 }  // end function
