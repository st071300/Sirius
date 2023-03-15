#include <iostream>
#include "Poly.hpp"
#include "Quaternion.hpp"
#include <vector>
#include <complex>

using namespace std;

int main(int argc, char * argv[]) {

    std::vector<std::complex<double> > vc{{1., 2.},
                                          {3., 4.},
                                          {5., 6.}};
    Polynom <std::complex<double>> p1(vc);
    p1.PrintPolynom(std::cout);

    std::vector<double> vd{1., 2., 3.};
    Polynom<double> p2(vd);
    p2.PrintPolynom(std::cout);

    std::vector<double> vd2{4., 5., 6.};
    Polynom<double> p22(vd2);
    p22.PrintPolynom(std::cout);

    Polynom<double> new_p;
    new_p = p2 + p22;
    new_p.PrintPolynom(std::cout);


    std::vector<Quaternion < double> > vq{{1., 2., 3., 4.},
                                          {5., 6., 7., 8.}};
    Polynom <Quaternion
            <double>> p3(vq);
    p3.PrintPolynom(std::cout);

    return 0;
}

