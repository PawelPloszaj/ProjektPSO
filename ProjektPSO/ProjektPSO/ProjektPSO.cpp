#include "PSO.h"

double AckleyFunction(std::vector<double> vals_inp)
{
    double x = vals_inp[0];
    double y = vals_inp[1];

    double obj_val = 20 + std::exp(1) - 20 * std::exp(-0.2 * std::sqrt(0.5 * (x * x + y * y))) - std::exp(0.5 * (std::cos(2 * M_PI * x) + std::cos(2 * M_PI * y)));

    return obj_val;
}

double BoothFunction(std::vector<double> vals_inp) 
{
    double p1 = pow(vals_inp[0] + 2 * vals_inp[1] - 7, 2);
    double p2 = pow(2 * vals_inp[0] + vals_inp[1] - 5, 2);

    double obj_val = p1 + p2;

    return obj_val;
}

double LeviFunction(std::vector<double> vals_inp)
{
    double x = vals_inp[0];
    double y = vals_inp[1];
    double pi = M_PI;

    double obj_val = std::pow(std::sin(3 * pi * x), 2) + std::pow(x - 1, 2) * (1 + std::pow(std::sin(3 * M_PI * y), 2)) + std::pow(y - 1, 2) * (1 + std::pow(std::sin(2 * M_PI * x), 2));

    return obj_val;
}

int main()
{
    PSO pso;

    int function, particles, iteration, minRand, maxRand;
    double FunctionValueError, W, C1, C2;

    // MENU
    std::cout << "1. Ackley Function" << std::endl;
    std::cout << "2. Booth Function" << std::endl;
    std::cout << "3. Levi Function" << std::endl;
    std::cout << "Wybierz numer funkcji: ";
    std::cin >> function;
    std::cout << "Podaj liczbe czastek: ";
    std::cin >> particles;
    std::cout << "Podaj liczbe max iteracji: ";
    std::cin >> iteration;
    std::cout << "Podaj minRand: ";
    std::cin >> minRand;
    std::cout << "Podaj maxRand: ";
    std::cin >> maxRand;
    std::cout << "Podaj maksymalny blad funkcji: ";
    std::cin >> FunctionValueError;
    std::cout << "Podaj wartosc wspolczynnika bezwladnosci: ";
    std::cin >> W;
    std::cout << "Podaj wartosc wspolczynnika dazenia do najlepszego lokalnego rowiazania: ";
    std::cin >> C1;
    std::cout << "Podaj wartosc wspolczynnika dazenia do najlepszego globalnego rozwiazania: ";
    std::cin >> C2;
    std::cout << "------------------------------------" << std::endl;

    // SET VALUES
    pso.set_numberOfParticles(particles); // Ilość cząstek
    pso.set_numberOfDimension(2); // Ilość wymiarów
    pso.set_maximumOfIteration(iteration); // Maksymalna liczba iteracji
    pso.set_minRand(minRand); // Minimalna wartość losowa
    pso.set_maxRand(maxRand); // Maksymalna wartość losowa
    pso.set_errorCon(FunctionValueError); // Wartość maksymalna błędu wyliczonej funkcji przed zakończeniem optymalizacji
    pso.set_w(W); // Wartość współczynnika bezwładności
    pso.set_c1(C1); // Wartość współczynnika dążenia do najlepszego lokalnego rozwiązania
    pso.set_c2(C2); // Wartość współczynnika dążenia do najlepszego globalnego rozwiązania

    if (function == 1)
        auto result = pso.optimize(AckleyFunction);
    else if(function == 2)
        auto result = pso.optimize(BoothFunction);
    else if (function == 3)
        auto result = pso.optimize(LeviFunction);

    pso.Results();

    return 0;
}
