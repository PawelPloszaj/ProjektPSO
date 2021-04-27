#define _USE_MATH_DEFINES
#include <cmath>
#include <algorithm>
#include <functional>
#include <utility>
#include <iostream>
#include <vector>
#include <random>
#include <fstream>
#include <math.h>

class PSO
{
private:

    int numberOfParticles; // Liczba cz¹stek
    int numberOfDimensions; // Liczba wymiarów
    int maximumOfIteration;  // Maksymalna liczba iteracji

    double minRand; // Minimalna wartoœæ zmiennych losowanych
    double maxRand; // Maksymalna wartoœæ zmiennych losowanych
    double FunctionValueError; // Wartoœæ maksymalna b³êdu wyliczonej funkcji przed zakoñczeniem optymalizacji

    double W; // Wspó³czynnik bezw³adnoœci
    double C1; // Wspó³czynnik d¹¿enia do najlepszego lokalnego rozwi¹zania
    double C2; // Wspó³czynnik d¹¿enia do najlepszego globalnego rozwi¹zania

    int iteracja = 1;

    std::pair<std::vector<double>, double> result;

public:

    PSO(int NOP = 500, int NOD = 1, int MOI = 1000,
        double MIR = 0, double MAR = 640, double ERC = 0.01,
        double w = 0.9, double c1 = 0.6, double c2 = 1.4)
        : numberOfParticles(NOP),
        numberOfDimensions(NOD),
        maximumOfIteration(MOI),
        minRand(MIR),
        maxRand(MAR),
        FunctionValueError(ERC),
        W(w),
        C1(c1),
        C2(c2)
    {}

    auto set_numberOfParticles(int NOP) 
    { 
        this->numberOfParticles = NOP; 
    };

    auto set_numberOfDimension(int NOD) 
    { 
        this->numberOfDimensions = NOD;
    };

    auto set_maximumOfIteration(int MOI) 
    { 
        this->maximumOfIteration = MOI; 
    };

    auto set_minRand(double MIR) 
    { 
        this->minRand = MIR; 
    };

    auto set_maxRand(double MAR) 
    { 
        this->maxRand = MAR; 
    };

    auto set_errorCon(double ERC) 
    { 
        this->FunctionValueError = ERC;
    };

    auto set_w(double w) 
    { 
        this->W = w; 
    };

    auto set_c1(double c1) 
    { 
        this->C1 = c1; 
    };

    auto set_c2(double c2) 
    { 
        this->C2 = c2; 
    };

    auto optimize(std::function<double(std::vector<double>)> fitFunc)
    {

        std::mt19937 random(static_cast<unsigned int>(time(nullptr)));
        std::uniform_real_distribution<> rand1(minRand, maxRand);
        std::uniform_real_distribution<> rand2(0, 1);

        std::vector<double> localBest(numberOfDimensions);
        generate(localBest.begin(), localBest.end(), [&] { // Generowanie wartoœci localBest w zale¿noœci od iloœci wymiarów
            return rand1(random); 
            });

        std::vector<double> globalBest(numberOfDimensions);
        globalBest = localBest; // Przypisanie wartoœci localBest do globalBest

        std::vector <std::vector<double> > particles; // Deklaracja cz¹steczek
        particles.resize(numberOfParticles, std::vector<double>(numberOfDimensions, 0));

        std::vector<double> velocity; // Deklaracja prêdkoœci cz¹steczki
        velocity.resize(numberOfDimensions, 0);

        for (auto& row : particles)
        {
            std::vector<double> temp(numberOfDimensions);
            temp.clear();

            for (auto& col : row)
            {
                col = rand1(random);
                temp.push_back(col);
            }

            localBest.clear();
            localBest = temp;

            if (fitFunc(localBest) < fitFunc(globalBest))
            {
                globalBest.clear();
                globalBest = localBest;
            }
        }

        int iterator = 0;
        while (fitFunc(globalBest) > FunctionValueError and iterator < maximumOfIteration)
        {
            ++iterator;
            for (int i = 0; i < numberOfParticles; ++i)
            {
                for (int j = 0; j < numberOfDimensions; ++j)
                {
                    double r1 = rand2(random);
                    double r2 = rand2(random);

                    velocity[j] = W * velocity[j]
                        + (C1 * r1 * (localBest[j] - particles[i][j]))
                        + (C2 * r2 * (globalBest[j] - particles[i][j]));

                    particles[i][j] = particles[i][j] + velocity[j];
                }

                if (fitFunc(particles[i]) < fitFunc(localBest))
                {
                    localBest.clear();
                    localBest = particles[i];
                }
            }

            if (fitFunc(localBest) < fitFunc(globalBest))
            {
                globalBest.clear();
                globalBest = localBest;
            }

            // Logs
            PrintInfo(particles);

        }

        this->result = make_pair(globalBest, fitFunc(globalBest));
        return make_pair(globalBest, fitFunc(globalBest)); // Zwracanie wyników

    }

    void PrintInfo(std::vector <std::vector<double> > particles)
    {
        std::ofstream out;
        out.open("E:\\Dane\\dane.txt", std::ios::app);
        out << "ITERACJA:" << iteracja << std::endl;
        std::cout << "ITERACJA:" << iteracja << std::endl;
        for (int i = 0; i < particles.size(); i++)
        {
            out << "nr:" << i + 1 << " \tX:" << particles[i][0] << "\t Y:" << particles[i][1] << std::endl;
            std::cout << "nr:" << i + 1 << " \tX:" << particles[i][0] << "\t Y:" << particles[i][1] << std::endl;
        }
        iteracja++;
        std::cout << "------------------------------------" << std::endl;
        out << "------------------------------------" << std::endl;
        out.close();
    }

    void Results()
    {
        std::cout << std::endl;
        std::cout << "Najlepsze wyliczone rozwiazanie: " << result.first[0] << " , " << result.first[1] << std::endl;
        std::cout << "Przyblizony blad funkcji: " << result.second << std::endl;
        std::ofstream out;
        out.open("E:\\Dane\\dane.txt", std::ios::app);
        out << "Najlepsze wyliczone rozwiazanie: " << result.first[0] << " , " << result.first[1] << std::endl;
        out << "Przyblizony blad funkcji: " << result.second << std::endl;
        out << "------------------------------------" << std::endl << std::endl;;
        out.close();
    }

};