#include <iostream>
#include <cmath>
#include <vector>
#include "random.h"
#include <iomanip>

// g++ -c random.cpp -o random.o
// g++ -o sample_new main.cpp random.o -std=c++0x

using namespace std;

double waveFunctionMany(double alpha, int numberOfParticles, double xs[]){
    double xSum;
    for(int n=0; n<numberOfParticles; n++){
        xSum += xs[n]*xs[n];
    }
    return exp(-alpha*xSum);
}


double localEnergyAnalyticalMany(double alpha, int numberOfParticles, double xs[]){
    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    double xSum2 = 0.0;

    for(int n=0; n<numberOfParticles; n++){
        xSum2 += xs[n]*xs[n];
    }
    return (-hbar*hbar/(2*m))*(-2*alpha*numberOfParticles + 4*alpha*alpha*xSum2) + 0.5*m*omega*xSum2;
}

double localEnergyNumericalMany(double alpha, int numberOfParticles, double xs[]){
    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    double energy = 0.0;
    double dx = 0.001;
    double dpsidx2 = 0.0;
    double xSum2 = 0.0;

    double x_less[numberOfParticles] = {};
    double x_more[numberOfParticles] = {};

    for(int jj=0; jj<numberOfParticles ;jj++){
        x_less[jj] = xs[jj];
        x_more[jj] = xs[jj];
        xSum2 += xs[jj]*xs[jj];
    }
    for(int n=0; n<numberOfParticles; n++){
        x_less[n] = xs[n]-dx;
        x_more[n] = xs[n]+dx;

        dpsidx2 += (waveFunctionMany(alpha, numberOfParticles,x_more)+waveFunctionMany(alpha, numberOfParticles,x_less))/(dx*dx);
    }
    
    dpsidx2 += -2*numberOfParticles*waveFunctionMany(alpha,numberOfParticles,xs)/(dx*dx);

    energy = (1.0/waveFunctionMany(alpha,numberOfParticles,xs))*(-hbar*hbar/(2*m)*dpsidx2) + 0.5*m*omega*xSum2;
    
    return energy;
}

void MonteCarloSampling(int maxVariations, int numberOfParticles, double stepSize, int numberOfDimensions){
    Random random;

    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    int numberOfMCcycles = 1e6;

    double positionOld[numberOfParticles] = {};
    double positionNew[numberOfParticles] = {};

    double waveFunctionOld = 0.0;
    double waveFunctionNew = 0.0;

    double deltaEnergy_ana = 0.0;
    double deltaEnergy_num = 0.0;

    double variance = 0.0;
    double error = 0.0;

    double alpha = 0.3;

    cout << "alpha" << "\t" <<"energy (ana): " << "\t" << "energy (num): " << "\t" << "variance (ana): " << "\t" << "error: " << endl;

    // Evaluating the energy for different alpha (i.e. different wavefunctions)
    for (int i=0; i < maxVariations; i++){
        alpha += 0.05;
        
        double energy_ana = 0.0;
        double energy2_ana = 0.0;
        double energy_num = 0.0;
        double energy2_num = 0.0;

        // Pick a position:
        for(int m=0; m<numberOfParticles;m++){
            positionOld[m] = stepSize*(random.nextDouble()-0.5);
        }

        waveFunctionOld = waveFunctionMany(alpha, numberOfParticles, positionOld);

        // Moving one particle at the time
        for(int whichParticle=0; whichParticle<numberOfParticles; whichParticle++){
            
            for (int j=0; j<numberOfMCcycles; j++){

                // Suggest a move to a new position:
                positionNew[whichParticle] = positionOld[whichParticle] + stepSize*(random.nextDouble()-0.5);
            
                waveFunctionNew = waveFunctionMany(alpha, numberOfParticles,positionNew);

                // Check if new position is accepted:
                if (random.nextDouble() <= (waveFunctionNew*waveFunctionNew)/(waveFunctionOld*waveFunctionOld)){
                    // If accepted, make the move:
                    for (int ii=0; ii<numberOfParticles;ii++){
                        positionOld[ii] = positionNew[ii];
                    }
                    waveFunctionOld = waveFunctionNew;
                }
                // Calculate the local energy (in position x with alpha) analytically and numerically:
                deltaEnergy_ana = localEnergyAnalyticalMany(alpha, numberOfParticles, positionOld);
                deltaEnergy_num = localEnergyNumericalMany(alpha, numberOfParticles, positionOld);
                // Calculating the sum of all local energies:
                energy_ana += deltaEnergy_ana;
                energy2_ana += deltaEnergy_ana*deltaEnergy_ana;

                energy_num += deltaEnergy_num;
                energy2_num += deltaEnergy_num*deltaEnergy_num; 
            }
        }
        // Evaluating the expectation value - the average of all local energies:
        energy_ana /= (numberOfMCcycles*numberOfParticles);
        energy2_ana /= (numberOfMCcycles*numberOfParticles);
        energy_num /= (numberOfMCcycles*numberOfParticles);
        energy2_num /= (numberOfMCcycles*numberOfParticles);

        // Calculating the variance:
        variance = energy2_ana - energy_ana*energy_ana;
        
        // Estimating an error:
        error = sqrt(variance/numberOfMCcycles);
        cout << alpha << "\t" << energy_ana << "\t" << energy_num << "\t"  << variance << "\t" << error << endl;

    }

}


int main() {
    Random random;

    int maxVariantions = 10;
    int numberOfParticles = 1;
    int numberOfDimensions = 2;
    double stepSize = 1;

    std::cout << std::setprecision(6) << std::fixed;

    MonteCarloSampling(maxVariantions, numberOfParticles, stepSize, numberOfDimensions);
    // MonteCarloSamplingOne(maxVariantions);

    


    return 0;
}