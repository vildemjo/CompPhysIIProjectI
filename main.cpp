#include <iostream>
#include <cmath>
#include <vector>
#include "random.h"
#include <iomanip>

// g++ -c random.cpp -o random.o
// g++ -o sample_new main.cpp random.o -std=c++0x

using namespace std;

const int numberOfDimensions = 3;
const int numberOfParticles = 2;

double waveFunctionMany(double alpha, double r[numberOfParticles][numberOfDimensions]){
    
    double rSum = 0.0;

    for(int i1=0; i1<numberOfParticles; i1++){
        for(int n1=0; n1<numberOfDimensions; n1++){
            rSum += r[i1][n1]*r[i1][n1];
        }
    }
    return exp(-alpha*rSum);
}


double localEnergyAnalyticalMany(double alpha, double r[numberOfParticles][numberOfDimensions]){
    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    double rSum2 = 0.0;

    for(int i2=0; i2<numberOfParticles; i2++){
        for(int n2=0; n2<numberOfDimensions; n2++){
            rSum2 += r[i2][n2]*r[i2][n2];
        }
        
    }
    return (-hbar*hbar/(2*m))*(-2*alpha*numberOfParticles*numberOfDimensions + 4*alpha*alpha*rSum2) + 0.5*m*omega*rSum2;
}

double localEnergyNumericalMany(double alpha, double r[numberOfParticles][numberOfDimensions]){
    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    double energy = 0.0;
    double step = 0.001;
    double dpsidr2 = 0.0;
    double rSum2 = 0.0;

    double r_less[numberOfParticles][numberOfDimensions] = {0.0};
    double r_more[numberOfParticles][numberOfDimensions] = {0.0};

    // Making copies of r
    for(int i3=0; i3<numberOfParticles; i3++){
        for(int n3=0; n3<numberOfDimensions;n3++){
            r_less[i3][n3] = r[i3][n3];
            r_more[i3][n3] = r[i3][n3];

            rSum2 += r[i3][n3]*r[i3][n3];
        }
    }

    // Changing the copies to obtain an array of "next step" and "previous step"
    for(int i4=0; i4<numberOfParticles; i4++){
         for(int n4=0; n4<numberOfDimensions;n4++){
            // finding the next step and previous step for one particle in one spesific dimension
            r_less[i4][n4] = r[i4][n4]-step;
            r_more[i4][n4] = r[i4][n4]+step;

            dpsidr2 += (waveFunctionMany(alpha, r_more)+waveFunctionMany(alpha,r_less))/(step*step);

            // Resetting the arrays so that a new particle and spesific dimension can be calculated
            r_less[i4][n4] = r[i4][n4];
            r_more[i4][n4] = r[i4][n4];
        }
    }
    
    dpsidr2 += -2*numberOfParticles*numberOfDimensions*waveFunctionMany(alpha,r)/(step*step);

    energy = (1.0/waveFunctionMany(alpha, r))*(-hbar*hbar/(2*m)*dpsidr2) + 0.5*m*omega*rSum2;
    
    return energy;
}

void MonteCarloSampling(int maxVariations, double stepSize){
    Random random;

    double hbar = 1.0;
    double m = 1.0;
    double omega = 1.0;

    int numberOfMCcycles = 1e6;

    double positionOld[numberOfParticles][numberOfDimensions] = {0.0};
    double positionNew[numberOfParticles][numberOfDimensions] = {0.0};

    double waveFunctionOld = 0.0;
    double waveFunctionNew = 0.0;

    double deltaEnergy_ana = 0.0;
    double deltaEnergy_num = 0.0;

    double variance = 0.0;
    double error = 0.0;

    double alpha = 0.3;

    cout << "alpha" << "\t" <<"energy (ana): " << "\t" << "energy (num): " << "\t" << "variance (ana): " << "\t" << "error: " << endl;

    // Evaluating the energy for different alpha (i.e. different wavefunctions)
    for (int j=0; j < maxVariations; j++){
        alpha += 0.05;
        
        double energy_ana = 0.0;
        double energy2_ana = 0.0;
        double energy_num = 0.0;
        double energy2_num = 0.0;

        // Pick a position:
        for(int i5=0; i5<numberOfParticles;i5++){
            for(int n5=0; n5<numberOfDimensions; n5++){
                positionOld[i5][n5] = stepSize*(random.nextDouble()-0.5);
            }
        }

        waveFunctionOld = waveFunctionMany(alpha, positionOld);

        // Moving one particle at the time
        for(int whichParticle=0; whichParticle<numberOfParticles; whichParticle++){
            
            for (int jj=0; jj<numberOfMCcycles; jj++){

                // Suggest a move to a new position:
                for(int n6=0; n6<numberOfDimensions; n6++){
                    positionNew[whichParticle][n6] = positionOld[whichParticle][n6] + stepSize*(random.nextDouble()-0.5);
                }
                waveFunctionNew = waveFunctionMany(alpha, positionNew);

                // Check if new position is accepted:
                if (random.nextDouble() <= (waveFunctionNew*waveFunctionNew)/(waveFunctionOld*waveFunctionOld)){
                    
                    // If accepted, make the move:
                    for(int n7=0; n7<numberOfDimensions;n7++){
                        positionOld[whichParticle][n7] = positionNew[whichParticle][n7];
                    }
                    
                    waveFunctionOld = waveFunctionNew;
                }
                // Calculate the local energy (in position x with alpha) analytically and numerically:
                deltaEnergy_ana = localEnergyAnalyticalMany(alpha, positionOld);
                deltaEnergy_num = localEnergyNumericalMany(alpha, positionOld);
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

    double stepSize = 1;
 
    

    std::cout << std::setprecision(6) << std::fixed;

    MonteCarloSampling(maxVariantions, stepSize);
    // MonteCarloSamplingOne(maxVariantions);



    return 0;
}