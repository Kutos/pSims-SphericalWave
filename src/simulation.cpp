#include <cstdint>
#include <iostream>
#include <chrono>
#include <cmath>
#include "omp.h"
#include <ctime>
#include <cstring>
#include "hdf5Helper.hpp"
#include "particle.hpp"
#include "simulation.hpp"

auto computeDV(
    particle* particles_array,
    const float& lBox,
    const float& dX,
    const float& dT,
    const uint16_t& p,
    const uint16_t& step,
    const uint16_t& nParticles,
    const uint16_t& nWave,
    const uint8_t& nNeighbors,
    const uint8_t& nDim,
    float* dV) -> void
{
    float calc0;
    int8_t V_sign;
    float arr_X[nDim];
    uint16_t k{(step <= nWave) * step + (step > nWave) * nWave};



    for(uint8_t d{0}; d < nDim; d++) dV[d] = 0.f;
    for(uint16_t j{0}; j < k; j++)
    {
        for(uint16_t i{0}; i < nParticles; i++)
        {
            V_sign = (particles_array[i].getChargeSign() == (particles_array[p].getChargeSign())) * 2 - 1;
            for(uint8_t d{0}; d < nDim; d++)
            {
                for(uint8_t dp{0}; dp < nDim; dp++)
                    arr_X[dp] = particles_array[p].getPosition(dp);
                arr_X[d] -= 0.5 * dX;
                calc0 = - particles_array[i].calWaveEnergieValue(j, arr_X, lBox, nNeighbors);
                arr_X[d] += dX;
                calc0 += particles_array[i].calWaveEnergieValue(j, arr_X, lBox, nNeighbors);
                #pragma omp critical
                dV[d] += (calc0 * V_sign) / dX;
            }
        }
        particles_array[p].propagateWave(dT, j);
    }
}


auto computeMotion(
    particle* pa,
    const float* dV,
    const float& wave_speed,
    const float& dT,
    const float& lBox,
    const uint8_t& nDim,
    const bool& walls
    ) -> void
{
    float temp = 1. / pa->getMass();

    float new_X[nDim];
    float new_V[nDim];

    for(uint8_t i{0}; i < nDim; i++)
    {
        new_V[i] = (temp * dV[i]) * dT + pa->getVelocity(i);
        new_X[0] += new_V[i] * new_V[i];
    }

    new_X[1] = wave_speed * wave_speed;
    new_X[0] = 1.f + (new_X[0] / new_X[1]);
    temp = 1.f / sqrt(new_X[0]);
    
    for(uint8_t i{0}; i < nDim; i++)
    {
        new_V[i] *= temp;
        new_X[i] = pa->getVelocity(i) * dT + pa->getPosition(i);
    }

    temp = .5f * lBox;

    if(walls)
    {
        for(uint8_t i{0}; i < nDim; i++)
        {
            if(new_X[i] < - temp)
            {
                new_X[i] = - temp - std::fmod(new_X[i], temp);
                new_V[i] *= -1;

            }
            else if(new_X[i] > temp)
            {
                new_X[i] = temp - std::fmod(new_X[i], temp);
                new_V[i] *= -1;
            }
        }

    }else{
        for(uint8_t i{0}; i < nDim; i++)
        {
            if(new_X[i] < - temp) new_X[0] = temp + std::fmod(new_X[i], temp);
            else if(new_X[i] > temp) new_X[0] = - temp + std::fmod(new_X[i], temp);
        }
    }
    
    // Updating the position in particles list
    pa->setPosition(new_X);
    // Updating the velocity in particles list
    pa->setVelocity(new_V);
}


auto simulate(
    const float& lBox,
    const float& dX,
    const float& dT,
    const float& wave_speed,
    const float& wave_width,
    const float& wave_range,
    const uint16_t& nStep,
    const uint16_t& nParticles,
    const uint8_t& nDim,
    const bool& walls) -> void
{
    std::time_t t = std::time(0);   // get time now
    struct tm * now = localtime( & t );
    char buffer [80];
    strftime (buffer,80,"results-%Y-%m-%d-%H-%M.h5",now);
    hdf5Helper output(buffer, lBox, dX, dT, wave_speed, wave_width, wave_range, nStep, nParticles, nDim, walls);
    particle particles_array[nParticles];
    uint16_t nWave;
    {
        // The maximum waves to store in the array
        const uint16_t nWave_max{(3.f*wave_range)/(wave_speed*dT) + 1};
        // The effective waves amount to store in the waves array
        nWave = nWave_max * (nWave_max < (nStep))+ (nStep) * (nWave_max > (nStep));
        srand((long int) time(NULL));
        bool chargeSign = rand();
        #pragma omp parallel for default(none) \
            shared(particles_array, lBox, wave_speed, wave_width, wave_range, \
                    nDim, nWave, chargeSign, nParticles)
        for(uint16_t p = 0; p < nParticles; p++)
        {
            particles_array[p] = particle(lBox, 1.f, &wave_speed, &wave_width, &wave_range,
                        &nDim, &nWave, chargeSign);
        }
    }

    {
        // Number of neighbor box
        float dV[nDim];
        const uint8_t nNeighbors{(1u - walls) * ((6.f*wave_range)/(lBox) + 1u)};
        uint16_t i;
        for(uint16_t step = 0; step < (nStep); step++)
        {
            std::cout << "Step: " << step << std::endl;
            i = step % nWave;

            #pragma omp parallel default(none) \
                shared (particles_array, lBox, dX, dT, i,\
                        step, nParticles, nWave, nNeighbors,\
                        walls, nDim, dV, wave_speed)
            {
                #pragma omp for
                for(uint16_t p = 0; p < nParticles; p++) 
                    particles_array[p].resetWave(i);
        
                #pragma omp for
                for(uint16_t p = 0; p < nParticles; p++)
                {
                    computeDV(particles_array, lBox, dX, dT, p,
                                step, nParticles, nWave, nNeighbors, nDim, dV);
                    computeMotion(&(particles_array[p]), dV, wave_speed, dT, lBox, nDim, walls);
                }
            }
            output.writeStep(particles_array, nParticles, step, nDim);
        }

    }

}