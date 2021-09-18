#include "particle.hpp"
#include "wave.hpp"
#include <iostream>
#include <cstdint>
#include <cmath>
#include <chrono>
#include "omp.h"

particle::particle() 
    : mass(0.f), nWave(nullptr), nDim(nullptr), chargeSign(false),
    position(nullptr), velocity(nullptr), waves_list(nullptr){}

particle::particle(const particle& p) 
    : mass(p.mass), nWave(p.nWave), nDim(p.nDim), chargeSign(p.chargeSign)
{
    position = new float(*p.position);
    velocity = new float(*p.velocity);
    waves_list = new wave(*p.waves_list);
}

particle::particle(particle&& p) 
    : mass(p.mass), nWave(p.nWave), nDim(p.nDim), chargeSign(p.chargeSign),
    position(p.position), velocity(p.velocity), waves_list(p.waves_list)
{
    p.nDim = nullptr;
    p.nWave = nullptr;
    p.position = nullptr;
    p.velocity = nullptr;
    p.waves_list = nullptr;
}

particle::particle(const float& lBox, const float& mass,
    const float* wave_speed, const float* wave_width, const float* wave_range,
    const uint8_t* nDim, const uint16_t* nWave, const bool& chargeSign)  
    : mass(mass), nDim(nDim), nWave(nWave), chargeSign(chargeSign)
{
    const float lBoxd5 = .5 * lBox;
    uint seed;
    position = new float[*nDim];
    velocity = new float[*nDim];
    waves_list = new wave[*nWave];

    srand(static_cast<long int>(time(NULL) * omp_get_num_threads()));

    for(uint8_t d{0}; d < *nDim; d++)
    {
        seed = rand() * (d+1);
        position[d] = float{rand_r(&seed) % int{lBox * 1000.f}} / 1000.f - lBoxd5;
    }

    for(uint16_t w{0}; w < *nWave; w++)
    {
        waves_list[w] = wave(position, wave_speed, wave_width, wave_range, nDim);
    }

}

particle::~particle()
{
    nDim = nullptr;
    nWave = nullptr;
    delete[] position;
    delete[] velocity;
    delete[] waves_list;
}

particle& particle::operator=(const particle& p)
{
    if(this == &p) return *this;
    this->~particle();

    mass = p.mass;
    chargeSign = p.chargeSign;
    nDim = p.nDim;
    nWave = p.nWave;
    position = new float(*p.position);
    velocity = new float(*p.velocity);
    waves_list = new wave(*p.waves_list);

    return *this;
}

particle& particle::operator=(particle&& p)
{
    if(this == &p) return *this;
    this->~particle();

    mass = p.mass;
    chargeSign = p.chargeSign;
    nDim = p.nDim;
    nWave = p.nWave;
    position = p.position;
    velocity = p.velocity;
    waves_list = p.waves_list;

    p.position = nullptr;
    p.velocity = nullptr;
    p.waves_list = nullptr;

    return *this;
}

auto particle::setPosition(const float* pos) -> void
{
    for(uint8_t d{0}; d < *nDim; d++)
        position[d] = pos[d];
}
auto particle::setVelocity(const float* vel) -> void
{
    for(uint8_t d{0}; d < *nDim; d++)
        velocity[d] = vel[d];
}

auto particle::resetWave(const uint16_t& indWave) -> void
{
    waves_list[indWave].resetCenter(position);
    waves_list[indWave].resetRadius();
}

auto particle::propagateWave(const float& dT, const uint16_t& indWave) -> void
{
    waves_list[indWave].propagate(dT);
}

auto particle::calWaveEnergieValue(const uint16_t& indWave, const float* position, const float& lBox, const uint8_t& nNeighbors) -> float
{
    switch (*nDim)
    {
    case 1:
        return waves_list[indWave].calEnergieValue1D(position, lBox, nNeighbors);
        break;
    case 2:
        return waves_list[indWave].calEnergieValue2D(position, lBox, nNeighbors);
        break;
    case 3:
        return waves_list[indWave].calEnergieValue3D(position, lBox, nNeighbors);
        break;

    default:
        throw "ERROR DIM";
        break;
    }
}