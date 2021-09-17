#ifndef __SIMULATION__
#define __SIMULATION__
#include <cstdint>
#include "particle.hpp"

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
    float* dV) -> void;

auto computeMotion(
    particle* pa,
    const float* dV,
    const float& wave_speed,
    const float& dT,
    const float& lBox,
    const uint8_t& nDim,
    const bool& walls
    ) -> void;

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
    const bool& walls) -> void;

#endif