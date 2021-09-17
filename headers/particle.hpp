#ifndef __PARTICLES__
#define __PARTICLES__
#include <cstdint>
#include "wave.hpp"
class particle
{
private:
    const uint8_t* nDim;
    const uint16_t* nWave;

    wave* waves_list;
    float* position;
    float* velocity;

    float mass;
    bool chargeSign;

public:
    particle();
    particle(const particle&);
    particle(particle&&);
    particle(const float& lBox, const float& mass,
    const float* wave_speed, const float* wave_width, const float* wave_range,
    const uint8_t* nDim, const uint16_t* nWave, const bool& chargeSign);
    ~particle();

    particle& operator=(const particle& p);
    particle& operator=(particle&& p);

    inline auto getMass() const -> float{return mass;}
    inline auto getChargeSign() const -> bool{return chargeSign;}
    inline auto getPosition(const uint8_t d) const -> float{return position[d];}
    inline auto getVelocity(const uint8_t d) const -> float{return velocity[d];}

    auto setPosition(const float* position) -> void;
    auto setVelocity(const float* velocity) -> void;

    auto calWaveEnergieValue(const uint16_t& indWave, const float* position, const float& lBox, const uint8_t& nNeighbors) -> float;
    auto propagateWave(const float& dT, const uint16_t& indWave) -> void;
    auto resetWave(const uint16_t& indWave) -> void;


};

#endif