#ifndef __WAVES__
#define __WAVES__
#include <cstdint>
class wave
{
private:
    const float m = 0.3606737602222408518399811702504730343566614885382464835338623517;
    const float* wave_speed;
    const float* wave_width;
    const float* wave_range;
    const uint8_t* nDim;

    float* center;
    float radius;

public:
    wave();
    wave(const wave&);
    wave(wave&&);
    wave(const float* center, const float* wave_speed, const float* wave_width, const float* wave_range, const uint8_t* nDim);
    ~wave();

    wave& operator=(const wave&);
    wave& operator=(wave&&);

    auto resetCenter(const float* ctr) -> void;
    inline auto resetRadius() -> void{radius = 0.f;};
    inline auto propagate(const float& dT) -> void{radius += *wave_speed * dT;};

    //auto calEnergieValue(const float* position, const float& lBox, const uint8_t& nNeighbors) -> float;
    auto calEnergieValue1D(const float* position, const float& lBox, const uint8_t& nNeighbors) -> float;
    auto calEnergieValue2D(const float* position, const float& lBox, const uint8_t& nNeighbors) -> float;
    auto calEnergieValue3D(const float* position, const float& lBox, const uint8_t& nNeighbors) -> float;

};

#endif