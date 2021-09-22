#include <cstdint>
#include <cmath>
#include "wave.hpp"

wave::wave()
    : center(nullptr), wave_speed(nullptr), wave_width(nullptr), wave_range(nullptr){}

wave::wave(const wave& w) 
    : wave_speed(w.wave_speed), wave_width(w.wave_width), wave_range(w.wave_range), radius(w.radius), nDim(nDim)
{
    center = new float[*nDim];
    for(uint8_t d = 0; d < *nDim; d++)
        center[d] = w.center[d];
}
wave::wave(wave&& w) 
    : wave_speed(w.wave_speed), wave_width(w.wave_width), wave_range(w.wave_range), radius(w.radius), nDim(nDim)
{
    center = w.center;
    w.center = nullptr;
}

wave::wave(const float* ctr, const float* wave_speed, const float* wave_width, const float* wave_range, const uint8_t* nDim)
    : wave_speed(wave_speed), wave_width(wave_width), wave_range(wave_range), radius(0.f), nDim(nDim)
{
    center = new float[*nDim];
    for(uint8_t d = 0; d < *nDim; d++)
        center[d] = ctr[d];
}
wave::~wave()
{
    nDim = nullptr;
    wave_speed = nullptr;
    wave_width = nullptr;
    wave_range = nullptr;
    delete[] center;
}


wave& wave::operator=(const wave& w)
{
    if(this == &w) return *this;
    this->~wave();

    radius = w.radius;
    nDim = w.nDim;
    wave_speed = w.wave_speed;
    wave_width = w.wave_width;
    wave_range = w.wave_range;
    center = new float(*w.center);

    return *this;
}

wave& wave::operator=(wave&& w)
{
    if(this == &w) return *this;
    this->~wave();

    radius = w.radius;
    nDim = w.nDim;
    wave_speed = w.wave_speed;
    wave_width = w.wave_width;
    wave_range = w.wave_range;
    center = w.center;

    w.center = nullptr;

    return *this;
}


auto wave::resetCenter(const float* ctr) -> void
{
    for(uint8_t d = 0; d < *nDim; d++)
        center[d] = ctr[d];
}

auto wave::calEnergieValue1D(const float* position, const float& lBox, const uint8_t& nNeighbors) -> float
{
    if(*nDim != 1) throw "ERROR DIM";

    float rValue = 0.;
    float a, b;
    float invRange2{1.f / ( (*wave_range )* (*wave_range))};
    float width_factor{1.f / ((*wave_width) * (*wave_width) * m)};

    for(int8_t i = - static_cast<int8_t>(nNeighbors); i < nNeighbors + 1; i++)
    {
        a = (position[0] - (center[0] - i * lBox)) - radius;
        b = a * invRange2;
        a = a * (-a * width_factor) - b;
        b = exp(a);
        rValue += b;

    }
    return rValue;
}

auto wave::calEnergieValue2D(const float* position, const float& lBox, const uint8_t& nNeighbors) -> float
{
    if(*nDim != 2) throw "ERROR DIM";

    float rValue = 0.;
    float a, b, c;
    float invRange2{1.f / ( (*wave_range )* (*wave_range))};
    float width_factor{1.f / ((*wave_width) * (*wave_width) * m)};

    for(int8_t i = - static_cast<int8_t>(nNeighbors); i < nNeighbors + 1; i++)
    {
        a = (position[0] - (center[0] - i * lBox));
        a *= a;
        for(int8_t j = - static_cast<int8_t>(nNeighbors); j < nNeighbors + 1; j++)
        {
            b = (position[1] - (center[1] - j * lBox));
            b = (b * b) + a;
            c = sqrt(b) - radius,
            b *= invRange2;
            c = c * (- c * width_factor) - b;
            b = exp(c);
            rValue += b;
        }
    }
    return rValue;
}

auto wave::calEnergieValue3D(const float* position, const float& lBox, const uint8_t& nNeighbors) -> float
{
    if(*nDim != 3) throw "ERROR DIM";

    float rValue = 0.;
    float a, b, c, d;
    float invRange2{1.f / ( (*wave_range )* (*wave_range))};
    float width_factor{1.f / ((*wave_width) * (*wave_width) * m)};

    for(int8_t i = - static_cast<int8_t>(nNeighbors); i < nNeighbors + 1; i++)
    {
        a = (position[0] - (center[0] - i * lBox));
        a *= a;
        for(int8_t j = - static_cast<int8_t>(nNeighbors); j < nNeighbors + 1; j++)
        {
            b = (position[1] - (center[1] - j * lBox));
            b = (b * b) + a;
            for(int8_t k = - static_cast<int8_t>(nNeighbors); k < nNeighbors + 1; k++)
            {
                c = (position[2] - (center[2] - k * lBox));
                c = (c * c) + b;
                d = sqrt(c) - radius,
                c *= invRange2;
                d = d * (- d * width_factor) - d;
                c = exp(d);
                rValue += c;
            }
        }
    }
    return rValue;
}