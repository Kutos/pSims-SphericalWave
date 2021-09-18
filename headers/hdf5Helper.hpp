#ifndef __H5__
#define __H5__
#include "particle.hpp"
#include "H5Cpp.h"
using namespace H5;
class hdf5Helper
{
private:
    H5File file;
    DataSet simDataset;
    hsize_t simDatasetSize[1];

    struct sParams
    {
    const float lBox;
    const float dX;
    const float dT;
    const float wave_speed;
    const float wave_width;
    const float wave_range;
    const int nStep;
    const int nParticles;
    const int nDim;
    const bool walls;
    };

    hid_t sParticle_id;
    DataSpace simMSpace;
    struct sSim1D
    {
        int t;
        int p;
        int x;
    };
    struct sSim2D
    {
        int t;
        int p;
        int x;
        int y;
    };
    struct sSim3D
    {
        int t;
        int p;
        int x;
        int y;
        int z;
    };

public:
    hdf5Helper(
    const char* fileName,
    const float& lBox,
    const float& dX,
    const float& dT,
    const float& wave_speed,
    const float& wave_width,
    const float wave_range,
    const int& nStep,
    const int& nParticles,
    const int& nDim,
    const bool& walls);
    ~hdf5Helper();

    auto writeStep(
        particle* particles_array,
        const int& nParticles,
        const int& step,
        const int& nDim) -> void;


};
#endif