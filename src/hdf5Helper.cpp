#include <iostream>
#include <string>
#include "H5Cpp.h"
#include "particle.hpp"
#include "hdf5Helper.hpp"
using namespace H5;
hdf5Helper::hdf5Helper(
    const char* flname,
    const float& lBox,
    const float& dX,
    const float& dT,
    const float& wave_speed,
    const float& wave_width,
    const float wave_range,
    const int& nStep,
    const int& nParticles,
    const int& nDim,
    const bool& walls)
{
    try
    {
        file = H5Fcreate(flname, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
        hid_t sParams_id = H5Tcreate(H5T_COMPOUND, sizeof(sParams));
        H5Tinsert(sParams_id, "lBox", HOFFSET(sParams, lBox), H5T_NATIVE_FLOAT);
        H5Tinsert(sParams_id, "dx", HOFFSET(sParams, dX), H5T_NATIVE_FLOAT);
        H5Tinsert(sParams_id, "dt", HOFFSET(sParams, dT), H5T_NATIVE_FLOAT);
        H5Tinsert(sParams_id, "wave_speed", HOFFSET(sParams, wave_speed), H5T_NATIVE_FLOAT);
        H5Tinsert(sParams_id, "wave_witdh", HOFFSET(sParams, wave_width), H5T_NATIVE_FLOAT);
        H5Tinsert(sParams_id, "wave_range", HOFFSET(sParams, wave_range), H5T_NATIVE_FLOAT);
        H5Tinsert(sParams_id, "nStep", HOFFSET(sParams, nStep), H5T_NATIVE_INT);
        H5Tinsert(sParams_id, "nParticles", HOFFSET(sParams, nParticles), H5T_NATIVE_INT);
        H5Tinsert(sParams_id, "nDim", HOFFSET(sParams, nDim), H5T_NATIVE_INT);
        H5Tinsert(sParams_id, "walls", HOFFSET(sParams, walls), H5T_NATIVE_HBOOL);

        hsize_t paramsDims[1] = {1};
        hsize_t paramsMaxDims[1] = {1};
        DataSpace paramsMSpace(1, paramsDims, paramsMaxDims);
        DataSet paramsDataset = file.createDataSet("/params", sParams_id, paramsMSpace, H5P_DEFAULT);
        DataSpace pSpace = paramsDataset.getSpace();
        sParams data[1]{ {lBox, dX, dT, wave_speed, wave_width, wave_range, nStep, nParticles, nDim, walls} };
        paramsDataset.write(data, sParams_id, paramsMSpace, pSpace);
        paramsMSpace.close();
        pSpace.close();
        paramsDataset.close();


        hsize_t simDims[1] = {nParticles};
        hsize_t simMDims[1] = {H5S_UNLIMITED};
        DataSpace simMSpace(1, simDims, simMDims);
        sParticle_id = H5Tcreate(H5T_COMPOUND, 8 + 4*nDim);
        H5Tinsert(sParticle_id, "Time", HOFFSET(sSim, t), H5T_NATIVE_INT);
        H5Tinsert(sParticle_id, "pID", HOFFSET(sSim, p), H5T_NATIVE_INT);
        H5Tinsert(sParticle_id, "x", HOFFSET(sSim, x), H5T_NATIVE_FLOAT);
        if(nDim > 1)
        {
            H5Tinsert(sParticle_id, "y", HOFFSET(sSim, y), H5T_NATIVE_FLOAT);
            /*if(nDim > 2)
            {
                H5Tinsert(sParticle_id, "z", HOFFSET(sSim, z), H5T_NATIVE_FLOAT);
            } */
        }
        DSetCreatPropList cparms;

        hsize_t chunk_dims[1] = {nParticles};
        cparms.setChunk(1, chunk_dims);
        simDataset = file.createDataSet("/simulation", sParticle_id, simMSpace, cparms);
        simDatasetSize[0] = 0;
    }
    // catch failure caused by the H5File operations
    catch (FileIException error) {
        throw error;
    }

    // catch failure caused by the DataSet operations
    catch (DataSetIException error) {
        throw error;
    }

    // catch failure caused by the DataSpace operations
    catch (DataSpaceIException error) {
        throw error;
    }

    // catch failure caused by the DataSpace operations
    catch (DataTypeIException error) {
        throw error;
    }
}

hdf5Helper::~hdf5Helper(){}


auto hdf5Helper::writeStep(
        particle* particles_array,
        const int& nParticles,
        const int& step,
        const int& nDim) -> void
{
    simDatasetSize[0] += nParticles;
    simDataset.extend(simDatasetSize);
    DataSpace fspace = simDataset.getSpace();
    hsize_t offset[1] = {nParticles * step};
    hsize_t dims[1] = {nParticles};
    fspace.selectHyperslab(H5S_SELECT_SET, dims, offset);
    DataSpace mspace(1, dims);

    sSim data[nParticles];
    for(uint i{0}; i < nParticles; i++)
    {
        data[i].t = step;
        data[i].p = i;
        data[i].x = particles_array[i].getPosition(0);
        if(nDim > 1) data[i].y = particles_array[i].getPosition(1);
        //if(nDim > 2) data[i].z = particles_array[i].getPosition(2);
    }
    simDataset.write(data, sParticle_id, mspace, fspace);
}