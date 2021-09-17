#include <cstdint>
#include <cstring>
#include <fstream>
#include <iostream>
#include <chrono>
#include "omp.h"
#include "simulation.hpp"

/*
./exec --nthread 8 --nstep 100 --nparticles 50 --sbox 1.0 --dx 1e-4 --dt 1e-2 --lightspeed 0.5 --width 0.5 --range 5e-2
*/
auto main(int argc, char* argv[]) -> int
{
    float lBox{1.f};
    float dX{1e-4f};
    float dT{1e-2f};
    float wave_speed{1.f};
    float wave_witdh{1e-1f};
    float wave_range{5e-2f};
    uint16_t nStep{100u};
    uint16_t nParticles{50u};
    uint8_t nDim{2u};
    uint8_t nThreads{4u};
    bool walls{false};

    for(auto i{1}; i < argc; i++)
    {
		if(argv[i] == "-h" || argv[i] == "--help")
        {
			std::cout << "Usage: App <options>\nOptions are: \n";
			std::cout << "--nthread: define the number of thread used during parallel computation.\n";
			std::cout << "--nstep: define the number of step in the simulation.\n";
			std::cout << "--nparticles: define the number of particles in the simulation.\n";
			std::cout << "--lbox: define the size of the simulated box.\n";
			std::cout << "--dx: define the differential spatial element of the simulation.\n";
			std::cout << "--dt: define the time step.\n";
			std::cout << "--speed: define the wave's speed.\n";
			std::cout << "--witdh: define the wave's width.\n";
			std::cout << "--range: define the wave's range.\n";
			std::cout << "--walls: Will enabled the bouncing at box limit otherwise the box will be cyclic.\n";
            return EXIT_SUCCESS;
        }else if (argv[i] == "--nthread")
        {
            i++;
            nThreads = std::stoi(argv[i]);
        }else if (argv[i] == "--nstep")
        {
            i++;
            nStep = std::stoi(argv[i]);
        }else if (argv[i] == "--nparticles")
        {
            i++;
            nParticles = std::stoi(argv[i]);
        }else if (argv[i] == "--lbox")
        {
            i++;
            lBox = std::stof(argv[i]);
        }else if (argv[i] == "--dx")
        {
            i++;
            dX = std::stof(argv[i]);
        }else if (argv[i] == "--dt")
        {
            i++;
            dT = std::stof(argv[i]);
        }else if (argv[i] == "--speed")
        {
            i++;
            wave_speed = std::stof(argv[i]);
        }else if (argv[i] == "--witdh")
        {
            i++;
            wave_witdh = std::stof(argv[i]);
        }else if (argv[i] == "--range")
        {
            i++;
            wave_range = std::stof(argv[i]);
        }else if (argv[i] == "--walls")
        {
            walls = true;
        }
    }

    std::cout << "Parameters review..." << std::endl;
    std::cout << "Threads amount: " << static_cast<int>(nThreads) << std::endl;
    std::cout << "Steps amount: " << nStep << std::endl;
    std::cout << "Particles number: " << nParticles << std::endl;
    std::cout << "Box size: " << lBox << std::endl;
    std::cout << "Differential space unit: " << dX << std::endl;
    std::cout << "Differential time unit: " << dT << std::endl;
    std::cout << "Interacting wave speed: " << wave_speed << std::endl;
    std::cout << "Interacting wave witdh: " << wave_witdh << std::endl;
    std::cout << "Interacting wave range: " << wave_range << std::endl;
    std::cout << "Dimension: " << static_cast<int>(nDim) << std::endl;
    std::cout << "Walls: " << std::boolalpha << walls << std::endl;
    std::cout << "Enter \"ok\" to continue or write something else to exit..." << std::endl;

    std::string paramsValidation{""};
    std::cin >> paramsValidation;

    if(paramsValidation != "ok")
        return EXIT_SUCCESS;

    omp_set_dynamic(false);
    omp_set_num_threads(static_cast<int>(nThreads));
    auto start = std::chrono::high_resolution_clock::now();
    simulate(lBox, dX, dT, wave_speed, wave_witdh, wave_range, nStep, nParticles, nDim, walls);
    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float> elapsed{ finish - start };
    std::cout << "Elapsed time: " << elapsed.count() << " s\n";
    std::ofstream plot_output{"plot.sh"};
    plot_output << "python plot.py -nstep " << nStep << " -lbox " << lBox;
    plot_output.close();

    
    return EXIT_SUCCESS;
}