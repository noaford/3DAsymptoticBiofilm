//
//  main.cpp
//  Asymptotic3DBiofilm
//
//  Created by Noah Ford on 1/16/18.
//  Copyright © 2018 Noah Ford. All rights reserved.
//

#include <iostream>
#include "Simulator.hpp"

int main(int argc, const char * argv[]) {
    Simulator simulator;
    simulator.ReadInputFile();
    simulator.Initialize();
    simulator.RunSimulation();
    simulator.Output();
    std::cout<<"Finished"<< std::endl;
    
}

