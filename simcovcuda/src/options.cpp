#include "options.hpp"

Options parseArgs(int argc, char** argv) {
    Options o;

    std::unordered_map<std::string, int*> argMapInts;
    std::unordered_map<std::string, long long*> argMapLongLongs;
    std::unordered_map<std::string, float*> argMapFloats;
    std::unordered_map<std::string, bool*> argMapBools;
    std::unordered_map<std::string, std::string*> argMapStrings;

    argMapInts.insert(std::make_pair("--blocks", &o.blocks));
    argMapInts.insert(std::make_pair("--threadsPerBlock", &o.threadsPerBlock));
    argMapInts.insert(std::make_pair("--dimX", &o.dimX));
    argMapInts.insert(std::make_pair("--dimY", &o.dimY));
    argMapInts.insert(std::make_pair("--dimZ", &o.dimZ));
    argMapInts.insert(std::make_pair("--tileX", &o.tileX));
    argMapInts.insert(std::make_pair("--tileY", &o.tileY));
    argMapInts.insert(std::make_pair("--tileZ", &o.tileZ));
    argMapInts.insert(std::make_pair("--rankDimX", &o.rankDimX));
    argMapInts.insert(std::make_pair("--rankDimY", &o.rankDimY));
    argMapInts.insert(std::make_pair("--rankDimZ", &o.rankDimZ));
    argMapInts.insert(std::make_pair("--timeSteps", &o.timeSteps));
    argMapInts.insert(std::make_pair("--numInfections", &o.numInfections));
    argMapInts.insert(std::make_pair("--initialInfections", &o.initialInfections));
    argMapInts.insert(std::make_pair("--incubationPeriod", &o.incubationPeriod));
    argMapInts.insert(std::make_pair("--apoptosisPeriod", &o.apoptosisPeriod));
    argMapInts.insert(std::make_pair("--expressingPeriod", &o.expressingPeriod));
    argMapInts.insert(std::make_pair("--tcellGenerationRate", &o.tcellGenerationRate));
    argMapInts.insert(std::make_pair("--tcellInitialDelay", &o.tcellInitialDelay));
    argMapInts.insert(std::make_pair("--tcellVascularPeriod", &o.tcellVascularPeriod));
    argMapInts.insert(std::make_pair("--tcellTissuePeriod", &o.tcellTissuePeriod));
    argMapInts.insert(std::make_pair("--tcellBindingPeriod", &o.tcellBindingPeriod));
    argMapInts.insert(std::make_pair("--samplePeriod", &o.samplePeriod));
    argMapInts.insert(std::make_pair("--reducePeriod", &o.reducePeriod));
    argMapInts.insert(std::make_pair("--checkActivesRate", &o.checkActivesRate));
    argMapInts.insert(std::make_pair("--seed", &o.seed));

    argMapFloats.insert(std::make_pair("--infectivity", &o.infectivity));
    argMapFloats.insert(std::make_pair("--virionProduction", &o.virionProduction));
    argMapFloats.insert(std::make_pair("--virionClearance", &o.virionClearance));
    argMapFloats.insert(std::make_pair("--virionDiffusion", &o.virionDiffusion));
    argMapFloats.insert(std::make_pair("--inflammationProduction", &o.inflammationProduction));
    argMapFloats.insert(std::make_pair("--virionDiffusion", &o.virionDiffusion));
    argMapFloats.insert(std::make_pair("--inflammationProduction", &o.inflammationProduction));
    argMapFloats.insert(std::make_pair("--inflammationDecay", &o.inflammationDecay));
    argMapFloats.insert(std::make_pair("--inflammationDiffusion", &o.inflammationDiffusion));
    argMapFloats.insert(std::make_pair("--minInflammation", &o.minInflammation));
    argMapFloats.insert(std::make_pair("--minVirions", &o.minVirions));
    argMapFloats.insert(std::make_pair("--maxBindingProb", &o.maxBindingProb));

    argMapLongLongs.insert(std::make_pair("--wholeLungX", &o.wholeLungX));
    argMapLongLongs.insert(std::make_pair("--wholeLungY", &o.wholeLungY));
    argMapLongLongs.insert(std::make_pair("--wholeLungZ", &o.wholeLungZ));

    argMapStrings.insert(std::make_pair("--output", &o.output));
    argMapStrings.insert(std::make_pair("--config", &o.config));
    argMapStrings.insert(std::make_pair("--gpuconfig", &o.gpuconfig));

    argMapBools.insert(std::make_pair("--commTest", &o.commTest));

    //set arguments from the command line
    for(int i = 1; i < argc; i++) {
        std::string argument = argv[i];
        std::string indicator = "--";

        if(argument.find(indicator) != std::string::npos) {
            if(argMapInts.find(argument) != argMapInts.end()){
                *(argMapInts.at(argument)) = atoi(argv[i+1]);
            } else if (argMapLongLongs.find(argument) != argMapLongLongs.end()) {
                *(argMapLongLongs.at(argument)) = atoll(argv[i+1]);
            } else if (argMapFloats.find(argument) != argMapFloats.end()) {
                *(argMapFloats.at(argument)) = atof(argv[i+1]);
            } else if (argMapBools.find(argument) != argMapBools.end()) {
                int val = atoi(argv[i+1]);
                if (val == 0) {
                    *(argMapBools.at(argument)) = false;
                } else {
                    *(argMapBools.at(argument)) = true;
                }
            } else if (argMapStrings.find(argument) != argMapStrings.end()){
                *(argMapStrings.at(argument)) = argv[i+1];
            }
        }

        if(argument == "--config"){
            //read config file
            std::ifstream configFile(argv[i+1]);
            std::string parameter;
            std::string value;
            while(configFile >> parameter >> value){
                if(argMapInts.find(parameter) != argMapInts.end()){
                    *(argMapInts.at(parameter)) = atoi(value.c_str());
                } else if (argMapLongLongs.find(parameter) != argMapLongLongs.end()) {
                    *(argMapLongLongs.at(parameter)) = atoll(value.c_str());
                } else if (argMapFloats.find(parameter) != argMapFloats.end()) {
                    *(argMapFloats.at(parameter)) = atof(value.c_str());
                } else if (argMapBools.find(parameter) != argMapBools.end()) {
                    int val = atoi(value.c_str());
                    if (val == 0) {
                        *(argMapBools.at(parameter)) = false;
                    } else {
                        *(argMapBools.at(parameter)) = true;
                    }
                } else if (argMapStrings.find(parameter) != argMapStrings.end()){
                    *(argMapStrings.at(parameter)) = value.c_str();
                }
            }
        }

    }

    return o;
}