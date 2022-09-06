/**
 * @file inputoutput.h
 * Header file for class \ref inputoutput
 */

#pragma once

#include <vector>
#include <string>
#include <ostream>
#include <fstream>
#include "yaml-cpp/yaml.h"

class domain;

using namespace std;

////////////////////////////////////////////////////////////////////////////////

/** Class implementing inputoutput object
 *
 *  @author David O. Lignell
 */

class inputoutput {

    public:

    //////////////////// DATA MEMBERS //////////////////////

        domain                   *domn;          ///< pointer to domain object

        ostream                  *ostrm;         ///< Runtime: points to cout or to a file

        string                   caseName;       ///< input file directory
        string                   inputFileDir;   ///< input file directory
        string                   dataDir;        ///< data directory (output)

        YAML::Node               inputFile;      ///< yaml input file object base node
        YAML::Node               params;         ///< yaml sub node
        YAML::Node               streamProps;    ///< yaml sub node
        YAML::Node               initParams;     ///< yaml sub node

        YAML::Node               dvParams;       ///< yaml sub node
        YAML::Node               dTimes;         ///< yaml sub node
        YAML::Node               dumpTimesGen;   ///< yaml sub node
 
        YAML::Node               scalarSc;       ///< yaml sub node: for hips
      //  YAML::Node               pPos;           ///< yaml sub node

      //  vector<double>           probePos;       ///< vector of probe positions
        vector<double>           dumpTimes;      ///< vector of dump times
        int                      iNextDumpTime;  ///< index of next dump time
        bool                     LdoDump;        ///< flag for whether we are dumping a file

        ofstream                 gnufile;        ///< gnuplot script file

        int                      nShift;         ///< for file naming


    //////////////////// MEMBER FUNCTIONS /////////////////

    void outputProperties(const string fname,
                          const double time);    ///< actually write the data file
    void dumpDomainIfNeeded();                   ///< calls outputProperties for dumpTimes
    void writeDataFile(const string fnameRaw,
                       const double time);       ///< writes the gnuplot file and calls outputProperties
    void outputHeader();                         ///< output header info during odt solution
       
     void loadVarsFromRestartFile();

    void set_iNextDumpTime(double time);

        private:


    //////////////////// CONSTRUCTOR FUNCTIONS /////////////////

    public:

        inputoutput(domain *p_domn, const string p_inputFileDir, const int nShift);
        //void init(domain *p_domn);
        ~inputoutput();

};


////////////////////////////////////////////////////////////////////////////////


