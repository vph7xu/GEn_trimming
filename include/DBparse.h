#ifndef DBparse_H
#define DBparse_H


namespace DBparse {

  std::map<TString, TString> DBFileMap {
    {"He3 Polarization", "He3_pol.csv"},
    {"Beam Polarization", "Beam_pol.csv"},
    {"Helicity Quality", "Helicity_quality.csv"},
    {"Moller Quality", "Moller_quality.csv"},
    {"Field Measurement", "Field_Meas.csv"},
    {"Asymmetry Correction", "corr.csv"}
  };
  
  struct DBrequest{
    TString var_names;   // Variable name
    TString info;        // More description about this variable
    bool    mandatory;   // Is variable mandatory?
  };

  struct DBInfo{
    TString                  cfg;
    vector<DBrequest>        var_req;
    map<TDatime,pair<double,double>>      He3Pol;
    vector<vector<TDatime*>> BeamPolTime;
    vector<vector<double>>   BeamPolValue;
    map<int,double>          GoodHel;
    map<int,double>          GoodMoller;
    double                   AccidentalAsymmetry;
    double                   AccidentalAsymmetryErr;
    double                   AccidentalFraction;
    double                   AccidentalFractionErr;
    double                   PionAsymmetry;
    double                   PionAsymmetryErr;
    double                   PionFraction;
    double                   PionFractionErr;
    double                   InelasticAsymmetry;
    double                   InelasticAsymmetryErr;
    double                   InelasticFraction;
    double                   InelasticFractionErr;
    double                   NitrogenFraction;
    double                   NitrogenFractionErr;
  };


  void DB_load(DBInfo &request);

 
}

#endif
