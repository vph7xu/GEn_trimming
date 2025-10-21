
#include "../../include/gen-ana.h"

map<int,double> fHe3Pol;
map<int,double> fBeamPol;
map<int,double> fGoodHel;

void getDB(TString DB_file){

  // Entries should follow this form:
  //{var name, variable pointer, vairable description, 1/0 (mandatory/not mandatory variable)}
  DBparse::DBRequest request[] = {
    {"He3 Polarization", &fHe3Pol, "He3 Polarization", 0},
    {"Beam Polarization", &fBeamPol, "Beam Polarization", 0},
    {"Good Helicity", &fGoodHel, "Is the helicity readback good (0/1 = no/yes)", 1}
  };
  
  const int nvar = sizeof(request) / sizeof(request[0]);
  
  DB_load(DB_file,request,nvar);
}


void DB_test(TString conf = "GEN2"){

  //TString DB_file = "../../DB/GEN2.csv";
  TString DB_file = "../../DB/Helicity_quality.csv";
 
  getDB(DB_file);

  cout<<fGoodHel[2034]<<endl;

}
