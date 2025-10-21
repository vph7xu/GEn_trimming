#include "../include/DBparse.h"
#include "../include/Analysis.h"

namespace DBparse {


  void DB_load(DBInfo &request){

    TString DB_dir = "/work/halla/sbs/vimukthi/GEn_analysis/DB/"; 

    for(int ivar = 0; ivar < request.var_req.size(); ivar++){
      DBrequest var_info = request.var_req[ivar];
      bool found_var = false;

      auto it = DBFileMap.find(var_info.var_names);
      TString file = it->second;
      if(it->first == "Asymmetry Correction")
	file = request.cfg + "_" + file;
      
      fstream file_csv; file_csv.open(DB_dir + file);
      
      // Here we check that the DB files are there
      if (it == DBFileMap.end()) { // This is not in the list of DB files
        if(!var_info.mandatory){
	  cout<<"DB variable NOT found but NOT mandatory: "<<var_info.var_names<<endl;
	  continue;
	}
	if(var_info.mandatory){
	  cout<<"!!!!WARNING, ERROR, DBparse::DB_load: DB variable NOT found and IS mandatory: "<<var_info.var_names<<endl;
	  exit(0);
	}
      }
      else if (!file_csv.is_open()) {// This is in the list but the file is not there
	if(!var_info.mandatory){
	  cout<<"DB file NOT found but NOT mandatory: "<<file<<endl;
	  continue;
	}
	if(var_info.mandatory){
	  cout<<"!!!!WARNING, ERROR, DBparse::DB_load: File "<<file<<" is not found"<<endl;
	  exit(0);
	}
      }
      else {// The file is found and everything is fine
	cout<<"DB file found: "<<file<<endl;
      }
      
      string line;
      int iline = 0;
      
      while (getline(file_csv, line)) {
	// Create a stringstream to parse the line
	stringstream ss(line);
	string cell;
	iline++;
	
	vector<string> val;

	// Split the line into cells using a comma as a delimiter
	while (getline(ss, cell, ',')) {
	  val.push_back(cell);  // Put one line into vectros
	}
	if(it->first == "He3 Polarization"){// He3 polarization DB
	  if(iline == 1) continue;
	  
	  if (val.size() < 3) {
            std::cerr << "DBparse warning: Skipping malformed line ► "
              	      << line << "\n";
            continue;                   // <-- avoid seg-fault here
	  }

	  try {
            double pol  = std::stod(val[1]);
            double err  = std::stod(val[2]);
            //time_t tkey = Utilities::SetTime(val[0]);

            // map<time_t, pair<double,double>>
            request.He3Pol.emplace(Utilities::SetTime(val[0]), std::make_pair(pol, err));
          }
          catch (const std::exception& e) {
            std::cerr << "DBparse error: " << e.what()
                      << " while parsing line ► " << line << '\n';
          }
	  
	  //request.He3Pol.emplace((Utilities::SetTime(val[0])),make_pair(stod(val[1]),stod(val[2])));
	
	}// end He3 Pol
	else if(it->first == "Beam Polarization"){// Beam polarization DB
	  if(iline == 1) continue;	 
	  vector<TDatime*> beam_time;
	  vector<double> beam_val;
	  
	  for(int irow = 0; irow < val.size(); irow++){ 
	    if(irow < 2) { // Time numbers
	      TDatime *time = new TDatime(val[irow].c_str());
	      beam_time.push_back(time);
	    }
	    if(irow >= 2) { // Values
	      beam_val.push_back(stod(val[irow]));
	    }
	  }
	  
	  request.BeamPolTime.push_back(beam_time);
	  request.BeamPolValue.push_back(beam_val); 
	} // End Beam Pol
	else if(it->first == "Helicity Quality"){// Helicity quality DB
	  if(iline == 1) continue;	  
	  request.GoodHel.insert(std::make_pair(stoi(val[0]),stod(val[1])));
	}
	else if(it->first == "Moller Quality"){// Moller quality DB
	  if(iline == 1) continue;	 
	  request.GoodMoller.insert(std::make_pair(stoi(val[0]),stod(val[1])));
	}
	else if(it->first == "Field Measurement"){// Field Measurement DB
	  if(iline == 1) continue;	  
	  if(request.cfg == val[0]){
	    field_hor_ang = stod(val[1]);
	    field_vert_ang = stod(val[2]);
	  }
	}// End Field Measurement
	else if(it->first == "Field Measurement"){// Field Measurement DB
	  if(iline == 1) continue;	  
	  if(request.cfg == val[0]){
	    field_hor_ang = stod(val[1]);
	    field_vert_ang = stod(val[2]);
	  }
	}// End Field Measurement
	else if(it->first == "Asymmetry Correction"){// Asym Correction DB
	  if(val.size() == 0) continue;

	  if(val[0] == "pion"){// Pion Corrections
	    request.PionAsymmetry = stod(val[1]);
	    request.PionAsymmetryErr = stod(val[2]);
	    request.PionFraction = stod(val[3]);
	    request.PionFractionErr = stod(val[4]);
	  } 
	  else if(val[0] == "Inelastic"){// Inelastic Corrections
	    request.InelasticAsymmetry = stod(val[1]);
	    request.InelasticAsymmetryErr = stod(val[2]);
	    request.InelasticFraction = stod(val[3]);
	    request.InelasticFractionErr = stod(val[4]);
	  } 
	  else if(val[0] == "Nitrogen"){// Nitrogen Corrections
	    request.NitrogenFraction = stod(val[1]);
	    request.NitrogenFractionErr = stod(val[2]);
	  } 
	  else if(val[0] == "time accidentals"){// Accidental Corrections
	    request.AccidentalAsymmetry = stod(val[1]);
	    request.AccidentalAsymmetryErr = stod(val[2]);
	    request.AccidentalFraction = stod(val[3]);
	    request.AccidentalFractionErr = stod(val[4]);
	  } 

	}// End Asym Correction 
      }
      
    }
  }

}
