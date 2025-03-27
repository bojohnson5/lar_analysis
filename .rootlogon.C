{

  cout << "Loading CENNSMC Root Utilities" << endl;

  //Load libraries needed for CENNSMC analysis                                                       
  const char* includepath = gSystem->GetIncludePath();
  TString temp(includepath);
  temp += "/N/u/bojohn/Quartz/cohar750_sim/analysis/include ";
  gSystem->SetIncludePath(temp.Data());
  //char* foo = ".include ";
  //TString temp2(foo);
  //temp2 += gSystem->Getenv("SCIBATHANLINSTALL");
  //temp2 += "/include";
  //gROOT->ProcessLine(temp2.Data());
  //TString temp3(foo);
  //gROOT->ProcessLine(temp3.Data());

  //ROOT LIBS
  //gSystem->Load("libGeom.so");
  //gSystem->Load("libEG.so"); //to get TDatabasePDG                        
  //gSystem->Load("libRint.so");
  //gSystem->Load("libRGL.so");
  //gSystem->Load("libMathCore.so");
  //gSystem->Load("libHist.so");

  //CENNSMC LIBS
  gSystem->Load("/N/u/bojohn/Quartz/cohar750_sim/builddir/lib/libcenns_io.so");
}
