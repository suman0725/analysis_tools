
{
  gInterpreter->AddIncludePath("/usr/local/include");

  gSystem->Load("libhipo4.so");
  gSystem->Load("libROOTDataFrame.so");
  gSystem->Load("libHipoDataFrame.so");

  gInterpreter->ProcessLine("#include \"hipo4/RHipoDS.hxx\"");
}
