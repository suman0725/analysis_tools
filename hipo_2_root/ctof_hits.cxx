#include "TCanvas.h"
#include "ROOT/RVec.hxx"

using RVecS = typename ROOT::VecOps::RVec<short>;
using RVecI = typename ROOT::VecOps::RVec<int>;

void ctof_hits(const char* fname = "ahdc_waveforms.root") {
  ROOT::DisableImplicitMT();

  //auto df = MakeHipoDataFrame("results/alert_data.root");
  ROOT::RDataFrame df("data",fname);
  //auto names =  df.GetColumnNames() ;
  //  Shows what is in the file.
  //for( const auto& n : names) {
  //  std::cout << n << "\n";
  //}
  //df.Describe().Print();

  auto h1 = df.Histo1D( {"", "; CTOF TDC ", 0, 0, 100000}, "CTOF_tdc_TDC");

  TCanvas *c = new TCanvas();
  h1->DrawCopy();
  c->SaveAs("results/ctof_hits.png");
}
