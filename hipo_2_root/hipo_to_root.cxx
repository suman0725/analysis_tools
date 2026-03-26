#include "TCanvas.h"
#include "hipo4/RHipoDS.hxx"
#include "ROOT/RVec.hxx"

using RVecS = typename ROOT::VecOps::RVec<short>;

void hipo_to_root(const char* fname = "alert_test.hipo") {
  ROOT::DisableImplicitMT();

  auto df = MakeHipoDataFrame(fname);
  auto names =  df.GetColumnNames() ;
  //  Shows what is in the file.
  for( const auto& n : names) {
    std::cout << n << "\n";
  }
  // df.Describe().Print();
  // auto c1 = df.Count();
  // std::cout << *c1 << " events\n";
  auto convert_to_int = [](const std::vector<short> &s) {
    RVecI res;
    for(const auto& v:s) {
      res.push_back(int(v));
    }
    return res;
  };
  auto df1 = df.Define("ATOF_s", convert_to_int, {"ATOF_tdc_sector"})
               .Define("ATOF_l", convert_to_int, {"ATOF_tdc_layer"})
               .Define("ATOF_c", convert_to_int, {"ATOF_tdc_component"})
               .Define("ATOF_o", convert_to_int, {"ATOF_tdc_order"})
               ;

  //
  auto calc_wedge = [](const RVecI &s,
                       const RVecI &l) {
    RVecI res;
    int i = 0;
    for (const auto &si : s) {
      res.push_back(int(si) *4 + int(l[i]));
      i++;
    }
    return res;
  };

  auto df2 = df1.Filter("ATOF_tdc_layer.size() > 0")
                .Define("ATOF_w", calc_wedge, {"ATOF_s", "ATOF_l"})
                ;

  auto disp = df2.Display({"ATOF_w","ATOF_tdc_TDC"},100,100);

  df2.Snapshot("data","results/alert_data.root");

}
