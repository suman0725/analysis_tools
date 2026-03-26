#include "TCanvas.h"
#include "ROOT/RVec.hxx"

using RVecS = typename ROOT::VecOps::RVec<short>;
using RVecI = typename ROOT::VecOps::RVec<int>;

void atof_tdcs(const char* fname = "alert_test.root") {
  ROOT::DisableImplicitMT();

  //auto df = MakeHipoDataFrame("results/alert_data.root");
  ROOT::RDataFrame df("data",fname);
  auto names =  df.GetColumnNames() ;
  //  Shows what is in the file.
  for( const auto& n : names) {
    std::cout << n << "\n";
  }
  df.Describe().Print();
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
               .Alias("tdc", "ATOF_tdc_TDC")
               .Alias("tot", "ATOF_tdc_ToT");


  // there is probably a better way to do this. I think shorts make it annoying.
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
  //auto is_bar = [](const RVecI& c) {
  //  return  c[c==10];
  //};
  auto t_diffs = [](const std::vector<std::tuple<int,int>> hits,const RVecI &tdc){
    RVecI res;
    for(const auto& [i,j] : hits) {
      res.push_back(tdc[i] -tdc[j]);
    }
    return res;
  };
  auto dual_mod = [](const std::vector<std::tuple<int,int>> hits,const RVecI &s){
    RVecI res;
    for (const auto &[i, j] : hits) {
      res.push_back(s[i]);
    }
    return res;
  };

  auto t_sums = [](const std::vector<std::tuple<int,int>> hits,const RVecI &tdc){
    RVecI res;
    for(const auto& [i,j] : hits) {
      res.push_back(tdc[i] +tdc[j]);
    }
    return res;
  };

  // Loop of hits in same bar
  auto dual_hits =
      [](const RVecI &s, const RVecI &l,
         const RVecI &c, const RVecI &o) {
        std::vector<std::tuple<int,int>> res;
        int i = 0;
        for (const auto &si : s) {
          int j = 0;
          if (o[i] == 0) {
            for (const auto &sj : s) {
              if ((si == sj) && (l[i] == l[j])) {
                if ( o[j] == 1) {
                  if ((c[i] == 10) && (c[j] == 10)) {
                    res.push_back(std::make_tuple(i,j));
                  }
                }
              }
              j++;
            }
          }
          i++;
        }
        return res;
      };

  auto df2 = df1.Filter("ATOF_tdc_layer.size() > 0")
                //.Define("is_Bar",is_bar,{"c"})
                //.Define("is_Wedge","!is_Bar")
                .Define("ATOF_w", calc_wedge, {"ATOF_s", "ATOF_l"})
                .Define("dual_hits",dual_hits,{"ATOF_s","ATOF_l","ATOF_c","ATOF_o"})
                .Define("dual_module", dual_mod,{"dual_hits","ATOF_s"})
                .Define("dual_layer", dual_mod,{"dual_hits","ATOF_l"})
                .Define("t_diff", t_diffs,{"dual_hits","ATOF_tdc_TDC"})
                .Define("t_sum", t_sums,{"dual_hits","ATOF_tdc_TDC"})
                .Define("t_diff_0","t_diff[dual_layer == 0]")
                .Define("t_diff_1","t_diff[dual_layer == 1]")
                .Define("t_diff_2","t_diff[dual_layer == 2]")
                .Define("t_diff_3","t_diff[dual_layer == 3]")
                .Define("t_sum_0","t_sum[dual_layer == 0]")
                .Define("t_sum_1","t_sum[dual_layer == 1]")
                .Define("t_sum_2","t_sum[dual_layer == 2]")
                .Define("t_sum_3","t_sum[dual_layer == 3]")
                .Filter("dual_hits.size()>0")
                ;

  auto disp = df2.Display({"ATOF_w","ATOF_tdc_TDC"},50,100);
  //.Display({"s","l","c","o","tdc","tot"},50,48);
  //.Display({"ATOF_tdc_sector", "ATOF_tdc_layer", "ATOF_tdc_component",
  //          "ATOF_tdc_order", "ATOF_tdc_TDC", "ATOF_tdc_ToT"},
  //         50, 48);
  //

  auto h1 = df2.Histo1D({"Module", "Module;Module Number", 15, 0, 15}, "ATOF_s");
  auto h1_wedge =
      df2.Histo1D({"Wedge", "; global wedge number", 60, 0, 60}, "ATOF_w");

  auto h1_diff = df2.Histo1D( {"", "; TDC diff", 200, -1000, 1000}, "t_diff");
  auto h1_diff0 = df2.Histo1D({"", "; TDC diff", 200, -1000, 1000}, "t_diff_0");
  auto h1_diff1 = df2.Histo1D({"", "; TDC diff", 200, -1000, 1000}, "t_diff_1");
  auto h1_diff2 = df2.Histo1D({"", "; TDC diff", 200, -1000, 1000}, "t_diff_2");
  auto h1_diff3 = df2.Histo1D({"", "; TDC diff", 200, -1000, 1000}, "t_diff_3");

  auto h1_sum = df2.Histo1D( {"", "; TDC sum", 200, 20e3  , 30e3}, "t_sum");
  auto h1_sum0 = df2.Histo1D({"", "; TDC sum", 200, 20e3  , 30e3}, "t_sum_0");
  auto h1_sum1 = df2.Histo1D({"", "; TDC sum", 200, 20e3  , 30e3}, "t_sum_1");
  auto h1_sum2 = df2.Histo1D({"", "; TDC sum", 200, 20e3  , 30e3}, "t_sum_2");
  auto h1_sum3 = df2.Histo1D({"", "; TDC sum", 200, 20e3  , 30e3}, "t_sum_3");

  TCanvas *c = new TCanvas();

  disp->Print();
  h1->DrawCopy();
  c->SaveAs("results/atof_hits.png");

  h1_wedge->DrawCopy();
  c->SaveAs("results/atof_wedge.png");

  TCanvas *c2 = new TCanvas("c2","c2",1200,1000);
  c2->Divide(1,2);

  c2->cd(1);
  TH1* h = h1_diff->DrawCopy();
  h->SetLineColor(1);
  h = h1_diff0->DrawCopy("same");
  h->SetLineColor(2);
  h = h1_diff1->DrawCopy("same");
  h->SetLineColor(4);
  h = h1_diff2->DrawCopy("same");
  h->SetLineColor(7);
  h = h1_diff3->DrawCopy("same");
  h->SetLineColor(8);

  c2->cd(2);
  h = h1_sum->DrawCopy();
  h->SetLineColor(1);
  h = h1_sum0->DrawCopy("same");
  h->SetLineColor(2);
  h = h1_sum1->DrawCopy("same");
  h->SetLineColor(4);
  h = h1_sum2->DrawCopy("same");
  h->SetLineColor(7);
  h = h1_sum3->DrawCopy("same");
  h->SetLineColor(8);


  c2->SaveAs("results/atof_diffsum.png");
}
