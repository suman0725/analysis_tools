#include "TCanvas.h"
#include "TMultiGraph.h"
#include "ROOT/RVec.hxx"
#include <string>
#include <algorithm>
#include <TInterpreter.h>
#include "TCanvas.h"
//#include "hipo4/RHipoDS.hxx"
#include "ROOT/RVec.hxx"

using RVecS = typename ROOT::VecOps::RVec<short>;
using RVecI = typename ROOT::VecOps::RVec<int>;

struct AHDCWaveform {
  int index;
  int samples[30];
};

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
//auto t_diffs = [](const std::vector<std::tuple<int,int>> hits,const RVecI &tdc){
//  RVecI res;
//  for(const auto& [i,j] : hits) {
//    res.push_back(tdc[i] -tdc[j]);
//  }
//  return res;
//};

void ahdc_waveforms(const char* fname = "alert_test.root") {

  std::cout << "derp\n";

  gInterpreter->GenerateDictionary("vector<ROOT::VecOps::RVec<short> >","vector;ROOT/RVec.hxx");
  gInterpreter->GenerateDictionary("vector<tuple<int,int> >","vector;tuple");

  std::vector<string> wf_cols;

  for(int i = 1; i<=30; i++ ) {
    wf_cols.push_back("AHDC_wf_s"+std::to_string(i));
  }

  ROOT::RDataFrame df("data",fname);
  auto names =  df.GetColumnNames() ;
  //  Shows what is in the file.
  for( const auto& n : names) {
    std::cout << n << "\n";
  }
  df.Describe().Print();
   auto c1 = df.Count();
   std::cout << *c1 << " events\n";

  auto wf_to_array = [](RVecS s1, RVecS s2, RVecS s3,RVecS s4, RVecS s5,RVecS s6, RVecS s7, RVecS s8,RVecS s9, RVecS s10,
                        RVecS s11, RVecS s12, RVecS s13,RVecS s14, RVecS s15,RVecS s16, RVecS s17, RVecS s18,RVecS s19, RVecS s20,
                        RVecS s21, RVecS s22, RVecS s23,RVecS s24, RVecS s25,RVecS s26, RVecS s27, RVecS s28,RVecS s29, RVecS s30
                        ) {
    std::vector<RVecS> res;
    int i = 0;
    for(short v1 : s1) {
      RVecS a1 = {s1[i],s2[i],s3[i],s4[i],s5[i],s6[i],s7[i],s8[i],s9[i],s10[i],
                 s11[i],s12[i],s13[i],s14[i],s15[i],s16[i],s17[i],s18[i],s19[i],s20[i],
                 s21[i],s22[i],s23[i],s24[i],s25[i],s26[i],s27[i],s28[i],s29[i],s30[i]
      };
      res.push_back(a1);
    }
    return res;
  };

  auto convert_to_int = [](const std::vector<short> &s) {
    RVecI res;
    for(const auto& v:s) {
      res.push_back(int(v));
    }
    return res;
  };
  auto df1 = df.Define("ATOF_s", convert_to_int, {"ATOF_tdc_sector"})
  .Define("n_samp",[]()-> RVecS {return {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};})
  .Define("AHDC_arrays",wf_to_array, wf_cols)
  .Define("ATOF_l", convert_to_int, {"ATOF_tdc_layer"})
  .Define("ATOF_c", convert_to_int, {"ATOF_tdc_component"})
  .Define("ATOF_o", convert_to_int, {"ATOF_tdc_order"})
  ;

  // there is probably a better way to do this. I think shorts make it annoying.
  //
  //auto dual_mod = [](const std::vector<std::tuple<int,int>> hits,const RVecI &s){
  //  RVecI res;
  //  for (const auto &[i, j] : hits) {
  //    res.push_back(s[i]);
  //  }
  //  return res;
  //};

  //auto t_sums = [](const std::vector<std::tuple<int,int>> hits,const RVecI &tdc){
  //  RVecI res;
  //  for(const auto& [i,j] : hits) {
  //    res.push_back(tdc[i] +tdc[j]);
  //  }
  //  return res;
  //};

  //// Loop of hits in same bar
  //auto dual_hits =
  //    [](const RVecI &s, const RVecI &l,
  //       const RVecI &c, const RVecI &o) {
  //      std::vector<std::tuple<int,int>> res;
  //      int i = 0;
  //      for (const auto &si : s) {
  //        int j = 0;
  //        if (o[i] == 0) {
  //          for (const auto &sj : s) {
  //            if ((si == sj) && (l[i] == l[j])) {
  //              if ( o[j] == 1) {
  //                if ((c[i] == 10) && (c[j] == 10)) {
  //                  res.push_back(std::make_tuple(i,j));
  //                }
  //              }
  //            }
  //            j++;
  //          }
  //        }
  //        i++;
  //      }
  //      return res;
  //    };

  auto df2 =
      df1.Filter("ATOF_tdc_layer.size() > 0")
          .Filter("AHDC_wf_s1.size() > 0")
          .Define("first_pulse", "AHDC_arrays[0]");
      //              //.Define("is_Bar",is_bar,{"c"})
      //              //.Define("is_Wedge","!is_Bar")
      //              .Define("ATOF_w", calc_wedge, {"ATOF_s", "ATOF_l"})
      //              .Define("dual_hits",dual_hits,{"ATOF_s","ATOF_l","ATOF_c","ATOF_o"})
      //              .Define("dual_module", dual_mod,{"dual_hits","ATOF_s"})
      //              .Define("dual_layer", dual_mod,{"dual_hits","ATOF_l"})
      //              .Define("t_diff", t_diffs,{"dual_hits","ATOF_tdc_TDC"})
      //              .Define("t_sum", t_sums,{"dual_hits","ATOF_tdc_TDC"})
      //              .Define("t_diff_0","t_diff[dual_layer == 0]")
      //              .Define("t_diff_1","t_diff[dual_layer == 1]")
      //              .Define("t_diff_2","t_diff[dual_layer == 2]")
      //              .Define("t_diff_3","t_diff[dual_layer == 3]")
      //              .Define("t_sum_0","t_sum[dual_layer == 0]")
      //              .Define("t_sum_1","t_sum[dual_layer == 1]")
      //              .Define("t_sum_2","t_sum[dual_layer == 2]")
      //              .Define("t_sum_3","t_sum[dual_layer == 3]")
      //              .Filter("dual_hits.size()>0")
      //              ;

      // auto disp = df2.Display({"ATOF_w","ATOF_tdc_TDC"},50,100);
      //.Display({"s","l","c","o","tdc","tot"},50,48);
      //.Display({"ATOF_tdc_sector", "ATOF_tdc_layer", "ATOF_tdc_component",
      //           "ATOF_tdc_order", "ATOF_tdc_TDC", "ATOF_tdc_ToT"},
      //          50, 48);
      //
      auto h1_ahdc_integral = df2.Histo1D(
          {"", "; pulse integral", 200, 0, 50000}, "AHDC_adc_integral");
  auto h1_ahdc_ToT       = df2.Histo1D({"", "; ToT", 200, -10000, 10000}, "AHDC_adc_timeOverThreshold");
  auto h1_ahdc_ToT1       = df2.Define("derp1","AHDC_adc_timeOverThreshold[AHDC_wf_layer>9 && AHDC_wf_layer<20]").Histo1D({"", "; ToT", 200, -10000, 10000}, "derp1");
  auto h1_ahdc_ToT2       = df2.Define("derp2","AHDC_adc_timeOverThreshold[AHDC_wf_layer>19 && AHDC_wf_layer<30]").Histo1D({"", "; ToT", 200, -10000, 10000}, "derp2");
  auto h1_ahdc_ToT3       = df2.Define("derp3","AHDC_adc_timeOverThreshold[AHDC_wf_layer>29 && AHDC_wf_layer<40]").Histo1D({"", "; ToT", 200, -10000, 10000}, "derp3");
  auto h1_ahdc_ToT4       = df2.Define("derp4","AHDC_adc_timeOverThreshold[AHDC_wf_layer>39 && AHDC_wf_layer<50]").Histo1D({"", "; ToT", 200, -10000, 10000}, "derp4");
  auto h1_ahdc_ToT5       = df2.Define("derp5","AHDC_adc_timeOverThreshold[AHDC_wf_layer>49 && AHDC_wf_layer<60]").Histo1D({"", "; ToT", 200, -10000, 10000}, "derp5");
  auto h1_ahdc_cf_time   = df2.Histo1D({"", "; Const. Fraction Time", 200, -20000, 20000}, "AHDC_adc_constantFractionTime");
  auto h1_ahdc_edge_time = df2.Histo1D({"", "; leading edge time", 200, -10000, 10000}, "AHDC_adc_leadingEdgeTime");

  auto h2_ahdc_integral_vs_Tot  = df2.Histo2D({"", "; pulse integral; ToT", 200, 0, 50000,200, 0, 10000}, "AHDC_adc_integral","AHDC_adc_timeOverThreshold");
  auto h2_ahdc_cf_vs_edge  = df2.Histo2D({"", "; Const. Fraction Time; edge time", 200, -20000, 20000,200, -10000, 10000}, "AHDC_adc_constantFractionTime","AHDC_adc_leadingEdgeTime");

  TCanvas *c = new TCanvas();
  TH1 *h = 0;
  // This tutorial is ran with multithreading enabled. The order in which points are inserted is not known, so to have a meaningful representation points are sorted.
  //graph->Sort();
  auto N  = df2.Count();
  int nevents = *N;

  int start_event = 0;
  TMultiGraph * mg = new TMultiGraph();
  for(int j = 0; j<std::min(10,nevents/10); j++){
    for(int i = 0; i<10;i++) {
      auto graph = df2.Range(start_event+i+j*10,start_event+j*10+i+1).Graph("n_samp", "first_pulse");
      auto gr = (TGraph*)graph->DrawClone("APL");
      gr->SetLineColor(i+1); mg->Add(gr);
    }

    ////disp->Print();
    //h1->DrawCopy();
  }
  mg->Draw("al");
  c->SaveAs("results/ahdc_waveforms.png");

  ////h1_wedge->DrawCopy();
  ////c->SaveAs("results/atof_wedge.png");

  //TCanvas *c2 = new TCanvas("c2","c2",1200,1000);
  //c2->Divide(2,2);
  //c2->cd(1);
  //gPad->SetLogy(true);
  //h = h1_ahdc_integral->DrawCopy();
  //c2->cd(2);
  //gPad->SetLogy(true);
  //h = h1_ahdc_ToT->DrawCopy();
  //h = h1_ahdc_ToT1->DrawCopy("same");
  //h->SetLineColor(2);
  //h = h1_ahdc_ToT1->DrawCopy("same");
  //h->SetLineColor(4);
  //h = h1_ahdc_ToT2->DrawCopy("same");
  //h->SetLineColor(6);
  //h = h1_ahdc_ToT3->DrawCopy("same");
  //h->SetLineColor(7);
  //h = h1_ahdc_ToT4->DrawCopy("same");
  //h->SetLineColor(8);
  //h = h1_ahdc_ToT5->DrawCopy("same");
  //h->SetLineColor(9);
  //c2->cd(3);
  //gPad->SetLogy(true);
  //h = h1_ahdc_cf_time->DrawCopy();
  //c2->cd(4);
  //gPad->SetLogy(true);
  //h = h1_ahdc_edge_time->DrawCopy();
  //c2->SaveAs("results/ahdc_plots.png");


  //TCanvas *c3 = new TCanvas("c3","c3",1200,1000);
  //h2_ahdc_integral_vs_Tot->DrawCopy("colz");
  //c3->SaveAs("results/integral_vs_tot.png");

  //c3 = new TCanvas("c3","c3",1200,1000);
  //h2_ahdc_cf_vs_edge->DrawCopy("colz");
  //c3->SaveAs("results/cf_vs_edge.png");


  ////h->SetLineColor(1);
  ////h = h1_diff0->DrawCopy("same");
  ////h->SetLineColor(2);
  ////h = h1_diff1->DrawCopy("same");
  ////h->SetLineColor(4);
  ////h = h1_diff2->DrawCopy("same");
  ////h->SetLineColor(7);
  ////h = h1_diff3->DrawCopy("same");
  ////h->SetLineColor(8);

  ////c2->cd(2);
  ////h = h1_sum->DrawCopy();
  ////h->SetLineColor(1);
  ////h = h1_sum0->DrawCopy("same");
  ////h->SetLineColor(2);
  ////h = h1_sum1->DrawCopy("same");
  ////h->SetLineColor(4);
  ////h = h1_sum2->DrawCopy("same");
  ////h->SetLineColor(7);
  ////h = h1_sum3->DrawCopy("same");
  ////h->SetLineColor(8);


  ////c2->SaveAs("results/atof_diffsum.png");

  ////auto names2 =  df2.GetColumnNames() ;
  //////  Shows what is in the file.
  ////for( const auto& n : names2) {
  ////  std::cout << n << "\n";
  ////}
  ////df2.Describe().Print();
  ////// auto c1 = df.Count();
  ////df2.Snapshot("data","ahdc_waveforms");

}

