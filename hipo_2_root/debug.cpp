#include "hipo4/RHipoDS.hxx"

void debug() {

 auto df = MakeHipoDataFrame("alert_test.hipo");
 //auto names =  df.GetColumnNames() ;
 // Shows what is in the file.
 //for( const auto& n : names) {
 //  std::cout << n << "\n";
 //}
 df.Describe().Print();
 auto c1 = df.Count();
 std::cout << *c1 << " events\n";


 //auto df1 = df.Alias("s",  "ATOF_tdc_sector"   )
 //             .Alias("l",  "ATOF_tdc_layer"    )
 //             .Alias("c",  "ATOF_tdc_component")
 //             .Alias("o",  "ATOF_tdc_order"    )
 //             .Alias("tdc","ATOF_tdc_TDC"      )
 //             .Alias("tot","ATOF_tdc_ToT"      );

 auto disp =
     df.Filter("ATOF_tdc_layer.size() > 0")
         //.Display({"s","l","c","o","tdc","tot"},50,48);
         .Display({"ATOF_tdc_sector", "ATOF_tdc_layer", "ATOF_tdc_component",
                   "ATOF_tdc_order", "ATOF_tdc_TDC", "ATOF_tdc_ToT"},
                  50, 48);
 disp->Print();


}
