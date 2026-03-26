//******************************************************************
//*  ██╗  ██╗██╗██████╗  ██████╗     ██╗  ██╗    ██████╗
//*  ██║  ██║██║██╔══██╗██╔═══██╗    ██║  ██║   ██╔═████╗
//*  ███████║██║██████╔╝██║   ██║    ███████║   ██║██╔██║
//*  ██╔══██║██║██╔═══╝ ██║   ██║    ╚════██║   ████╔╝██║
//*  ██║  ██║██║██║     ╚██████╔╝         ██║██╗╚██████╔╝
//*  ╚═╝  ╚═╝╚═╝╚═╝      ╚═════╝          ╚═╝╚═╝ ╚═════╝
//************************ Jefferson National Lab (2017) ***********
//******************************************************************
//* Example program for reading HIPO-4 Files..
//* Reads the file created by writeFile program
//*--
//* Author: G.Gavalian
//*

#include <cstdlib>
#include <iostream>
#include "hipo4/reader.h"



int main(int argc, char** argv) {

   std::cout << " reading file example program (HIPO) "  << __cplusplus << std::endl;

   char inputFile[256];

   if(argc>1) {
     snprintf(inputFile,256,"%s",argv[1]);
      //sprintf(outputFile,"%s",argv[2]);
   } else {
      std::cout << " *** please provide a file name..." << std::endl;
     exit(1);
   }

   hipo::reader  reader;
   reader.open(inputFile);
   hipo::dictionary  factory;
   reader.readDictionary(factory);

   factory.show();

   hipo::bank  atof_tdc(factory.getSchema("ATOF::tdc"));

   hipo::event      event;

   int counter = 0;

   reader.gotoEvent(2);

   while(reader.next()==true){
      reader.read(event);

      event.getStructure(atof_tdc);

      atof_tdc.show();

      int nrows = atof_tdc.getRows();
      printf("---------- ATOF TDC BANK CONTENT -------\n");
      for(int row = 0; row < nrows; row++){
         int  s   = atof_tdc.getShort("sector",row);
         int  l   = atof_tdc.getShort("layer",row);
         int  c   = atof_tdc.getShort("component",row);
         int  o   = atof_tdc.getShort("order",row);
         int tdc  = atof_tdc.getInt("TDC",row);
         if(tdc != 0)
         printf("%3d %6d %6d %6d : %8.4f\n",s,l,c,o,tdc);
      }
      printf("---------- END OF ATOF TDC BANK -------\n");

      counter++;
   }
   printf("processed events = %d\n",counter);
}
//### END OF GENERATED CODE
