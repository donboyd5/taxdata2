
#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
DataFrame matchrecs(DataFrame adf, DataFrame bdf) {
  
  // Extracting columns
  IntegerVector ida = adf["ida"];
  NumericVector weighta = adf["weighta"];
  NumericVector ranka = adf["ranka"];
  
  IntegerVector idb = bdf["idb"];
  NumericVector weightb = bdf["weightb"];
  NumericVector rankb = bdf["rankb"];
  
  int ia = 0; // 0-based in C++
  int ib = 0;

  double dweighta = weighta[0]; // seed with info from the first record
  double dweightb = weightb[0];
  
  std::vector<int> records_ida, records_idb;
  std::vector<double> records_weighta, records_weightb, records_ranka, records_rankb;
  
  while (ia < ida.size() && ib < idb.size()) {
    if (dweightb > dweighta) { // done with the a record so write result
      records_ida.push_back(ida[ia]);
      records_idb.push_back(idb[ib]);
      records_weighta.push_back(dweighta);
      records_weightb.push_back(dweighta);
      records_ranka.push_back(ranka[ia]);
      records_rankb.push_back(rankb[ib]);
      
      dweightb = dweightb - dweighta;
      ia++;
      dweighta = weighta[ia]; // get the next weighta
      
    } else {
      records_ida.push_back(ida[ia]);
      records_idb.push_back(idb[ib]);
      records_weighta.push_back(dweightb);
      records_weightb.push_back(dweightb);
      records_ranka.push_back(ranka[ia]);
      records_rankb.push_back(rankb[ib]);
      
      dweighta = dweighta - dweightb;
      ib++;
      dweightb = weightb[ib]; // get the next weightb
    }
  }
  
  return DataFrame::create(_["ida"]=records_ida, _["idb"]=records_idb, _["weighta"]=records_weighta, 
                          _["weightb"]=records_weightb, _["ranka"]=records_ranka, _["rankb"]=records_rankb);
}


