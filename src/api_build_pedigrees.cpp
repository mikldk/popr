#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include "popr_types.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
void wipe_pedigrees(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees) {
  std::vector<Pedigree*>* peds = pedigrees;
  
  for (auto it = pedigrees->begin(); it != pedigrees->end(); ++it) {
    delete *it;
  }
  
  delete peds;
}

bool pedigree_size_comparator(Pedigree* p1, Pedigree* p2) { 
  return (p1->get_all_individuals()->size() > p2->get_all_individuals()->size());
}

//' Build pedigrees
//' 
//' @export
// [[Rcpp::export]]
Rcpp::XPtr< std::vector<Pedigree*> > build_pedigrees(Rcpp::XPtr<Population> population, bool progress = true) {
  //Rcpp::Rcout << "Bulding pedigrees!" << std::endl;
  
  // Construct pedigrees
  //std::cout << "Starts giving pedigrees ids..." << std::endl;

  std::vector<Pedigree*>* pedigrees = new std::vector<Pedigree*>();
  Rcpp::XPtr< std::vector<Pedigree*> > res(pedigrees, true);
  //Rcpp::XPtr< std::vector<Pedigree*> > res(pedigrees, false);
  res.attr("class") = CharacterVector::create("popr_abort", "popr_pedigreelist_abort", "externalptr");

  int pedigree_id = 1;
  Pedigree* ped;
  
  std::unordered_map<int, Individual*> pop = *(population->get_population());
  int N = pop.size();
  int k = 0;
  Progress p(N, progress);

  ped = new Pedigree(pedigree_id);
  
  for (auto it = pop.begin(); it != pop.end(); ++it) {
    if (it->second->pedigree_is_set()) {
      continue;
    }
    
    int ped_size = 0;
    it->second->set_pedigree_id(pedigree_id, ped, &ped_size);

    pedigree_id += 1;
    
    pedigrees->push_back(ped);
    ped = new Pedigree(pedigree_id);
    
    if (k % CHECK_ABORT_EVERY == 0 && Progress::check_abort() ) {
      return res;
    }
    
    if (progress) {
      p.increment();
    }
    
    ++k;
  }
  
  if (ped->get_all_individuals()->size() > 0) {
    pedigrees->push_back(ped);
  } else {
    delete ped;
  }
  
  std::sort(pedigrees->begin(), pedigrees->end(), pedigree_size_comparator);
  
  res.attr("class") = CharacterVector::create("popr_pedigreelist", "externalptr");

  return res;
}


/*
// [[Rcpp::export]]
IntegerVector build_pedigrees(Rcpp::XPtr<Population> population) {
  // Construct pedigrees
  //std::cout << "Starts giving pedigrees ids..." << std::endl;
  
  std::vector<Pedigree*>* pedigrees = new std::vector<Pedigree*>();

  int pedigree_id = 1;
  Pedigree* ped;
  
  std::unordered_map<int, Individual*> pop = *(population->get_population());
  int N = pop.size();
  int k = 0;
  Progress p(N, true);
  
  IntegerVector pedigree_ids(N);

  ped = new Pedigree(pedigree_id);
  
  for (auto it = pop.begin(); it != pop.end(); ++it) {
    if (it->second->pedigree_is_set()) {
      continue;
    }
    
    int ped_size = 0;
    it->second->set_pedigree_id(pedigree_id, ped, &ped_size);
    
    // prepare for next pedigree
    pedigree_id += 1;
    
    pedigrees->push_back(ped);
    ped = new Pedigree(pedigree_id);
    
    if (k % CHECK_ABORT_EVERY == 0 && Progress::check_abort() ) {
      IntegerVector res(0);
      return res;
    }
    
    p.increment();
    ++k;
  }
  
  k = 0;
  for (auto it = pop.begin(); it != pop.end(); ++it) {
    pedigree_ids[k] = it->second->get_pedigree_id();
    ++k;
  }
  
  if (ped->get_all_individuals()->size() > 0) {
    pedigrees->push_back(ped);
  } else {
    delete ped;
  }
  
  //std::sort(pedigrees->begin(), pedigrees->end(), pedigree_size_comparator);
  
  // FIXME:
  //population->set_pedigrees();
  
  return pedigree_ids;
}
*/


