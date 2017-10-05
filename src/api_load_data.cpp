#include <sys/resource.h>
#include <sstream>

#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include "popr_types.hpp"

#include "nanoflann.hpp"

using namespace Rcpp;

// [[Rcpp::export]]
void wipe_population(Rcpp::XPtr<Population> population) {
  Population* pop = population;
  delete pop;
}

//' Construct a population
//' 
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<Population> load_individuals(IntegerVector pid, 
                                        LogicalVector is_male, 
                                        IntegerVector pid_mom, 
                                        IntegerVector pid_dad, 
                                        LogicalVector is_alive = LogicalVector::create(), 
                                        IntegerVector birth_year = IntegerVector::create(),
                                        NumericVector etrs89e = NumericVector::create(), 
                                        NumericVector etrs89n = NumericVector::create(),
                                        bool progress = true, 
                                        bool error_on_gender_mismatch = true) {
  std::unordered_map<int, Individual*>* pop = new std::unordered_map<int, Individual*>(); // pid's are garanteed to be unique

  Population* population = new Population(pop);
  Rcpp::XPtr<Population> res(population, true);
  res.attr("class") = CharacterVector::create("popr_abort", "popr_population_abort", "externalptr");
  
  
  int N = pid.size();
  
  if (is_male.size() != N || pid_mom.size() != N || pid_dad.size() != N || 
      (is_alive.size() != 0 && is_alive.size() != N) ||
      (birth_year.size() != 0 && birth_year.size() != N) ||  
      (etrs89e.size() != 0 && etrs89e.size() != N) || 
      (etrs89n.size() != 0 && etrs89n.size() != N)) {
    stop("all vectors (pid, pid_mom, pid_dad, is_alive*, birth_year*, etrs89e*, etrs89n*) must have same length (or length 0 for *-marked)");
  }
  
  if (etrs89e.size() != etrs89n.size()) {
    stop("etrs89e and etrs89n must have same length");
  }
  
  // N for building hashmap, and N for building relations (adding children)
  //Progress p(2*N, true);
  Progress p(2*N, progress);
  
  

  // Build hashmap
  for (int k = 0; k < N; ++k) {
    int i_pid = pid[k];
    bool i_is_male = is_male[k];
    
    Individual* i = new Individual(i_pid, i_is_male);

    if (is_alive.size() == N && !LogicalVector::is_na(is_alive[k])) {
      i->set_alive_status(is_alive[k]);
    }
    
    if (birth_year.size() == N && !IntegerVector::is_na(birth_year[k])) {
      i->set_birth_year(birth_year[k]);
    }
        
    if (etrs89e.size() == N && etrs89n.size() == N && !NumericVector::is_na(etrs89e[k]) && !NumericVector::is_na(etrs89n[k])) {
      i->set_location(etrs89e[k], etrs89n[k]);
    }
    
    (*pop)[i->get_pid()] = i;
    
    if (k % CHECK_ABORT_EVERY == 0 && Progress::check_abort() ) {
      return res;
    }
    
    if (progress) {
      p.increment();
    }
  }
  
  // Fill out parents
  for (int k = 0; k < N; ++k) {
    int i_pid = pid[k];
    int i_pid_mom = pid_mom[k];
    int i_pid_dad = pid_dad[k];
    bool i_is_male = is_male[k];
    
    Individual* i = (*pop)[i_pid];
    
    if (i_pid_mom > 0) {
      std::unordered_map<int, Individual*>::const_iterator m = pop->find(i_pid_mom);
      
      if (m == pop->end()) {
        std::ostringstream err;
        err << "NOT FOUND: pid_mom = " << i_pid_mom << " for pid = " << i_pid << " was not found as a pid itself!";
        stop(err.str());
      } else {
        if (m->second->get_is_male()) {
          std::ostringstream err;
          err << "NOT FEMALE: pid_mom = " << i_pid_mom << " for pid = " << i_pid << " was a male!";
          
          if (error_on_gender_mismatch) {
            stop(err.str());
          } else {
            warning(err.str());
          }
        }
        
        // FIXME: Check gender and age (younger than child?) of mother?
        i->set_mother(m->second);
        m->second->add_child(i);
      }
    }
    
    if (i_pid_dad > 0) {
      std::unordered_map<int, Individual*>::const_iterator f = pop->find(i_pid_dad);    
      if (f == pop->end()) {
        std::ostringstream err;
        err << "NOT FOUND: pid_dad = " << i_pid_dad << " for pid = " << i_pid << " was not found as a pid itself!";
        stop(err.str());
      } else {
        if (!(f->second->get_is_male())) {
          std::ostringstream err;
          err << "NOT MALE: pid_dad = " << pid_dad << " for pid = " << i_pid << " was not a male!";
          
          if (error_on_gender_mismatch) {
            stop(err.str());
          } else {
            warning(err.str());
          }
        }
        
        i->set_father(f->second);
        f->second->add_child(i);
      }
    }
    
    if (k % CHECK_ABORT_EVERY == 0 && Progress::check_abort() ) {      
      return res;
    }
    
    if (progress) {
      p.increment();
    }
  }
  
  if (etrs89e.size() != 0 && etrs89n.size() != 0) {
    population->build_pointcloud();
  }

  res.attr("class") = CharacterVector::create("popr_population", "externalptr");
  
  return res;
}

