#include <Rcpp.h>
#include "popr_types.hpp"

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

/*
==========================================
Individual
==========================================
*/
Population::Population(std::unordered_map<int, Individual*>* population) {
  m_population = population;  
}

Population::~Population() {
  //delete m_den_kd_tree;
  //delete m_cloud;
  
  std::unordered_map<int, Individual*> pop = *m_population;
  
  for (auto it = pop.begin(); it != pop.end(); ++it) {
    delete (it->second);
  }
  
  delete m_population;
}

void Population::build_pointcloud() {
  this->m_points.clear();
  this->m_points.resize(m_population->size());
  
  for (auto i : *m_population) {
    // If location is not set, do not include it in the point cloud
    if (!(i.second->location_is_set())) {
      continue;
    }
    
    //Point p = { .x = i.second->get_etrs89e(), .y = i.second->get_etrs89n(), .ind = i.second};
    //this->m_points.push_back(p);
    this->m_points.push_back({ .x = i.second->get_etrs89e(), .y = i.second->get_etrs89n(), .ind = i.second});
  }
}

std::unordered_map<int, Individual*>* Population::get_population() const {
  return m_population;
}

Individual* Population::get_individual(int pid) const {
  std::unordered_map<int, Individual*>::const_iterator got = m_population->find(pid);
  
  if (got == m_population->end()) {
    Rcpp::Rcerr << "Individual with pid = " << pid << " not found!" << std::endl;
    Rcpp::stop("Individual not found");
  }
  
  return got->second;
}

int Population::get_population_size() const {
  return m_population->size();
}
