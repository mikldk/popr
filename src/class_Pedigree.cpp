#include "popr_types.hpp"

/*
==========================================
Pedigree
==========================================
*/

Pedigree::Pedigree(int id) {
  //Rcpp::Rcout << "Pedigree with id = " << id << " created" << std::endl;
  
  m_pedigree_id = id;
  
  m_all_individuals = new std::vector<Individual*>();
  m_relations = new std::vector< std::pair<Individual*, Individual*>* >();
}

Pedigree::~Pedigree() {
  delete m_all_individuals;
  
  for (auto it = m_relations->begin(); it != m_relations->end(); ++it) {
    delete *it;
  }
  
  delete m_relations;
  
  this->m_points.clear(); // FIXME??
}

void Pedigree::add_member(Individual* i) {
  m_all_individuals->push_back(i); // so that ones without location are included
  
  // If location is not set, do not include it in the point cloud
  if (!(i->location_is_set())) {
    return;
  }
  
  //Point p = { .x = i->get_etrs89e(), .y = i->get_etrs89n(), .ind = i};
  //this->m_points.push_back(p);
  this->m_points.push_back({ .x = i->get_etrs89e(), .y = i->get_etrs89n(), .ind = i});
}

void Pedigree::add_relation(Individual* lhs, Individual* rhs) {
  //std::pair<Individual*, Individual*>* pair = new std::pair(lhs, rhs);
  std::pair<Individual*, Individual*>* pair = new std::pair<Individual*, Individual*>(lhs, rhs);
  m_relations->push_back(pair);
}

std::vector<Individual*>* Pedigree::get_all_individuals() const {
  return m_all_individuals;
}

std::vector< std::pair<Individual*, Individual*>* >* Pedigree::get_relations() const {
  return m_relations;
}

