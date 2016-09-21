#include "popr_types.hpp"

#include <stdexcept>

/*
==========================================
Individual
==========================================
*/
Individual::Individual(int pid, bool is_male) {
  m_pid = pid;
  
  m_is_male = is_male;
  
  m_children = new std::vector<Individual*>();
}

Individual::~Individual() {
  delete m_children;
}

int Individual::get_pid() const {
  return m_pid;
}

/*
int Individual::get_pid_m() const {
  return m_pid_m;
}


int Individual::get_pid_f() const {
  return m_pid_f;
}
*/

bool Individual::get_is_male() const {
  return m_is_male;
}

void Individual::add_child(Individual* child) {
  m_children->push_back(child);
}

void Individual::set_mother(Individual* i) {
  // FIXME: Check sex of i?
  m_mother = i;
}

void Individual::set_father(Individual* i) {
  // FIXME: Check sex of i?
  m_father = i;
}

Individual* Individual::get_mother() const {
  return m_mother;
}

Individual* Individual::get_father() const {
  return m_father;
}

std::vector<Individual*>* Individual::get_children() const {
  return m_children;
}

bool Individual::pedigree_is_set() const {
  return (m_pedigree_id != 0);
}

int Individual::get_pedigree_id() const {
  return m_pedigree_id;
}

Pedigree* Individual::get_pedigree() const {
  return m_pedigree;
}

void Individual::set_pedigree_id(int id, Pedigree* ped, int* pedigree_size) {
  if (this->pedigree_is_set()) {
    return;
  }
  
  m_pedigree = ped;
  m_pedigree_id = id;
  *pedigree_size += 1;
  ped->add_member(this);
  
  if (m_mother != NULL) {
    m_mother->set_pedigree_id(id, ped, pedigree_size);
  }
  
  if (m_father != NULL) {  
    m_father->set_pedigree_id(id, ped, pedigree_size);
  }
  
  for (auto &child : (*m_children)) {
    ped->add_relation(this, child);
    child->set_pedigree_id(id, ped, pedigree_size);
  }
}

void Individual::set_alive_status(bool is_alive) {
  m_is_alive = is_alive;
}

bool Individual::get_alive_status() const {
  return m_is_alive;
}

void Individual::set_birth_year(int birth_year) {
  m_birth_year = birth_year;
}

int Individual::get_birth_year() const {
  return m_birth_year;
}

void Individual::set_location(double etrs89e, double etrs89n) {
  m_etrs89e = etrs89e;
  m_etrs89n = etrs89n;
  m_etrs89set = true;
}

double Individual::get_etrs89e() const {
  return m_etrs89e;
}

double Individual::get_etrs89n() const {
  return m_etrs89n;
}

bool Individual::location_is_set() const {
  return m_etrs89set;
}

void Individual::dijkstra_reset() {
  m_dijkstra_visited = false;
  m_dijkstra_distance = 0;
}

void Individual::dijkstra_tick_distance(int step) {
  m_dijkstra_distance += step;
}

void Individual::dijkstra_set_distance_if_less(int dist) {
  if (m_dijkstra_distance < dist) {
    m_dijkstra_distance = dist;
  }
}

void Individual::dijkstra_mark_visited() {
  m_dijkstra_visited = true;
}

int Individual::dijkstra_get_distance() const {
  return m_dijkstra_distance; 
}

bool Individual::dijkstra_was_visited() const {
  return m_dijkstra_visited; 
}

/*
void Individual::set_pedigree_id_f(int id, Pedigree* ped, int* pedigree_size) {
  if (!(this->get_is_male())) {
    return;
  } 
  
  if (m_pedigree_id_f != 0) {
    return;
  }
  
  m_pedigree_id_f = id;
  *pedigree_size += 1;
  ped->add_member(this);
  
  if (m_father != NULL) {
    m_father->set_pedigree_id_f(id, ped, pedigree_size);
  }
  
  for (auto &child : (*m_children)) {
    ped->add_relation(this, child);
    child->set_pedigree_id_f(id, ped, pedigree_size);
  }
}
*/


// ASSUMES TREE!
//FIXME: Heavily relies on it being a tree, hence there is only one path connecting every pair of nodes
void Individual::meiosis_dist_tree_internal(Individual* dest, int* dist) const {
  if (this->get_pid() == dest->get_pid()) {
    //FIXME: Heavily relies on it being a tree, hence there is only one path connecting every pair of nodes
    *dist = dest->dijkstra_get_distance();
    return;
  }
  
  if (dest->get_mother() != NULL) {
    throw std::invalid_argument("meiosis_dist_tree_internal assumes tree (e.g. Ychr)!");
  }
  
  if (dest->dijkstra_was_visited()) {
    return;
  }
  
  dest->dijkstra_mark_visited();
  dest->dijkstra_tick_distance(1);
  int m = dest->dijkstra_get_distance();
  
  // FIXME: If not tree, then distance must be somehow checked if shorter and then adjusted
  
  Individual* father = dest->get_father();
  if (father != NULL) {  
    //tree: ok
    father->dijkstra_tick_distance(m);

    // general? FIXME Correct?
    //father->dijkstra_set_distance_if_less(m);
    
    this->meiosis_dist_tree_internal(father, dist); 
  }
  
  std::vector<Individual*>* children = dest->get_children();
  for (auto child : *children) {
    //tree: ok
    child->dijkstra_tick_distance(m);

    // general? FIXME Correct?
    //child->dijkstra_set_distance_if_less(m);
    
    this->meiosis_dist_tree_internal(child, dist);
  }
}

// ASSUMES TREE!
int Individual::meiosis_dist_tree(Individual* dest) const {
  if (!(this->pedigree_is_set())) {
    throw std::invalid_argument("!(this->pedigree_is_set())");
  }
  
  if (dest == NULL) {
    throw std::invalid_argument("dest is NULL");
  }
  
  if (!(dest->pedigree_is_set())) {
    throw std::invalid_argument("!(dest->pedigree_is_set())");
  }
  
  if (this->get_pedigree_id() != dest->get_pedigree_id()) {
    return -1;
  }
  
  std::vector<Individual*>* inds = this->get_pedigree()->get_all_individuals();
  for (auto child : *inds) {
    child->dijkstra_reset();
  }
  
  // At this point, the individuals this and dest belong to same pedigree
  int dist = 0;
  this->meiosis_dist_tree_internal(dest, &dist);
  return dist;
}


