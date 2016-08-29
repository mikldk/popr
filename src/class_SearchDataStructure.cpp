#include <Rcpp.h>
#include "popr_types.hpp"

SearchDataStructure::SearchDataStructure(kd_tree* tree, PointCloud<double>* cloud) {
  m_tree = tree;
  m_cloud = cloud;
}

SearchDataStructure::~SearchDataStructure() {
  
}

kd_tree* SearchDataStructure::get_tree() const {
  return m_tree;
}

PointCloud<double>* SearchDataStructure::get_cloud() const {
  return m_cloud;
}
