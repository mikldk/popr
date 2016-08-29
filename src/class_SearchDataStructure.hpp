#include <Rcpp.h>
#include "popr_types.hpp"

class SearchDataStructure {
private:
  kd_tree* m_tree;
  PointCloud<double>* m_cloud;
public:
  SearchDataStructure(kd_tree* tree, PointCloud<double>* cloud);
  ~SearchDataStructure();
  kd_tree* get_tree() const;
  PointCloud<double>* get_cloud() const;
};
