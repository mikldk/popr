#include "popr_types.hpp"

#include <vector>
#include <unordered_map>

// FIXME: PointCloud<double> vs PointCloud<int> vs PointCloud<int64> ?
class Population : public PointCloud<double> {
private:
  std::unordered_map<int, Individual*>* m_population = NULL;

public:
  Population(std::unordered_map<int, Individual*>* population);
  ~Population();
  void build_pointcloud();
  std::unordered_map<int, Individual*>* get_population() const;
  int get_population_size() const;
  Individual* get_individual(int pid) const;
};

