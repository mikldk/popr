#include "popr_types.hpp"

#include <vector>

// FIXME: PointCloud<double> vs PointCloud<int> vs PointCloud<int64> ?
class Pedigree : public PointCloud<double> {
private:
  int m_pedigree_id;
  std::vector<Individual*>* m_all_individuals = NULL;
  std::vector< std::pair<Individual*, Individual*>* >* m_relations = NULL; 
  
public:
  Pedigree(int id);
  ~Pedigree();
  int get_id() const;
  void add_member(Individual* i);
  void add_relation(Individual* lhs, Individual* rhs);
  std::vector<Individual*>* get_all_individuals() const;
  std::vector< std::pair<Individual*, Individual*>* >* get_relations() const;
  
  void populate_father_haplotypes(int loci, double mutation_rate);
};

