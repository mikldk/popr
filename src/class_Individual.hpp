#include "popr_types.hpp"

#include <vector>

class Individual {
private:
  int m_pid; 
  bool m_is_male; // CANNOT BE 0/NA!
  int m_birth_year = 0;
  bool m_is_alive = false;
  
  std::vector<Individual*>* m_children = NULL;
  Individual* m_mother = NULL;
  Individual* m_father = NULL;
  
  Pedigree* m_pedigree = NULL;
  int m_pedigree_id = 0;
  
  double m_etrs89e = std::numeric_limits<double>::quiet_NaN();
  double m_etrs89n = std::numeric_limits<double>::quiet_NaN();
  bool m_etrs89set = false;
  
  void meiosis_dist_tree_internal(Individual* dest, int* dist) const;
  
  bool m_dijkstra_visited = false;
  int m_dijkstra_distance = 0;

  std::vector<int> m_father_haplotype;
  bool m_father_haplotype_set = false;
  bool m_father_haplotype_mutated = false;
  void father_haplotype_mutate(double mutation_rate);
  
public:
  Individual(int pid, bool is_male);
  ~Individual();
  int get_pid() const;
  bool get_is_male() const;
  void add_child(Individual* child);
  void set_mother(Individual* i);
  void set_father(Individual* i);
  Individual* get_mother() const;
  Individual* get_father() const;
  std::vector<Individual*>* get_children() const;
  bool pedigree_is_set() const;
  Pedigree* get_pedigree() const;
  int get_pedigree_id() const;
  
  void set_pedigree_id(int id, Pedigree* ped, int* pedigree_size);
  
  void set_alive_status(bool is_alive);
  bool get_alive_status() const;  

  void set_birth_year(int birth_year);
  int get_birth_year() const;

  void set_location(double etrs89e, double etrs89n);
  
  double get_etrs89e() const;
  double get_etrs89n() const;
  bool location_is_set() const;
  
  int meiosis_dist_tree(Individual* dest) const;
  
  void dijkstra_reset();
  void dijkstra_tick_distance(int step);
  void dijkstra_set_distance_if_less(int dist);
  void dijkstra_mark_visited();
  int dijkstra_get_distance() const;
  bool dijkstra_was_visited() const;
  
  bool is_father_haplotype_set() const;
  void set_father_haplotype(std::vector<int> h);
  std::vector<int> get_father_haplotype() const;
  void pass_haplotype_to_children(bool recursive, double mutation_rate);
};

