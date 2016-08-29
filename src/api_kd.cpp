#include <Rcpp.h>

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>


#include "popr_types.hpp"

SearchDataStructure* build_kdtree_from_pointcloud(PointCloud<double>* pointcloud, int max_leaf_size) {
  kd_tree* index = new kd_tree(2, *pointcloud, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size /* max leaf */));
  index->buildIndex();
  
  SearchDataStructure* search_tree = new SearchDataStructure(index, pointcloud);
  
  return search_tree;
}

//' Build kd-tree with all individuals from population
//' 
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<SearchDataStructure> build_kdtree_from_population(Rcpp::XPtr<Population> population, int max_leaf_size) {
  SearchDataStructure* search_tree = build_kdtree_from_pointcloud(population, max_leaf_size);
  Rcpp::XPtr<SearchDataStructure> res(search_tree, false);
  res.attr("class") = CharacterVector::create("popr_search_data_structure", "externalptr");
  
  return res;
}

//' Build kd-tree with all individuals from pedigree
//' 
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<SearchDataStructure> build_kdtree_from_pedigree(Rcpp::XPtr<Pedigree> pedigree, int max_leaf_size) {
  SearchDataStructure* search_tree = build_kdtree_from_pointcloud(pedigree, max_leaf_size);
  Rcpp::XPtr<SearchDataStructure> res(search_tree, false);
  res.attr("class") = CharacterVector::create("popr_search_data_structure", "externalptr");
  
  return res;
}


//' Build kd-tree from a list of pids
//' 
//' @export
// [[Rcpp::export]]
Rcpp::XPtr<SearchDataStructure> build_kdtree_from_pids(IntegerVector pids, Rcpp::XPtr<Population> population, int max_leaf_size) {
  Population* pop = population;
  int N = pids.size();
  int k = 0;
  Progress p(N, true);
  PointCloud<double>* cloud = new PointCloud<double>();
  cloud->m_points.resize(N);
  
  kd_tree* index = new kd_tree(2, *cloud, nanoflann::KDTreeSingleIndexAdaptorParams(max_leaf_size /* max leaf */));
  
  SearchDataStructure* search_tree = new SearchDataStructure(index, cloud);
  Rcpp::XPtr<SearchDataStructure> res(search_tree, false);
  res.attr("class") = CharacterVector::create("popr_search_data_structure", "externalptr");
  
  for (size_t i = 0; i < N; ++i) {
    Individual* ind = pop->get_individual(pids[i]);
    
    if (!(ind->location_is_set())) {
      continue;
    }
    
    cloud->m_points[k].x = ind->get_etrs89e();
    cloud->m_points[k].y = ind->get_etrs89n();
    cloud->m_points[k].ind = ind;
    ++k;
    
    if (k % CHECK_ABORT_EVERY == 0 && Progress::check_abort() ) {
      return res;
    }
    
    p.increment();
  }
  
  index->buildIndex();
  
  return res;
}

// Old, searches kd for each radius
//IntegerMatrix analyse_meioses_internal(Rcpp::XPtr<Population> population, IntegerVector pids, NumericVector radii) {


// [[Rcpp::export]]
IntegerMatrix analyse_meioses_search_tree(Rcpp::XPtr<Population> population, Rcpp::XPtr<SearchDataStructure> search_tree, IntegerVector pids, NumericVector radii, bool report_progress = false) {
  if (!std::is_sorted(radii.begin(), radii.end(), std::greater<double>())) {
    stop("Expected radii to be sorted in decreasing order");
  }
  
  Population* pop = population;

  /* **Important note:** 
  If L2 norms are used, notice that search radius and all passed and 
  returned distances are actually *squared distances*.
  */
  const double largest_radius = radii[0]; // sorted desc
  const double search_radius = static_cast<double>(largest_radius*largest_radius);
  nanoflann::SearchParams params;
  params.sorted = true; // Important! Results sorted by ascending distances  
  
  PointCloud<double>* cloud = search_tree->get_cloud();
  kd_tree* tree = search_tree->get_tree();
  
  if (cloud == NULL || tree == NULL) {
    stop("Be sure to build_kd_tree before trying to search it...");
  }
  
  //Rcout << "cloud size = " <<  search_tree->get_cloud()->m_points.size() << std::endl;
  
  int n_radii = radii.size();
  int n_pids = pids.size();
  
  std::vector< std::vector< std::unordered_map<int, int> > > res1;
  res1.reserve(n_pids);
  
  Progress p(2*n_pids, report_progress);
  
  size_t rows = 0;
  
  for (int i1 = 0; i1 < n_pids; ++i1) {
    Individual* ind = pop->get_individual(pids[i1]);

    if (!(ind->location_is_set())) {
      Rcpp::stop("Individual searched for does not have a location.");
    }

    int search_pedid = ind->get_pedigree_id();
          
    if (search_pedid == 0) {
      Rcpp::stop("Individual searched for is not in a pedigree, not even its own. Make sure that build_pedigrees has been called.");
    }

    std::vector< std::unordered_map<int, int> > res2;
    res2.reserve(n_radii);
          
    const double etrs89e = ind->get_etrs89e();
    const double etrs89n = ind->get_etrs89n();  
    const double query_pt[2] = { etrs89e, etrs89n };
    std::vector<std::pair<size_t, double>> ret_matches;

    // Results sorted by ascending distances  
    const size_t n_matches = tree->radiusSearch(&query_pt[0], search_radius, ret_matches, params);
    
    // treat i2 = 0 explicitely as it is known that all individuals are within this radius
    std::unordered_map<int, int> dist_tab_0;
    for (size_t i = 0; i < n_matches; ++i) {
      size_t index = ret_matches[i].first;
      Individual* match = cloud->m_points[index].ind;
      // FIXME: if (match == NULL)?
      
      const double match_dist = sqrt(ret_matches[i].second);
      int dist = -1;
      
      if (search_pedid == match->get_pedigree_id()) {
        dist = ind->meiosis_dist_tree(match);
      }
      
      dist_tab_0[dist] += 1;
    }
    
    res2.push_back(dist_tab_0);
    rows += dist_tab_0.size();
    
    // i2 = 0 treated explicitely above      
    for (int i2 = 1; i2 < n_radii; ++i2) {
      const double radius = radii[i2];
      const double radius_sq = radius*radius;
      
      std::unordered_map<int, int> dist_tab_i2;
      // search results are in ascinding order, break if a match is too far away (then the rest will also be too far away)
      for (size_t i = 0; i < n_matches; ++i) {
        //const double match_dist = sqrt(ret_matches[i].second);
        
        //if (match_dist >= radius) {
        if (ret_matches[i].second >= radius_sq) {
          break; // match distance is greater than radius currently considered, stop
        }
        
        size_t index = ret_matches[i].first;
        Individual* match = cloud->m_points[index].ind;
        // FIXME: if (match == NULL)?

        int dist = -1;
        
        if (search_pedid == match->get_pedigree_id()) {
          dist = ind->meiosis_dist_tree(match);
        }
        
        dist_tab_i2[dist] += 1;
      }
      
      res2.push_back(dist_tab_i2);
      rows += dist_tab_i2.size();
    }
    
    res1.push_back(res2);
    
    if (Progress::check_abort()) {
      IntegerMatrix empty(1, 1);
      empty(0, 0) = NA_INTEGER;
      return empty;
    }
    
    if (report_progress) {
      p.increment();
    }
  }
  
  IntegerMatrix res_all(rows, 4);
  
  int k = 0;
  
  for (int i1 = 0; i1 < n_pids; ++i1) {  
    std::vector< std::unordered_map<int, int> > res2 = res1.at(i1);

    for (int i2 = 0; i2 < n_radii; ++i2) {
      int radius_index = i2 + 1; // radius index with R 1-indexed vectors
      std::unordered_map<int, int> res = res2.at(i2);
      
      for (auto res_entry : res) {
        res_all(k, 0) = pids[i1];
        res_all(k, 1) = radius_index;
        res_all(k, 2) = res_entry.first;
        res_all(k, 3) = res_entry.second;
        ++k;
      }
    }
    
    if (Progress::check_abort()) {
      IntegerMatrix empty(1, 1);
      empty(0, 0) = NA_INTEGER;
      return empty;
    }
    
    if (report_progress) {
      p.increment();
    }
  }
  
  return res_all;
}



























// [[Rcpp::export]]
IntegerVector get_pids_within_radius(Rcpp::XPtr<Individual> ind, Rcpp::XPtr<SearchDataStructure> search_tree, double radius) {
  /* **Important note:** 
  If L2 norms are used, notice that search radius and all passed and 
  returned distances are actually *squared distances*.
  */
  const double search_radius = static_cast<double>(radius*radius);
  nanoflann::SearchParams params;
  
  PointCloud<double>* cloud = search_tree->get_cloud();
  kd_tree* tree = search_tree->get_tree();
  
  if (!(ind->location_is_set())) {
    Rcpp::stop("Individual searched for does not have a location.");
  }

  int search_pedid = ind->get_pedigree_id();
        
  if (search_pedid == 0) {
    Rcpp::stop("Individual searched for is not in a pedigree, not even its own. Make sure that build_pedigrees has been called.");
  }
  
  const double etrs89e = ind->get_etrs89e();
  const double etrs89n = ind->get_etrs89n();  
  const double query_pt[2] = { etrs89e, etrs89n };
  std::vector<std::pair<size_t, double>> ret_matches;

  // Results sorted by ascending distances  
  const size_t n_matches = tree->radiusSearch(&query_pt[0], search_radius, ret_matches, params);
  
  IntegerVector pids(n_matches);
    
  for (size_t i = 0; i < n_matches; ++i) {
    size_t index = ret_matches[i].first;
    Individual* match = cloud->m_points[index].ind;
    pids[i] = match->get_pid();
  }
  
  return pids;
}





























/**
 * A result-set class used when performing a radius based search.
 */
template <typename DistanceType, typename IndexType = size_t>
class RadiusResultSetCountOnly
{
public:
  const DistanceType radius;
  size_t m_hits;
  inline RadiusResultSetCountOnly(DistanceType radius_) : radius(radius_), m_hits(0) {  }
  inline ~RadiusResultSetCountOnly() { }
  inline size_t size() const { return m_hits; }
  inline bool full() const { return true; }
  inline void addPoint(DistanceType dist, IndexType index) {  if (dist<radius) ++m_hits; }
  inline DistanceType worstDist() const { return radius; }
};
/** @} */

/*
// [[Rcpp::export]]
IntegerMatrix analyse_meioses_internal_tree_new(Rcpp::XPtr<Population> population, Rcpp::XPtr<SearchDataStructure> search_tree, IntegerVector pids, NumericVector radii) {
  PointCloud<double>* cloud = search_tree->get_cloud();
  kd_tree* tree = search_tree->get_tree();
  
  double radius = 10;
  nanoflann::SearchParams params;
  
  RadiusResultSetCountOnly<double, size_t> result_set(radius, IndicesDists);
  const double query_pt[2] = { etrs89e, etrs89n };
  
  const size_t n_matches = tree->radiusSearchCustomCallback(&query_pt[0], result_set, params);
  
  
}
*/

// Note, that this includes the ones in the pedigree! Hence, remember to substract these afterwards!
// [[Rcpp::export]]
IntegerMatrix radius_search_count(Rcpp::XPtr<Population> population, Rcpp::XPtr<SearchDataStructure> search_tree, IntegerVector pids, NumericVector radii, bool report_progress = false) {
  Population* pop = population;
  
  nanoflann::SearchParams params;
  
  PointCloud<double>* cloud = search_tree->get_cloud();
  kd_tree* tree = search_tree->get_tree();
  
  if (cloud == NULL || tree == NULL) {
    stop("Be sure to build_kd_tree before trying to search it...");
  }
  
  int n_radii = radii.size();
  int n_pids = pids.size();
  
  std::vector< std::vector< std::unordered_map<int, int> > > res1;
  res1.reserve(n_pids);
  
  Progress p(2*n_pids, report_progress);
  
  size_t rows = 0;
  
  for (int i1 = 0; i1 < n_pids; ++i1) {
    Individual* ind = pop->get_individual(pids[i1]);
    
    if (!(ind->location_is_set())) {
      Rcpp::stop("Individual searched for does not have a location.");
    }
    
    int search_pedid = ind->get_pedigree_id();
    
    if (search_pedid == 0) {
      Rcpp::stop("Individual searched for is not in a pedigree, not even its own. Make sure that build_pedigrees has been called.");
    }
     
    std::vector< std::unordered_map<int, int> > res2;
    res2.reserve(n_radii);
    
    const double etrs89e = ind->get_etrs89e();
    const double etrs89n = ind->get_etrs89n();  
    const double query_pt[2] = { etrs89e, etrs89n };
    
    for (int i2 = 0; i2 < n_radii; ++i2) {
      const double radius = radii[i2];
      const double radius_sq = radius*radius;
      
      // Results sorted by ascending distances  
      RadiusResultSetCountOnly<double, size_t> result_set(radius_sq);
      const size_t n_matches = tree->radiusSearchCustomCallback(&query_pt[0], result_set, params);
      
      // treat i2 = 0 explicitely as it is known that all individuals are within this radius
      std::unordered_map<int, int> dist_tab_0;
      dist_tab_0[-1] += n_matches;
      
      res2.push_back(dist_tab_0);
      rows += dist_tab_0.size();
    }
    
    res1.push_back(res2);
    
    if (Progress::check_abort()) {
      IntegerMatrix empty(1, 1);
      empty(0, 0) = NA_INTEGER;
      return empty;
    }
    
    if (report_progress) {
      p.increment();
    }
  }
  
  IntegerMatrix res_all(rows, 4);
  
  int k = 0;
  
  for (int i1 = 0; i1 < n_pids; ++i1) {  
    std::vector< std::unordered_map<int, int> > res2 = res1.at(i1);
    
    for (int i2 = 0; i2 < n_radii; ++i2) {
      int radius_index = i2 + 1; // radius index with R 1-indexed vectors
      std::unordered_map<int, int> res = res2.at(i2);
      
      for (auto res_entry : res) {
        res_all(k, 0) = pids[i1];
        res_all(k, 1) = radius_index;
        res_all(k, 2) = res_entry.first;
        res_all(k, 3) = res_entry.second;
        ++k;
      }
    }
    
    if (Progress::check_abort()) {
      IntegerMatrix empty(1, 1);
      empty(0, 0) = NA_INTEGER;
      return empty;
    }
    
    if (report_progress) {
      p.increment();
    }
  }
  
  return res_all;
  
}

