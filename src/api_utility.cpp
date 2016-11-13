#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "popr_types.hpp"

//[[Rcpp::export]]
void popr_test() {
  Rcout << "mikl was here 1324" << std::endl;
}


//[[Rcpp::export]]
int pop_size(Rcpp::XPtr<Population> population) {
  return population->get_population_size();
}

/*
//[[Rcpp::export]]
int get_pid(Rcpp::XPtr<Individual> ind) {
  return ind->get_pid();
}
*/

//' Get number of children for each individual in the population
//' 
//' @export
// [[Rcpp::export]]
IntegerMatrix get_number_of_children(Rcpp::XPtr<Population> population, bool progress = true) {
  std::unordered_map<int, Individual*> pop = *(population->get_population());
  size_t N = pop.size();
  size_t i = 0;
  IntegerMatrix children_count(N, 3);
  colnames(children_count) = CharacterVector::create("pid", "boys_n", "girls_n");

  Progress p(N, progress);
  
  for (auto it = pop.begin(); it != pop.end(); ++it) {
    Individual* ind = it->second;
    children_count(i, 0) = ind->get_pid();
    
    std::vector<Individual*>* children = ind->get_children();
    
    for (auto &child : (*children)) {
      if (child->get_is_male()) {
        children_count(i, 1) += 1;
      } else {
        children_count(i, 2) += 1;
      }
    }
    
    ++i;
    
    if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      return children_count;
    }
    
    if (progress) {
      p.increment();
    }
  }
  
  if (i != N) {
    stop("Expected N indviduals...");
  }
  
  return children_count;
}

//' Get ages children
//'
//' @param pid_birthyear Matrix where 1st column is pid and 2nd column is birthyear
//'  
//' @export
// [[Rcpp::export]]
List get_male_children_pids_birthyear(Rcpp::XPtr<Population> population, IntegerVector pids, IntegerMatrix pid_birthyear, bool progress = true) {
  size_t N = pids.size();
  List children_pids(N);
  Progress p(N, progress);
  
  // Fill pid => birthyear mapping
  std::unordered_map<int, int> pid_birthyear_table;
  for (size_t i = 0; i < pid_birthyear.nrow(); ++i) {
    pid_birthyear_table[pid_birthyear(i, 0)] = pid_birthyear(i, 1);
  }
  
  // get children
  
  for (size_t i = 0; i < N; ++i) {
    int pid = pids[i];

    std::unordered_map<int, int>::const_iterator got = pid_birthyear_table.find(pid);
    if (got == pid_birthyear_table.end()) {
      Rcpp::Rcerr << "birthyear for individual with pid = " << pid << " not found!" << std::endl;
      Rcpp::stop("birthyear for individual not found");
    }
    int birthyear = got->second;
        
    Individual* indv = population->get_individual(pid);
    std::vector<Individual*>* children = indv->get_children();
    
    IntegerVector c_pids;
    IntegerVector c_birthyears;
    
    for (auto &child : (*children)) {
      if (child->get_is_male()) {
        int child_pid = child->get_pid();
        
        std::unordered_map<int, int>::const_iterator got_child = pid_birthyear_table.find(child_pid);
        if (got_child == pid_birthyear_table.end()) {
          Rcpp::Rcerr << "birthyear for individual with pid = " << pid << " not found!" << std::endl;
          Rcpp::stop("birthyear for individual not found");
        }
        int birthyear_child = got_child->second;

        c_pids.push_back(child_pid);
        c_birthyears.push_back(birthyear_child);                
      }
    }
    
    List i_info;
    i_info["parent_pid"] = pid;
    i_info["parent_birthyear"] = birthyear;
    i_info["children_pid"] = c_pids;
    i_info["children_birthyears"] = c_birthyears;
    
    children_pids[i] = i_info;
    
    if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort() ) {
      return children_pids;
    }
    
    if (progress) {
      p.increment();
    }
  }
  
  return children_pids;
}

/*
//[[Rcpp::export]]
Rcpp::List pedigrees_to_list(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees) {
  std::vector<Pedigree*>* peds = pedigrees;

  int N = peds->size();
  int k = 0;
  Progress progress(N, true);

  Rcpp::List l(N);
  
  for (auto it = peds->begin(); it != peds->end(); ++it) {
    //Rcpp::XPtr<Pedigree> p(*it, true);
    Rcpp::XPtr<Pedigree> p(*it, false);
    p.attr("class") = Rcpp::CharacterVector::create("popr_pedigree", "externalptr");
    l.push_back(p);
    
    if (k % (N / 10000) == 0 && Progress::check_abort() ) {
      Rcpp::List l_abort(0);
      l_abort.attr("class") = CharacterVector::create("popr_abort", "popr_pedigrees_to_list_abort", "list");
      return l;
    }
    
    progress.increment();
    ++k;
  }
  
  l.attr("class") = CharacterVector::create("popr_pedigrees_list", "list");
  
  return l;
}
*/

//' Get number of pedigrees
//' 
//' @export
// [[Rcpp::export]]
int pedigrees_count(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees) {
  return pedigrees->size();
}

//' Get pedigree size
//' 
//' @export
// [[Rcpp::export]]
int pedigree_size(Rcpp::XPtr<Pedigree> ped) {  
  return ped->get_all_individuals()->size();
}

//' Get distribution of pedigree sizes
//' 
//' @export
//[[Rcpp::export]]
std::unordered_map<int, int> pedigrees_table(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees) {
  std::vector<Pedigree*>* peds = pedigrees;
  std::unordered_map<int, int> tab;
  
  for (auto it = peds->begin(); it != peds->end(); ++it) {
    tab[(*it)->get_all_individuals()->size()] += 1;
  }
  
  return tab;
}

//[[Rcpp::export]]
Rcpp::XPtr<Pedigree> get_pedigree(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int index) {  
  std::vector<Pedigree*>* peds = pedigrees;
  Pedigree* p = peds->at(index);
  
  //Rcpp::XPtr<Pedigree> res(p, true);
  Rcpp::XPtr<Pedigree> res(p, false); // do NOT delete pedigree when not used any more, it still exists in list of pedigrees etc.!
  res.attr("class") = CharacterVector::create("popr_pedigree", "externalptr");
  
  return res;
}

//[[Rcpp::export]]
Rcpp::XPtr<Individual> get_individual(Rcpp::XPtr<Population> population, int pid) {  
  Population* pop = population;
  
  Individual* ind = population->get_individual(pid);
  //Rcpp::XPtr<Individual> res(ind, true);
  Rcpp::XPtr<Individual> res(ind, false); // do NOT delete individual when not used any more, it still exists in pedigree and population etc.!
  res.attr("class") = CharacterVector::create("popr_individual", "externalptr");
  
  return res;
}

//[[Rcpp::export]]
void print_individual(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;
  
  int pid_m = (i->get_mother() != NULL) ? i->get_mother()->get_pid() : 0;
  int pid_f = (i->get_father() != NULL) ? i->get_father()->get_pid() : 0;
  std::vector<Individual*>* children = i->get_children();
  
  Rcpp::Rcout << "  " << i->get_pid() << " with mother  = " << pid_m << " and father " << pid_f << " and children (n = " << children->size() << ")";
  
  if (children->size() == 0) {
    Rcpp::Rcout << std::endl;
  } else {
    Rcpp::Rcout << ": " << std::endl;

    for (auto child : *children) {    
      std::vector<Individual*>* child_children = child->get_children();
      
      Rcpp::Rcout << "  " << child->get_pid() << " with mother  = " << pid_m << " and father " << pid_f << " and " <<  child_children->size() << " children" << std::endl;
    }
  }
}


//[[Rcpp::export]]
void print_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  std::vector< std::pair<Individual*, Individual*>* >* rels = p->get_relations();
  
  Rcpp::Rcout << "Pedigree with " << p->get_all_individuals()->size() << " individuals:" << std::endl;
  
  for (auto i : *inds) {    
    int pid_m = (i->get_mother() != NULL) ? i->get_mother()->get_pid() : 0;
    int pid_f = (i->get_father() != NULL) ? i->get_father()->get_pid() : 0;
    
    Rcpp::Rcout << "  " << i->get_pid() << " with mother  = " << pid_m << " and father " << pid_f << std::endl;
  } 
}

//' get pids in pedigree
//' 
//' @export
// [[Rcpp::export]]
IntegerVector get_pids_in_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  IntegerVector res(inds->size());
  int i = 0;
  for (auto ind : *inds) {   
    res(i) = ind->get_pid();
    ++i;
  } 
  
  return res;
}

//' get pids in pedigree with certain criteria
//' 
//' @export
// [[Rcpp::export]]
IntegerVector get_pids_in_pedigree_criteria(Rcpp::XPtr<Pedigree> ped, bool must_be_alive, bool use_birth_year, int birth_year_min, int birth_year_max) {    
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  IntegerVector res;

  for (auto ind : *inds) {   
    bool skip = false;
    
    if (must_be_alive && ind->get_alive_status() == false) {
      skip = true;
    }
    
    if (!skip && use_birth_year) {
      int birth_year = ind->get_birth_year();
      
      if (birth_year < birth_year_min || birth_year > birth_year_max) {
        skip = true;
      }
    }
    
    if (skip) {
      continue;
    }
  
    res.push_back(ind->get_pid());
  } 
  
  return res;
}

//[[Rcpp::export]]
CharacterMatrix get_pedigree_edgelist(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector< std::pair<Individual*, Individual*>* >* rels = p->get_relations();
  
  CharacterMatrix edgelist(rels->size(), 2);
  int i = 0;
  
  for (auto pair: *rels) {
    edgelist(i, 0) = std::to_string(pair->first->get_pid());
    edgelist(i, 1) = std::to_string(pair->second->get_pid());
    ++i;
  }
  
  return edgelist;
}

/*
//[[Rcpp::export]]
List get_pedigree_as_graph(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  CharacterVector nodes(inds->size());
  List edges(inds->size());
  
  int i = 0;
  for (auto individual : *inds) {
    nodes(i) = std::to_string(individual->get_pid());
    
    std::vector<Individual*>* children = individual->get_children();
    
    CharacterVector edges_i(children->size());
    int j = 0;
    for (auto child : *children) {
      edges_i(j) = std::to_string(child->get_pid());
      ++j;
    }
    
    List tmp;
    tmp["edges"] = edges_i;    
    edges(i) = tmp; 
   
    ++i;
  }
  
  edges.attr("names") = nodes;
  
  List ret;
  ret["nodes"] = nodes;
  ret["edgelist"] = edges;
  
  return ret;
}
*/

//' Get pedigree information as graph (mainly intended for plotting)
//' 
//' @export
// [[Rcpp::export]]
List get_pedigree_as_graph(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  
  CharacterVector nodes(inds->size());
  
  int i = 0;
  for (auto individual : *inds) {
    nodes(i) = std::to_string(individual->get_pid());   
    ++i;
  }
  
  List ret;
  ret["nodes"] = nodes;
  ret["edgelist"] = get_pedigree_edgelist(ped);
  
  return ret;
}

//[[Rcpp::export]]
int meiosis_dist_tree(Rcpp::XPtr<Individual> src, Rcpp::XPtr<Individual> dest) {
  Individual* i1 = src;
  Individual* i2 = dest;
  int dist = i1->meiosis_dist_tree(i2);
  
  if (dist == -1) {
    Rcpp::stop("Individuals were not in the same pedigree");
  }
  
  return dist;
}

//[[Rcpp::export]]
IntegerMatrix meiosis_dist_tree_matrix(Rcpp::XPtr<Pedigree> ped) {  
  Pedigree* p = ped;
  
  std::vector<Individual*>* inds = p->get_all_individuals();
  size_t n = inds->size();
  
  CharacterVector nms(n);
  IntegerMatrix res(n, n);
  int i = 0;
  
  for (size_t i = 0; i < (n-1); ++i) {
    Individual* i1 = inds->at(i);
    nms[i] = std::to_string(i1->get_pid());
    res(i,i) = 0;
    
    for (size_t j = i+1; j < n; ++j) {
      Individual* i2 = inds->at(j);
      
      int dist = i1->meiosis_dist_tree(i2);
      res(i,j) = dist;
      res(j,i) = dist;
    }
  }
  
  nms[n-1] = std::to_string(inds->at(n-1)->get_pid());
  
  rownames(res) = nms;
  colnames(res) = nms;
  
  return res;
}

//' Get pedigree from individual
//' 
//' @export
//[[Rcpp::export]]
Rcpp::XPtr<Pedigree> get_pedigree_from_individual(Rcpp::XPtr<Individual> individual) {  
  Individual* i = individual;  
  Rcpp::XPtr<Pedigree> res(i->get_pedigree(), false); // do NOT delete pedigree when not used any more, it still exists in list of pedigrees etc.!
  res.attr("class") = CharacterVector::create("popr_pedigree", "externalptr");
  
  return res;
}

//' Get pedigree id from pid
//' 
//' @export
// [[Rcpp::export]]
IntegerVector get_pedigree_id_from_pid(Rcpp::XPtr<Population> population, IntegerVector pids) {  
  std::unordered_map<int, Individual*> pop = *(population->get_population());
  
  int N = pids.size();
  IntegerVector pedigree_ids(N);
  
  for (int i = 0; i < N; ++i) {
    Individual* ind = population->get_individual(pids[i]);
    pedigree_ids[i] = ind->get_pedigree_id();
  }
  
  return pedigree_ids;
}



/*


int get_total_brothers_count() {
  double u = unif_rand();
  if      (u <= 0.703286965481710) { return 1; }
  else if (u <= 0.956079340546110) { return 2; }
  else if (u <= 0.995327150953117) { return 3; }
  else if (u <= 0.999448737764039) { return 4; }
  else if (u <= 0.999891808346213) { return 5; }
  else if (u <= 0.999989696032973) { return 6; }
  else if (u <= 0.999994848016486) { return 7; }
  return 8;
  //else if (u <= 1) { return 8; }
}

int get_brother_age_diff() {
  double u = unif_rand();
  if     (u <= 0.0432283523192614) { return 0; }
  else if (u <= 0.119201846474574) { return 1; }
  else if (u <= 0.375429975429975) { return 2; }
  else if (u <= 0.633921524830616) { return 3; }
  else if (u <= 0.770188370188370) { return 4; }
  else if (u <= 0.847859429677612) { return 5; }
  else if (u <= 0.898548134911771) { return 6; }
  else if (u <= 0.931650658923386) { return 7; }
  else if (u <= 0.954627354627355) { return 8; }
  else if (u <= 0.969175787357606) { return 9; }
  else if (u <= 0.979390961209143) { return 10; }
  else if (u <= 0.987536296627206) { return 11; }
  else if (u <= 0.992673665400938) { return 12; }
  else if (u <= 0.995756086665178) { return 13; }
  else if (u <= 0.997617452162907) { return 14; }
  else if (u <= 0.998555580373762) { return 15; }
  else if (u <= 0.999404363040727) { return 16; }
  else if (u <= 0.999776636140272) { return 17; }
  else if (u <= 0.999925545380091) { return 18; }
  else if (u <= 0.999970218152036) { return 19; }
  else if (u <= 0.999985109076018) { return 20; }
  return 21;
  //else if (u <= 1) { return 21; }
}

int get_age_first_boy() {
  double u = unif_rand();
       if (u <= 2.13125039128425e-05) { return 12; }
  else if (u <= 9.79043148496203e-05) { return 13; }
  else if (u <= 0.000237767621777649) { return 14; }
  else if (u <= 0.000491519621489931) { return 15; }
  else if (u <= 0.000989699400452624) { return 16; }
  else if (u <= 0.00227377776120139) { return 17; }
  else if (u <= 0.00551527640319528) { return 18; }
  else if (u <= 0.0129446820640627) { return 19; }
  else if (u <= 0.0272407100793491) { return 20; }
  else if (u <= 0.0496221692665701) { return 21; }
  else if (u <= 0.0822622690090885) { return 22; }
  else if (u <= 0.12520363431473) { return 23; }
  else if (u <= 0.17728073762576) { return 24; }
  else if (u <= 0.237352693967096) { return 25; }
  else if (u <= 0.304276620316411) { return 26; }
  else if (u <= 0.374604553150055) { return 27; }
  else if (u <= 0.446408376879663) { return 28; }
  else if (u <= 0.517004048043712) { return 29; }
  else if (u <= 0.585188076186873) { return 30; }
  else if (u <= 0.648707330036111) { return 31; }
  else if (u <= 0.705415240837954) { return 32; }
  else if (u <= 0.7552512011594) { return 33; }
  else if (u <= 0.798089334024214) { return 34; }
  else if (u <= 0.83482809467547) { return 35; }
  else if (u <= 0.865947680466957) { return 36; }
  else if (u <= 0.891667876595274) { return 37; }
  else if (u <= 0.912807882429572) { return 38; }
  else if (u <= 0.930364057527776) { return 39; }
  else if (u <= 0.944628782802941) { return 40; }
  else if (u <= 0.956357320112477) { return 41; }
  else if (u <= 0.965597622590189) { return 42; }
  else if (u <= 0.972869182518951) { return 43; }
  else if (u <= 0.978710806638312) { return 44; }
  else if (u <= 0.983416207892819) { return 45; }
  else if (u <= 0.987069970282377) { return 46; }
  else if (u <= 0.989872564546916) { return 47; }
  else if (u <= 0.992097057142819) { return 48; }
  else if (u <= 0.993836690274705) { return 49; }
  else if (u <= 0.995193364351907) { return 50; }
  else if (u <= 0.996297618460891) { return 51; }
  else if (u <= 0.997143458459932) { return 52; }
  else if (u <= 0.997788827719043) { return 53; }
  else if (u <= 0.998291003592489) { return 54; }
  else if (u <= 0.998668634521195) { return 55; }
  else if (u <= 0.998965011528733) { return 56; }
  else if (u <= 0.999181466646597) { return 57; }
  else if (u <= 0.999370615118824) { return 58; }
  else if (u <= 0.999516472567477) { return 59; }
  else if (u <= 0.999632359307503) { return 60; }
  else if (u <= 0.99971627729166) { return 61; }
  else if (u <= 0.999793535118344) { return 62; }
  else if (u <= 0.999844818330885) { return 63; }
  else if (u <= 0.999894103496183) { return 64; }
  else if (u <= 0.999920744126074) { return 65; }
  else if (u <= 0.99994871678746) { return 66; }
  else if (u <= 0.999968031244131) { return 67; }
  else if (u <= 0.999982683590571) { return 68; }
  else if (u <= 0.999994005858274) { return 69; }
  return 70;
  //else if (u <= 1) { return 70; }
}

//FIXME: Scheme with e.g. 10 brothers constantly

//' Extend pedigrees by random sampling
//'
//' @param pid Vector of pids with pid_f == NA
//'
//' @export
// [[Rcpp::export]]
List extend_male_pedigrees_one_generation(IntegerVector pid, IntegerVector birthyear, NumericVector etrs89e, NumericVector etrs89n, bool progress = true) {
  size_t N = pid.size();
  
  IntegerVector pid_f(N);
  
  if (birthyear.size() != N) stop("birthyear.size() != N = pid.size()");
  if (etrs89e.size() != N) stop("etrs89e.size() != N = pid.size()");
  if (etrs89n.size() != N) stop("etrs89n.size() != N = pid.size()");
  if (regions.size() != N) stop("regions.size() != N = pid.size()");
  
  LogicalVector f_rnd(N);  
  
  std::vector<int> pids_no_father;
  pids_no_father.reserve(N); // worst case
  
  std::unordered_map<int, int> pid_index;
  
  for (size_t i = 0; i < N; ++i) {
    pid_index[pid[i]] = i;
    pids_no_father.push_back(pid[i]);
  }
  
  IntegerVector new_fathers_pids;
  IntegerVector new_fathers_birthyears; // FIXME: Use?
  
  for (size_t i = 0; i < 10; ++i) {
    int oldest_brother_pids_no_father_index = floor(unif_rand()*pids_no_father.size());
    int oldest_brother_pid = pids_no_father[oldest_brother_pids_no_father_index];
    int oldest_brother_index = pid_index[oldest_brother_pid]; // index in pids_no_father not same as index in pid because pids_no_father is getting smaller
    int oldest_brother_birthyear = birthyear[oldest_brother_index];
    
    //> formatC(range(db_males$pid), big.mark = ",")
    //[1] "1,900,006"  "21,023,836"
    
    int father_pid = 99000000 + i;    
    int father_birthyear = oldest_brother_birthyear - get_age_first_boy();
    
    pid_f[oldest_brother_index] = father_pid;

    new_fathers_pids.push_back(father_pid);
    new_fathers_birthyears.push_back(father_birthyear);
    
    // FIXME: Add new_father to pool of individuals without father?
        
    // Number of brothers total:
    int n_brothers = get_total_brothers_count();
    
    if (n_brothers == 1) {
      // One brother only
      pids_no_father.erase(oldest_brother_pids_no_father_index); // FIXME: Complexity??? O(1 + n_elements_efter_index)?
      continue;
    }
    
    std::vector<int> brothers_candidates(n_brothers);
    std::vector<int> brothers_pid(n_brothers);
    
    for (size_t j = 0; j < n_brothers; ++j) {
      brothers_pid.push_back();
    }
    
    size_t index = floor(unif_rand()*pids_no_father.size());
    
    if (IntegerVector::is_na(pid_f[i])) {
      pids_no_father.push_back(pid[i]);
    }
  }
  
  List res;
  res['new_fathers_pid'] = new_fathers_pids;
  res['new_fathers_birthyears'] = new_fathers_birthyears;
  res['pid_f'] = pid_f;

}


*/








//' @export
// [[Rcpp::export]]
void pedigree_populate_father_haplotypes(Rcpp::XPtr<Pedigree> ped, int loci, double mutation_rate) {  
  Pedigree* p = ped;
  ped->populate_father_haplotypes(loci, mutation_rate);
}

//' @export
// [[Rcpp::export]]
void pedigrees_all_populate_father_haplotypes(Rcpp::XPtr< std::vector<Pedigree*> > pedigrees, int loci, double mutation_rate, bool progress = true) {
  std::vector<Pedigree*> peds = (*pedigrees);
  
  size_t N = peds.size();
  Progress p(N, progress);
  
  for (size_t i = 0; i < N; ++i) {
    peds.at(i)->populate_father_haplotypes(loci, mutation_rate);
    
     if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort()) {
      stop("Aborted.");
    }
    
    if (progress) {
      p.increment();
    }
  }
}


//' @export
// [[Rcpp::export]]
List pedigree_get_father_haplotypes_pids(Rcpp::XPtr<Population> population, IntegerVector pids) {  
 
  size_t N = pids.size();
  List haps(N);
  
  for (size_t i = 0; i < N; ++i) {
    Individual* indv = population->get_individual(pids[i]);
    haps(i) = indv->get_father_haplotype();
  }
  
  return haps;
}

/*
//' @export
// [[Rcpp::export]]
List pedigree_get_father_haplotypes_pedigree(Rcpp::XPtr<Pedigree> ped) {  
  std::vector<Individual*> inds = *(ped->get_all_individuals());
  size_t N = inds.size();
  List haps(N);
  
  for (size_t i = 0; i < N; ++i) {
    Individual* indv = inds[i];
    haps(i) = indv->get_father_haplotype();
  }
  
  return haps;
}
*/



