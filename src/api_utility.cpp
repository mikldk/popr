#include <Rcpp.h>

// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>

#include <string>

#include "popr_types.hpp"

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
    
    if (i % CHECK_ABORT_EVERY == 0 && Progress::check_abort() ) {
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

//[[Rcpp::export]]
void popr_test() {
  Rcout << "mikl was here" << std::endl;
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

//' Build pedigrees
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

