#ifndef mcgrid_h
#define mcgrid_h

#include "key.h"
#include "bond.h"

#include <cairo/cairo.h>
#include <iostream>
#include <map>
#include <cmath>
#include <vector>
#include <stack>
#include <string>
#include <random>
#include <queue>

// Chayes-Machta Grid
class CMGrid {
public:
  float grid_size;
  float q;
  float pq;
  float fullness;
  float crit_prob;
  float nb_links;
  float final_averaged_density;
  int lc_size;
  int evolve_steps = 0;

  std::map<Bond, bool> bond_parameters;
  std::map<Key, int> active_sites;

  std::vector<Bond> bond_list;
  std::vector<Key> key_list;
  std::map<Key, std::vector<Key>> neighbours_list;
  std::map<Key, bool> sides;

  std::map<Key, int> clusters;
  std::map<int, std::vector<Key>> cluster_list;

  CMGrid(int gs, float number_of_states, float fullness);
  Bond make_valid_bond(Key k1, Key k2);
  bool bonded(Key k1, Key k2);
  std::vector<Key> bonded_neighbours(Key k);
  std::vector<Key> unbonded_neighbours(Key k);
  float bond_density();
  float largest_cluster_fraction();
  std::vector<Key> largest_cluster();
  int multi_bernoulli(double prob, int m);
  void map_clusters();
  void evolve(int nb_steps, float p);
  float fixed_interface_step(float p);
  void fixed_interface_evolve(int nb_steps, float p, bool progress=true);
  std::map<int, int> frac_sweep();
  void display_clusters(std::string name);
  bool bernoulli(const double &p);

  void display_left_cluster(std::string name);

};

#endif
