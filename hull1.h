#ifndef hull_h
#define hull_h

#include "CMGrid.h"
#include <algorithm>
#include <list>
#include <cmath>

/*
The Hull class takes a CMGrid in argument and creates a parallel lattice in order to
allow correct pathfinding. We differentiate the outer lattice (i.e. the set of original
lattice sites) and the inner lattice (which is a 0.5-offset parallel lattice).
*/

class Hull {
public:
  CMGrid g;
  float grid_size;

  std::vector<Key> fine_grid;
  std::map<Key, bool> fine_sides;
  std::map<Key, std::vector<Key>> fine_neighbours_list;
  std::map<Key, std::vector<Key>> outer_neighbours_list;

  std::map<Key, bool> cluster_map;
  std::map<Key, bool> dense_cluster_map;
  std::map<Bond, bool> dense_cluster_bond_parameters;

  // the white scaffolding
  std::vector<Bond> hull;
  std::vector<Bond> ap;

  // the individual points of the perimeters
  std::map<Key, bool> hull_points;
  std::map<Key, bool> ap_points;

  // the bonds that connect the points to each other
  std::vector<Bond> linked_ap;
  std::vector<Bond> linked_hull;

  // maps that assign to each perimeter point its index from the top
  std::map<Key, float> ordered_ap_keys;
  std::map<Key, float> ordered_hull_keys;

  std::vector<Key> ordered_hull, ordered_ap;

  Hull(CMGrid grid, std::vector<Key> c);
  Bond make_valid_bond(Key k1, Key k2);
  Bond make_valid_half_bond(Key k1, Key k2);
  std::vector<Key> square(Key k, float radius);
  float bond_length(Bond b);
  float eucl_dist(Key k1, Key k2);
  float angle(Key k1, Key k2, bool counterclockwise);
  bool parallel_aligned(Bond b1, Bond b2);
  bool perpendicular(Bond b1, Bond b2);
  bool on_segment(Bond b, Key k);
  int crossing_link(Bond b, Key center, float radius);
  int sign(float n);
  Key crossing_point(Bond b, Key center, float r);
  int intercrossing_link(Bond b, Key center, float r);
  Key intercrossing_point(Bond b, Key center, float r);

  // returns true if v1 and v2 have at least one common element
  bool common_element(std::vector<Key> v1, std::vector<Key> s2);

  /* inbetween_key() determines whether there exists a key in _list_ that lies
  between k1 and k2 (i.e. on the segment k1-k2) */
  bool inbetween_key(Key k1, Key k2, std::map<Key, bool> list);
  /* permitted_move() determines whether a move from one point of the inner lattice
   to another is blocked by a bond in the outer lattice */
  bool permitted_move(Key start, Key end);
  bool dense_permitted_move(Key start, Key end);

  std::vector<Bond> determine_hull();
  std::vector<Bond> new_determine_hull();
  std::vector<Bond> determine_ap();

  std::vector<Bond> link_ap();
  std::vector<Bond> old_link_ap();
  std::vector<Bond> link_hull();
  std::vector<Bond> old_link_hull();

  void write_hull();
  void write_ap();

  float hull_length();
  float hull_width();
  float ap_length();
  float ap_width();
  std::map<int, int> fractal_sweep(std::vector<Bond> h_or_ap);
  // bool h_or_ap: 0 for hull, 1 for accessible perimeter
  std::vector<std::map<float, float>> yardstick(std::vector<float> radii);
  std::vector<std::map<double, double>> yardstick_ap(std::vector<double> radii);
  std::vector<std::map<double, double>> yardstick_new(std::vector<double> radii);
  std::vector<Key> test_yardstick_h(float r);
  std::vector<Key> test_yardstick_interh(float r);
  std::vector<Key> test_yardstick_ap(float r);


  void display_hull(std::string name);
  void display_ap(std::string name);
  void display_linked_ap(std::string name, bool cluster);
  void display_linked_hull(std::string name, bool cluster);
  void display_yardstick_h(std::string name, float radius);
  void display_yardstick_ap(std::string name, float radius);

  std::map<double, double> yardstick_rilder(std::vector<double> radii);
  std::vector<Key> yardstick_rilder_test(double stick);
  void display_yardstick_rilder(std::string name, float radius);

};

#endif
