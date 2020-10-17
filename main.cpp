#include "utilities.h"

int main(int argc, char* argv[]) {

	// Command: ./main <q> <grid size> <steps>
	float f = 0.5, q = atoi(argv[1]), gs = atoi(argv[2]);
	int steps = atoi(argv[3]);

  CMGrid a(gs, q, f);
  a.fixed_interface_evolve(steps, a.crit_prob);

  a.display_clusters("images/unlinked_clusters.png");
  int left_cluster_index = a.clusters[Key(0, 0)], right_cluster_index = a.clusters[Key(0, a.grid_size-1)];
  std::vector<Key> left_cluster = a.cluster_list[left_cluster_index], right_cluster = a.cluster_list[right_cluster_index];


  Hull hl(a, left_cluster);

	hl.display_linked_hull("pres_images/left_hull_"+ std::to_string(static_cast<int>(gs)) +".png", 0);
	hl.display_linked_ap("pres_images/left_ap_"+ std::to_string(static_cast<int>(gs)) +".png", 0);

	Hull hr(a, right_cluster);

	hr.display_linked_hull("pres_images/right_hull_"+ std::to_string(static_cast<int>(gs)) +".png", 0);
	hr.display_linked_ap("pres_images/right_ap_"+ std::to_string(static_cast<int>(gs)) +".png", 0);

  return 0;
}
