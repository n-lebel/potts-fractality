#include "CMGrid.h"

std::vector<int> color(int a) {

	std::vector<int> rgb;

	while(a > 8) {
		a = a - 8;
	}
	switch(a) {
	case 1: rgb = {255, 255, 0}; break;
	case 2:  rgb = {255, 0, 255}; break;
	case 3:  rgb = {0, 255, 0}; break;
	case 4:  rgb = {0, 0, 255}; break;
	case 5:  rgb = {255, 0, 255}; break;
	case 6:  rgb = {255, 0, 0}; break;
	case 7:  rgb = {255, 255, 255}; break;
	case 8:  rgb = {0, 0, 0}; break;
	default: rgb = {0, 255, 255}; break;
	}
	return rgb;
}

CMGrid::CMGrid(int gs, float number_of_states, float fullness): grid_size(gs), q(number_of_states), fullness(fullness) {

		std::srand(time(NULL));     // setting RNG seed to current time
		std::random_device rd;     // only used once to initialise (seed) engine
		std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
		std::bernoulli_distribution distr_full(fullness);

    for (float i = 0; i < grid_size; i++) {
        for (float j = 0; j < grid_size; j++) {
            Key k(i, j);
            key_list.push_back(k);

            Bond b1(k, Key(i+1, j)), b2(k, Key(i, j+1));

            if (i != grid_size-1) {
                bond_parameters.insert(std::make_pair(b1, distr_full(rng)));
                bond_list.push_back(b1);
            }
            if (j != grid_size-1) {
                bond_parameters.insert(std::make_pair(b2, distr_full(rng)));
                bond_list.push_back(b2);
            }

            std::vector<Key> nb;
            if (i != 0) nb.push_back(Key(i-1, j));
            if (i != grid_size-1) nb.push_back(Key(i+1, j));
            if (j != 0) nb.push_back(Key(i, j-1));
            if (j != grid_size-1) nb.push_back(Key(i, j+1));

            neighbours_list.insert(std::make_pair(k, nb));
        }
    }

		nb_links = bond_list.size();
		crit_prob = sqrt(q)/(1+sqrt(q));
		pq = crit_prob/q;

		for (auto const& l: key_list) {
			if (l.i == 0 || l.i == grid_size-1 || l.j == 0 || l.j == grid_size-1) sides[l] = true;
			else sides[l] = false;
		}
}

Bond CMGrid::make_valid_bond(Key k1, Key k2) {
	if (k1.j < k2.j) return Bond(k1, k2);
	else if (k1.j == k2.j && k1.i < k2.i) return Bond(k1, k2);
	return Bond(k2, k1);
}

/*Bond CMGrid::make_valid_bond(Key k1, Key k2) {
	if (k1.j < k2.j) return Bond(k1, k2);
	return Bond(k2, k2);
}*/


std::vector<Key> CMGrid::bonded_neighbours(Key k) {
	std::vector<Key> nb = neighbours_list.at(k);
	std::vector<Key> final;
	for (auto const& x: nb) {
		if (bonded(k, x)) final.push_back(x);
	}
	return final;
}

std::vector<Key> CMGrid::unbonded_neighbours(Key k) {
	std::vector<Key> nb = neighbours_list.at(k);
	std::vector<Key> final;
	for (auto const& x: nb) {
		if (!bonded(k, x)) final.push_back(x);
	}
	return final;
}

bool CMGrid::bonded(Key k1, Key k2) {
	Bond b = make_valid_bond(k1, k2);
	// is the initial bond in the list?
	if (bond_parameters.at(b) == 1) return true;
	return false;
}

float CMGrid::bond_density() {
	float n = 0;
	for (auto const& b: bond_list) {
		if (bond_parameters.at(b) == 1) n++;
	}
	return n/nb_links;
}

float CMGrid::largest_cluster_fraction() {
	map_clusters();
	float largest_size = -1;
	for (auto const& [index, list]: cluster_list) {
		if (list.size() > largest_size) largest_size = list.size();
	}

	return largest_size/(grid_size*grid_size);
}

std::vector<Key> CMGrid::largest_cluster() {
	int biggest_cluster_index = 0;
	size_t biggest_cluster_size = 0;
	for (auto const& [index, cluster]: cluster_list) {
		if (cluster.size() > biggest_cluster_size) {
			biggest_cluster_index = index;
			biggest_cluster_size = cluster.size();
		}
	}
	return cluster_list[biggest_cluster_index];
}


int CMGrid::multi_bernoulli(double prob, int m) {
	// this function will give a random number between 1 and m with probability prob
	// or will return 0 with probability 1-m*prob

	std::srand(time(NULL));     // setting RNG seed to current time
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_int_distribution distr_active(1, m);
	std::bernoulli_distribution distr_inactive(1-static_cast<double>(m)/q);

	if (distr_inactive(rng)) return 0;
	else return distr_active(rng);
}

void CMGrid::map_clusters() {

	// reset cluster map
	for (auto const& k: key_list)	clusters[k] = 0;
	cluster_list.clear();


	// use BFS on each key to determine the cluster it belongs to
	int count = 0;
	for (auto const& k: key_list) {
		// check that node hasn't been assigned a cluster yet before BFS
		if (clusters[k] == 0) {
			count++; // means it's a new cluster
			std::queue<Key> q;
			q.push(k);

			while (!q.empty()) {
				Key vertex = q.front();
				q.pop();
				if (clusters[vertex] == 0) {
					clusters[vertex] = count;
					for (auto const& neighbour: bonded_neighbours(vertex)) {
						if (clusters[neighbour] == 0) { q.push(neighbour); }
					}
				}
			}
		}
	}

	for (auto const& k: key_list) {
		int index = clusters[k];
		cluster_list[index].push_back(k);
	}
}

void CMGrid::evolve(int nb_steps, float p) {

	std::srand(time(NULL));     // setting RNG seed to current time
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_real_distribution distr_active(0.0, 1.0);
	std::bernoulli_distribution distr_pc(p);

	int count = 0;

	for (int i = 0; i < nb_steps; i++) {
		float bd = bond_density();
		std::cout << i << ": " << bd << std::endl << std::flush;
		active_sites.clear();

		map_clusters();
		for (auto const& [index, list]: cluster_list) {
			int activity = multi_bernoulli(1/q, floor(q));
			//std::cout << activity << " ";
			for (auto const& k: list) active_sites[k] = activity;
		}
		for (auto const& b: bond_list) {
			if (active_sites[b.a] == active_sites[b.b] && active_sites[b.a] != 0) bond_parameters.at(b) = 0;
		}
		for (auto const& b: bond_list) {
			if (active_sites[b.a] == active_sites[b.b] && active_sites[b.a] != 0) {
				if (distr_pc(rng) == true) bond_parameters.at(b) = 1;
			}
		}

		if (i > nb_steps-floor(nb_steps/15)+2) {
			count++;
			final_averaged_density += bd;
		}
	}
	final_averaged_density = final_averaged_density/count;
}

float CMGrid::fixed_interface_step(float p) {

	std::srand(time(NULL));     // setting RNG seed to current time
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_real_distribution distr_active(0.0, 1.0);
	std::bernoulli_distribution distr_pc(p);

	active_sites.clear();
	map_clusters();

	int left_cluster_index, right_cluster_index;
	std::vector<Key> left_cluster, right_cluster;

	// first take care of the side clusters (only spinwise, except for the sides to which we always add bonds)
	left_cluster_index = clusters[Key(0, 0)], right_cluster_index = clusters[Key(0, grid_size-1)];
	left_cluster = cluster_list[left_cluster_index], right_cluster = cluster_list[right_cluster_index];

	int leftspin = multi_bernoulli(1/q, floor(q)), rightspin = multi_bernoulli(1/q, floor(q));
	while (leftspin == rightspin) rightspin = multi_bernoulli(1/q, floor(q));

	for (auto const& k: left_cluster) active_sites[k] = leftspin;
	for (auto const& k: right_cluster) active_sites[k] = rightspin;

	// then take care of the rest
	for (auto const& [index, list]: cluster_list) if (index != left_cluster_index && index != right_cluster_index){
		int activity = multi_bernoulli(1/q, floor(q));
		for (auto const& k: list) active_sites[k] = activity;
	}

	// take care of the bonds
	for (auto const& b: bond_list) {
		if (active_sites[b.a] == active_sites[b.b] && active_sites[b.a] != 0) bond_parameters.at(b) = 0;
	}
	for (auto const& b: bond_list) {
		if (active_sites[b.a] == active_sites[b.b] && active_sites[b.a] != 0) {
			if (distr_pc(rng) == true) bond_parameters.at(b) = 1;
		}
	}
	for (auto const& b: bond_list) {
		if (sides[b.a] && sides[b.b]) if (active_sites[b.b] == active_sites[b.a]) bond_parameters.at(b) = true;
	}
	return bond_density();
}

void CMGrid::fixed_interface_evolve(int nb_steps, float p, bool progress) {
	std::srand(time(NULL));     // setting RNG seed to current time
	std::random_device rd;     // only used once to initialise (seed) engine
	std::mt19937 rng(rd());    // random-number engine used (Mersenne-Twister in this case)
	std::uniform_real_distribution distr_active(0.0, 1.0);
	std::bernoulli_distribution distr_pc(p);

	// create the two fixed boundary conditions ONLY WORKS FOR INTEGER q
	int leftspin = multi_bernoulli(1/q, floor(q)), rightspin = multi_bernoulli(1/q, floor(q));
	while (leftspin == rightspin) rightspin = multi_bernoulli(1/q, floor(q));
	//std::cout << "Left spin:" << leftspin << ", Right spin: " << rightspin;

	for (auto const& [k, b]: sides) if (b == true) {
		if (k.j < floor(grid_size/2)) {
			active_sites[k] = leftspin;
			for (auto const& n: neighbours_list[k]) if (sides[n] == true && active_sites[k] == active_sites[n]) {
				bond_parameters.at(make_valid_bond(k, n)) = true;
			}
		}
		else {
			active_sites[k] = rightspin;
			for (auto const& n: neighbours_list[k]) if (sides[n] == true && active_sites[k] == active_sites[n]) {
				bond_parameters.at(make_valid_bond(k, n)) = true;
			}
		}
	}

	// truncates the grid, removes all bonds in the middle
	for (auto const& b: bond_list) {
		if (b.a.j == floor(grid_size/2)-1 && b.b.j == floor(grid_size/2)) bond_parameters.at(b) = 0;
	}

	int count = 0;

	time_t tstart, tend;
	tstart = time(0);

	for (int i = 0; i < nb_steps; i++) {
		float bd = fixed_interface_step(p);
		if (progress) std::cout << i << ": " << bd << std::endl;

		if (i > nb_steps-floor(nb_steps/15)+2) {
			count++;
			final_averaged_density += bd;
		}
	}
	map_clusters();
	final_averaged_density = final_averaged_density/count;
	tend = time(0);
	std::cout << "It took " << difftime(tend, tstart) << " to thermalize the grid." << std::endl;
	evolve_steps += nb_steps;
}


/*std::map<int, int> CMGrid::frac_sweep() {
	map_clusters();

	int biggest_cluster_index = 0;
	size_t biggest_cluster_size = 0;
	for (auto const& cluster: cluster_list) {
		if (cluster.second.size() > biggest_cluster_size) {
			biggest_cluster_index = cluster.first;
			biggest_cluster_size = cluster.second.size();
		}
	}

	std::vector<Key> perc_cluster = cluster_list[biggest_cluster_index];
	lc_size = biggest_cluster_size;

	// approximate the center of the cluster
	int leftmost(perc_cluster[0].j), rightmost(perc_cluster[0].j), upmost(perc_cluster[0].i), downmost(perc_cluster[0].i);
	for(auto const& k: perc_cluster) {
		if (k.i > downmost) downmost = k.i;
		if (k.i < upmost) upmost = k.i;
		if (k.j > rightmost) rightmost = k.j;
		if (k.j < leftmost) leftmost = k.j;
	}
	Key center(floor((downmost+upmost)/2), floor((leftmost+rightmost)/2));
	// checking if center belongs to the cluster
	bool belongs = false;
	for (auto const& k: perc_cluster) if (k == center) belongs = true;
	std::cout << "Center(p): " << center << std::endl;

	std::cout << belongs << std::endl;
	// if it doesn't, circle around it until you find one that does
	float a = 0, b = 0;
	if (belongs == false) {
		for (float r = 1; r < 5; r++) {
			for (float u = center.i-r; u <= center.i+r; u++) {
				for (float v = center.j-r; v <= center.j+r; v++) {
					for (auto const& k: perc_cluster) {
						if (k.i == u && k.j == v) {
							a = u;
							b = v;
							goto label; // only used to get out of the 3 nested loops
						}
					}
				}
			}
		}
label:          center.i = a;
		center.j = b;
	}

	std::map<int, int> cluster_chars;
	std::cout << "First element: " << perc_cluster[0] << ", nb of nodes " << perc_cluster.size() << std::endl;
	std::cout << "Center: " << center << std::endl;
	std::cout << "Upmost: " << upmost << "; downmost: " << downmost << "; leftmost: " << leftmost << "; rightmost: " << rightmost << std::endl;
	for (int r = 0; r < grid_size; r++) {
		cluster_chars.insert(std::pair<int, int>(r, 0));
		for (int u = center.i-r; u <= center.i+r; u++) {
			for (int v = center.j-r; v <= center.j+r; v++) {
				if (u < grid_size && v < grid_size && u >= 0 && v >= 0) {
					Key k(u, v);

					try {
						if (clusters.at(k) == biggest_cluster_index) cluster_chars[r]++;
					}
					catch (...) { std::cout << k << std::endl; }
				}
			}
		}
		if (r>1 && cluster_chars[r] == cluster_chars[r-1]) {
			cluster_chars.erase(r);
			break;
			// if the next radius box gives the same number, means the cluster's over so we stop
		}
		std::cout << "Radius: " << r << "; number of nodes: " << cluster_chars[r] << std::endl;
	}
	return cluster_chars;
}*/

void CMGrid::display_clusters(std::string name) {

	map_clusters();
	// Little function using cairo to display the graph
	int WIDTH = 4096;
	int HEIGHT = 4096;

	int u = HEIGHT/grid_size;
	int v = WIDTH/grid_size;

	cairo_surface_t *surface;
	cairo_t *cr;

	surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, HEIGHT, WIDTH);
	cr = cairo_create(surface);

	//cairo_set_source_rgb(cr, 0, 255, 255);
	cairo_set_line_width(cr, 4);

	for (auto const& x : bond_parameters)
	{
		if (x.second == 1) {
			Key k1 = x.first.a;
			Key k2 = x.first.b;
			int c = clusters.at(k1);

			std::vector<int> rgb(color(c));
			//std::cout << "Cluster " << c << " color " << rgb[0] << ", " << rgb[1] << ", " << rgb[2] << std::endl;
			cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
			cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
			cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
			cairo_stroke(cr);
		}
	}


		cairo_surface_write_to_png(surface, name.c_str()); // Output to PNG
		cairo_destroy(cr);
		cairo_surface_destroy(surface);
}

void CMGrid::display_left_cluster(std::string name) {
	int WIDTH = 4096;
	int HEIGHT = 4096;

	int u = HEIGHT/grid_size;
	int v = WIDTH/grid_size;

	cairo_surface_t *surface;
	cairo_t *cr;

	surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, HEIGHT, WIDTH);
	cr = cairo_create(surface);

	//cairo_set_source_rgb(cr, 0, 255, 255);
	cairo_set_line_width(cr, 4);
	int left_cluster_index = clusters[Key(0, 0)];
	std::vector<Key> left_cluster = cluster_list[left_cluster_index];

	for (auto const& k1 : left_cluster)
	{
		for (auto const& k2: neighbours_list[k1]) {
			if (bond_parameters.at(make_valid_bond(k1, k2)) == 1 ) {
				int c = clusters.at(k1);

				std::vector<int> rgb(color(c));
				//std::cout << "Cluster " << c << " color " << rgb[0] << ", " << rgb[1] << ", " << rgb[2] << std::endl;
				cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
				cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
				cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
				cairo_stroke(cr);
			}
		}
	}




		cairo_surface_write_to_png(surface, name.c_str()); // Output to PNG
		cairo_destroy(cr);
		cairo_surface_destroy(surface);
}

/*
void CMGrid::display_largest_cluster_hull(std::string name) {

	std::vector<Bond> hull = determine_hull();

	int biggest_cluster_index = 0;
	size_t biggest_cluster_size = 0;
	for (auto const& [index, cluster]: cluster_list) {
		if (cluster.size() > biggest_cluster_size) {
			biggest_cluster_index = index;
			biggest_cluster_size = cluster.size();
		}
	}
	std::vector<Key> perc_cluster = cluster_list[biggest_cluster_index];
	std::vector<Bond> cluster;
	for (auto const& b: bond_list) {
		if (bond_parameters[b] == 1 && clusters[b.a] == biggest_cluster_index && clusters[b.b] == biggest_cluster_index) {
			cluster.push_back(b);
		}
	}

	int WIDTH = 4096;
	int HEIGHT = 4096;

	int u = HEIGHT/grid_size;
	int v = WIDTH/grid_size;

	cairo_surface_t *surface;
	cairo_t *cr;

	surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, HEIGHT, WIDTH);
	cr = cairo_create(surface);

	cairo_set_line_width(cr, 20);

	for (auto const& b: cluster) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color(6));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	cairo_set_line_width(cr, 2);

	for (auto const& b: hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color(7));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}


		cairo_surface_write_to_png(surface, name.c_str()); // Output to PNG
		cairo_destroy(cr);
		cairo_surface_destroy(surface);
}
*/
