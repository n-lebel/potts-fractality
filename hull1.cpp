#include "hull1.h"
#define LEFT 1
#define RIGHT 0

std::string gen_name() {
	std::srand(time(NULL)); // setting RNG seed to current time
	std::random_device rd; // only used once to initialise (seed) engine
	std::mt19937 rng(rd()); // random-number engine used (Mersenne-Twister in this case)

	std::uniform_int_distribution<int> uni(0, 1000000000);
	int t = uni(rng);
	return std::to_string(t);
}

double distance(double x1, double y1, double x2, double y2){
  return sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2));
};

std::vector<int> color_hull(int a) {

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

Hull::Hull(CMGrid grid, std::vector<Key> c) : g(grid), grid_size(grid.grid_size) {
	for (float i = 0; i < grid_size; i += 0.5) {
		for (float j = 0; j < grid_size; j += 0.5) {
			Key k(i, j);
			fine_grid.push_back(k);

			fine_sides[k] = (i == 0 || i == grid_size - 1) + (j == 0 || j == grid_size - 1);

			std::vector<Key> nb;
			if (i != 0) nb.push_back(Key(i-0.5, j));
			if (j != 0) nb.push_back(Key(i, j-0.5));
			if (j != grid_size-1) nb.push_back(Key(i, j+0.5));
			if (i != grid_size-1) nb.push_back(Key(i+0.5, j));
			fine_neighbours_list[k] = nb;

			nb.clear();
			if (i != 0 && i != 0.5) nb.push_back(Key(i-1, j));
			if (i != grid_size-1 && i != grid_size-1.5) nb.push_back(Key(i+1, j));
			if (j != 0 && j != 0.5) nb.push_back(Key(i, j-1));
			if (j != grid_size-1 && j != grid_size-1.5) nb.push_back(Key(i, j+1));
			outer_neighbours_list[k] = nb;
		}
	}

	for (auto const& k: fine_grid) {
		for (auto const& n: fine_neighbours_list[k]) {
			g.bond_parameters[make_valid_bond(k, n)] = 0;
		}
	}

	// we extend the cluster map to the finer grid
	for (auto const& k: fine_grid) cluster_map[k] = false;
	for (auto const& k: c) {
		cluster_map[k] = true;
		for (auto const& n: g.neighbours_list[k]) {
			if (g.bond_parameters.at(make_valid_bond(k, n)) == 1 && cluster_map[n]) {
				Key u((k.i+n.i)/2, (k.j+n.j)/2);
				g.bond_parameters[make_valid_bond(k, u)] = 1;
				g.bond_parameters[make_valid_bond(n, u)] = 1;
				cluster_map[u] = true;
			}
		}
	}

	// we create a copy of bond_parameters where all pairs of points in the cluster are closed off
	dense_cluster_map = cluster_map;
	dense_cluster_bond_parameters = g.bond_parameters;
	for (auto const& k: fine_grid) for (auto const& n: fine_neighbours_list[k]) {
		dense_cluster_bond_parameters[make_valid_bond(k, n)] = false;
	}
	for (auto const& [k, in]: cluster_map) if (in) {
		for (auto const& n: outer_neighbours_list[k]) if (cluster_map[n])	{
			dense_cluster_bond_parameters[make_valid_bond(k, n)] = true;
			dense_cluster_bond_parameters[make_valid_half_bond(k, n)] = true;
			dense_cluster_bond_parameters[make_valid_half_bond(n, k)] = true;
			Key u((k.i+n.i)/2, (k.j+n.j)/2);
			dense_cluster_map[u] = true;
		}
 }

}
/*
Bond Hull::make_valid_bond(Key k1, Key k2) {
	return k1.j < k2.j || k1.i < k2.i ? Bond(k1, k2) : Bond(k2, k1);
}
*/

Bond Hull::make_valid_bond(Key k1, Key k2) {
	if (k1.j < k2.j) return Bond(k1, k2);
	else if (k1.j == k2.j && k1.i < k2.i) return Bond(k1, k2);
	return Bond(k2, k1);
}

Bond Hull::make_valid_half_bond(Key k1, Key k2) {
	Key n((k1.i+k2.i)/2, (k1.j+k2.j)/2);
	if (k1.i < k2.i || k1.j < k2.j) return Bond(k1, n);
	return Bond(n, k1);
}

std::vector<Key> Hull::square(Key k, float r) {
	if (r == 0) return fine_neighbours_list[k];
	std::vector<Key> points;
	for (float i = k.i - r; i <= k.i + r; i+=0.25) {
		for (float j = k.j - r; j <= k.j + r; j+=0.25) {
			if (i > 0 && i < grid_size && j > 0 && j < grid_size) points.push_back(Key(i, j));
		}
	}
	return points;
}

float Hull::bond_length(Bond b) {
	Key k1(b.a), k2(b.b);
	return sqrt((k1.i-k2.i)*(k1.i-k2.i)+(k1.j-k2.j)*(k1.j-k2.j));
}

float Hull::eucl_dist(Key k1, Key k2) {
	return sqrt((k1.i-k2.i)*(k1.i-k2.i)+(k1.j-k2.j)*(k1.j-k2.j));
}

float Hull::angle(Key k, Key n, bool counterclockwise) {
	float comp_vec_x = k.i - n.i, comp_vec_y = k.j - n.j;;
	float comp_angle;

	if (counterclockwise) {
		if (n.j < k.j) comp_angle = acos(comp_vec_x/(sqrt(comp_vec_x*comp_vec_x+comp_vec_y*comp_vec_y)));
		else comp_angle = 2*M_PI - acos(comp_vec_x/(sqrt(comp_vec_x*comp_vec_x+comp_vec_y*comp_vec_y)));
	}
	else {
		if (n.j > k.j) comp_angle = acos(comp_vec_x/(sqrt(comp_vec_x*comp_vec_x+comp_vec_y*comp_vec_y)));
		else comp_angle = 2*M_PI - acos(comp_vec_x/(sqrt(comp_vec_x*comp_vec_x+comp_vec_y*comp_vec_y)));
	}
	return comp_angle;
}

bool Hull::parallel_aligned(Bond b1, Bond b2) { // only works for totally vertical or horizontal pairs
	b1 = make_valid_bond(b1.a, b1.b), b2 = make_valid_bond(b2.a, b2.b);
	if (b1.a.i == b1.b.i && b2.a.i == b2.b.i) { // means it's horizontal
		if (b1.a.j == b2.a.j && b1.b.j == b2.b.j) return true;
	}
	else if (b1.a.j == b1.b.j && b2.a.j == b2.b.j) { // means it's vertical
		if (b1.a.i == b2.a.i && b1.b.i == b2.b.i) return true;
	}
	return false;
}

bool Hull::perpendicular(Bond b1, Bond b2) {
	float x1 = b1.a.i - b1.b.i, y1 = b1.a.j-b1.b.j;
	float x2 = b2.a.i - b2.b.i, y2 = b2.a.j-b2.b.j;

	if (x1*x2 + y1*y2 == 0) return true;
	return false;
}

bool Hull::on_segment(Bond b, Key k) {
	if (eucl_dist(b.a, k) <= bond_length(b) && eucl_dist(b.b, k) <= bond_length(b)) return 1;
	return 0;
}

int Hull::crossing_link(Bond b, Key center, float r) {
	// to take care of the case where the circle crosses a bond twice and the two endpoint keys are outside
	// we subdivide the bond with 3 points in between, and we test whether some of those are inside the circle
	Key k1 = b.a, k2 = b.b;
	float bl = bond_length(b),   d1 = eucl_dist(k1, center),   d2 = eucl_dist(k2, center);

	if (std::max(d1, d2) > r+bl) return 0; // preemptively exclude all bonds that CAN'T intersect the circle
	if (std::max(d1, d2) < r) return 0; // bonds with keys both inside the circle can't intersect

	int crossings = 0;
	std::vector<Key> subdivisions;
	for (float t = 0; t <= 1; t += 0.02) subdivisions.push_back(Key(t*k1.i+(1-t)*k2.i, t*k1.j+(1-t)*k2.j));

	for (int i = 1; i < subdivisions.size(); i++) {
		float dist1 = eucl_dist(subdivisions[i-1], center), dist2 = eucl_dist(subdivisions[i], center);
		if ((dist1 >= r && dist2 <= r) || (dist1 <= r && dist2 >= r)) crossings++;
	}

	return crossings;
}

Key Hull::crossing_point(Bond b, Key center, float r) {
	// divides b into 100 subsegments, and approximates the crossing point that way
	// in the unlikely case that the circle cuts the segment in 2 spots, we only return one of those

	Key k1 = b.a, k2 = b.b;

	std::vector<Key> subdivisions;
	for (float t = 0; t <= 1; t += 0.001) subdivisions.push_back(Key(t*k1.i+(1-t)*k2.i, t*k1.j+(1-t)*k2.j));
	subdivisions.push_back(k1);

	for (int i = 1; i < subdivisions.size(); i++) {
		float dist1 = eucl_dist(subdivisions[i-1], center), dist2 = eucl_dist(subdivisions[i], center);
		if ((dist1 >= r && dist2 <= r) || (dist1 <= r && dist2 >= r)) return subdivisions.at(i);
	}

}

int Hull::intercrossing_link(Bond b, Key center, float r) {
	// to take care of the case where the circle crosses a bond twice and the two endpoint keys are outside
	// we subdivide the bond with 3 points in between, and we test whether some of those are inside the circle
	Key k1 = b.a, k2 = b.b;
	float bl = bond_length(b),   d1 = eucl_dist(k1, center),   d2 = eucl_dist(k2, center);

	if (std::max(d1, d2) > r+bl+1) return 0; // preemptively exclude all bonds that CAN'T intersect the circle

	float dx = k2.j - k1.j,   dy = k1.i - k2.i,    dr2 = pow(dx,2)+pow(dy,2);
	float D = k2.j*k1.i-k1.j*k2.i,    Delta = pow(r, 2)*dr2 - pow(D, 2);

	if (std::abs(Delta) < 1.0e-5) { // if discriminant is (almost) zero, tangent -> one intersection point
		float x = D*dy/dr2,   y = -D*dx/dr2;
		if (on_segment(b, Key(center.i - y, center.j + x))) return 1;
	}

	int crossings = 0;
	std::vector<Key> subdivisions;
	for (float t = 0; t <= 1; t += 0.02) subdivisions.push_back(Key(t*k1.i+(1-t)*k2.i, t*k1.j+(1-t)*k2.j));

	for (int i = 1; i < subdivisions.size(); i++) {
		float dist1 = eucl_dist(subdivisions[i-1], center), dist2 = eucl_dist(subdivisions[i], center);
		if ((dist1 >= r && dist2 <= r) || (dist1 <= r && dist2 >= r)) crossings++;
	}

	return crossings;
}

Key Hull::intercrossing_point(Bond b, Key center, float r) {
	// divides b into 100 subsegments, and approximates the crossing point that way
	// in the unlikely case that the circle cuts the segment in 2 spots, we only return one of those

	Key k1 = b.a, k2 = b.b;

	float dx = k2.j - k1.j,   dy = k1.i - k2.i,    dr2 = pow(dx,2)+pow(dy,2);
	float D = k2.j*k1.i-k1.j*k2.i,    Delta = pow(r, 2)*dr2 - pow(D, 2);

	if (std::abs(Delta) < 1.0e-5) { // if discriminant is (almost) zero, tangent -> one intersection point
		float x = D*dy/dr2,   y = -D*dx/dr2;
		if (on_segment(b, Key(center.i - y, center.j + x))) return Key(center.i - y, center.j + x);
	}

	std::vector<Key> subdivisions;
	for (float t = 0; t <= 1; t += 0.001) subdivisions.push_back(Key(t*k1.i+(1-t)*k2.i, t*k1.j+(1-t)*k2.j));
	subdivisions.push_back(k1);

	for (int i = 1; i < subdivisions.size(); i++) {
		float dist1 = eucl_dist(subdivisions[i-1], center), dist2 = eucl_dist(subdivisions[i], center);
		if ((dist1 >= r && dist2 <= r) || (dist1 <= r && dist2 >= r)) return subdivisions.at(i);
	}

}
/*
std::vector<Key> Hull::crossing_test(Bond b, Key center, float r) {
	// keep in mind: i corresponds to y, and j corresponds to -x
	Key k1 = b.a, k2 = b.b;
	float bl = bond_length(b),   d1 = eucl_dist(k1, center),   d2 = eucl_dist(k2, center);

	if (std::max(d1, d2) > r+bl) return {}; // preemptively exclude all bonds that CAN'T intersect the circle
	if (std::max(d1, d2) < r) return {}; // bonds both inside the circle can't intersect

	float dx = k2.j - k1.j,   dy = k1.i - k2.i,    dr2 = pow(dx,2)+pow(dy,2);
	float D = k2.j*k1.i-k1.j*k2.i,    Delta = pow(r, 2)*dr2 - pow(D, 2);

	if (Delta < 0) return {}; // if discriminant is negative, doesn't intersect
	else if (Delta == 0) {// if discriminant is zero, tangent -> one intersection point
		float x = D*dy/dr2,   y = -D*dx/dr2;
		if (on_segment(b, Key(-y,x))) return {Key(center.i+(-1)*y, center.j+x)};
		else return {};
	}

	else if (Delta > 0 && std::min(d1, d2) < r) // if one bond key is inside the circle, it intersects it once
	{
		std::cout << "min<r " << b << std::endl;
		float x1 = (D*dy+sign(dy)*dx*sqrt(Delta))/dr2,   y1 = (-D*dx+std::abs(dy)*sqrt(Delta))/dr2;
		float x2 = (D*dy-sign(dy)*dx*sqrt(Delta))/dr2,   y2 = (-D*dx-std::abs(dy)*sqrt(Delta))/dr2;
		Key p1(center.i+(-1)*y1, center.j + x1),   p2(center.i+(-1)*y2, center.j+x2);

		if (on_segment(b, p1)) return {p1};
		else return {p2};
	}
	else if (Delta > 0 && std::min(d1, d2) >= r) // if both bond keys are outside (or on) the circle, it intersects it twice
	{
		float x1 = (D*dy+sign(dy)*dx*sqrt(Delta))/dr2,   y1 = (-D*dx+std::abs(dy)*sqrt(Delta))/dr2;
		float x2 = (D*dy-sign(dy)*dx*sqrt(Delta))/dr2,   y2 = (-D*dx-std::abs(dy)*sqrt(Delta))/dr2;
		Key p1(center.i+(-1)*y1, center.j + x1),   p2(center.i+(-1)*y2, center.j+x2);
		std::cout << "p1: " << p1 << ",    p2: " << p2 << std::endl;
		if (on_segment(b, p1) && on_segment(b, p2)) {
			// return the one that's closer in terms of hull-distance
			Key first = ordered_hull_keys[b.a] < ordered_hull_keys[b.b] ? b.a : b.b;
			Key p = eucl_dist(p1, first) < eucl_dist(p2, first) ? p1 : p2;
			return {p};
		}
		else return {};
	}
}

int Hull::sign(float n) {
	if (n < 0) return -1;
	return 1;
}
*/
/*
Key Hull::crossing_point(Bond bond, Key center, float r) {
	float x_a = bond.a.j, x_b = bond.b.j, x_c = center.j;
	float y_a = bond.a.i, y_b = bond.b.i, y_c = center.i;

	float a = -x_a*x_a + 2*x_a*x_b - x_b*x_b - y_a*y_a + 2*y_a*y_b - y_b*y_b;
	float b = -2*x_a*x_b + 2*x_a*x_c + 2*x_b*x_b - 2*x_b*x_c - 2*y_a*y_b + 2*y_a*y_c + 2*y_b*y_b - 2*y_b*y_c;
	float c = r*r - x_b*x_b + 2*x_b*x_c - x_c*x_c - y_b*y_b + 2*y_b*y_c - y_c*y_c;

	float t1 = (-b + sqrt(b*b - 4*a*c))/(2*a);
	float t2 = (-b - sqrt(b*b - 4*a*c))/(2*a);

	if (t1 >= 0 && t1 <= 1) return Key(t1*x_a + (1-t1)*x_b, t1*y_a + (1-t1)*y_b);
	else if (t2 >= 0 && t2 <= 1) return Key(t2*x_a + (1-t2)*x_b, t2*y_a + (1-t2)*y_b);

}
*/
bool Hull::common_element(std::vector<Key> v1, std::vector<Key> v2) {
	for (auto const& k1: v1) for (auto const& k2: v2) {
		if (k1 == k2) return true;
	}
	return false;
}

bool Hull::inbetween_key(Key k1, Key k2, std::map<Key, bool> list) {
	if (k1.i == k2.i) { //segment is horizontal
		for (auto const& [k, exists]: list) if (exists) {
			if (k.i == k1.i) {
				if ((k.j > k1.j && k.j < k2.j) || (k.j < k1.j && k.j > k2.j)) return true;
			}
		}
	}
	if (k1.j == k2.j) {
		for (auto const& [k, exists]: list) if (exists) {
			if (k.j == k1.j) {
				if ((k.i > k1.i && k.i < k2.i) || (k.i < k1.i && k.i > k2.i)) return true;
			}
		}
	}
	return false;
}

bool Hull::permitted_move(Key start, Key end) {
	if (cluster_map[end] && cluster_map[start]) return 0;
	else if (cluster_map[end]) return 0;
	else return 1;
}

bool Hull::dense_permitted_move(Key start, Key end) {
	if (dense_cluster_map[end] && dense_cluster_map[start]) return 0;
	else if (dense_cluster_map[end]) return 0;
	else return 1;
}

std::vector<Bond> Hull::determine_hull() {
	std::cout << "CALCULATING HULL..." << std::endl << std::flush;

	if (!hull.empty()) hull.clear();
	std::vector<Bond> h;

	Key start(0, grid_size/2);
	if (cluster_map[start]) start = Key(0, grid_size/2-1);

	std::map<Key, bool> visited;
	std::stack<Key> stack;
	stack.push(start);

	while (!stack.empty()) {
		// Pop a vertex from stack
		Key s = stack.top();
		stack.pop();

		if (visited[s] != true) visited[s] = true;

		std::vector<Key> nb = fine_neighbours_list[s];
		for (auto const& v: nb) {
			if (visited[v] != true && permitted_move(s, v)) stack.push(v);
			if (cluster_map[v]) {
				Key u((v.i+s.i)/2, (v.j+s.j)/2);
				h.push_back(make_valid_bond(u, v));
			}
		}
	}

	std::cout << "DONE" << std::endl << std::flush;

	hull = h;
	return h;
}

std::vector<Bond> Hull::determine_ap() {
	std::cout << "CALCULATING AP..." << std::endl << std::flush;
	Key start(0, floor(grid_size/2));
	if (cluster_map[start]) start = Key(0, floor(grid_size/2)-1);

	std::map<Key, bool> visited;
	std::stack<Key> stack;
	stack.push(start);

	while (!stack.empty()) {
		// Pop a vertex from stack
		Key s = stack.top();
		stack.pop();

		if (visited[s] != true) visited[s] = true;

		std::vector<Key> nb = fine_neighbours_list[s];
		for (auto const& v: nb) {
			if (visited[v] != true && dense_permitted_move(s, v)) stack.push(v);
			if (dense_cluster_map[v]) {
				Key u((v.i+s.i)/2, (v.j+s.j)/2);
				ap.push_back(make_valid_bond(u, v));
			}
		}
	}
		std::cout << "DONE" << std::endl << std::flush;
	return ap;
}

std::vector<Bond> Hull::old_link_ap() {

	std::cout << "STARTING AP LINKAGE..." << std::endl << std::flush;

	time_t tstart, tend;
	tstart = time(0);

	if (ap.empty()) determine_ap();
	std::vector<Bond> ap_graph;
	std::map<Bond, bool> ap_map;
	std::map<Bond, bool> visited;
	std::map<Key, std::vector<Key>> support; // takes an AP key and returns its corresponding supporting key(s)

	// find the non-cluster points
	for (auto const& b: ap) {
		if (!dense_cluster_map[b.a] == true) {
			ap_points[b.a] = true;
			support[b.a].push_back(b.b);
		}
		else if (!dense_cluster_map[b.b] == true) {
			ap_points[b.b] = true;
			support[b.b].push_back(b.a);
		}
		ap_map[b] = true;
	}

	for (auto const& [k, v_s]: support) {
		std::vector<Key> neighbourhood;
		std::vector<Bond> candidates;

		for (auto const& n: square(k, 0.75)) if (ap_points[n] && k != n) neighbourhood.push_back(n);

		for (auto const& n: neighbourhood) {

			std::vector<Key> v_s_n = support[n];
			for (auto const& s: v_s) for (auto const& s_n: v_s_n) {
				if (dense_cluster_bond_parameters[make_valid_bond(s_n, s)]) {
					if (!visited[make_valid_bond(k, n)] && parallel_aligned(make_valid_bond(s_n, s), make_valid_bond(k, n))) {
						if (inbetween_key(k, n, ap_points) == false) {
							if (!common_element(v_s, v_s_n)) {
								ap_graph.push_back(make_valid_bond(k, n));
								visited[make_valid_bond(k, n)] = true;
							}
						}
					}
				}
				else if (!visited[make_valid_bond(k, n)] && perpendicular(make_valid_bond(s_n, n), make_valid_bond(s, k))) {
					if (s_n == s && k != n) {
						ap_graph.push_back(make_valid_bond(k, n));
						visited[make_valid_bond(k, n)] = true;
					}
					else {
						std::vector<Key> clust_k, clust_n;
						for (auto const& u: fine_neighbours_list[s]) if (dense_cluster_map[u]) clust_k.push_back(u);
						for (auto const& u: fine_neighbours_list[s_n]) if (dense_cluster_map[u]) clust_n.push_back(u);

						float dist = eucl_dist(k, n);
						if (common_element(clust_k, clust_n) && dist < 0.36) { // only diagonal links we want are the 0.25*sqrt(2)-length ones
							ap_graph.push_back(make_valid_bond(k, n));
							visited[make_valid_bond(k, n)] = true;
						}
					}
				}
			}
		}
	}

	Key first(0, 0), last(0, 0), current(0, 0);
	int index = 1;
	for (auto const& [k, exists]: ap_points) if (exists) {
		if (k.i == 0) {
			first = k, current = first;
			ordered_ap_keys[k] = index;
			index++;
		}
		else if (k.i == grid_size-1) last = k;
	}

	while (current != last) {
		for (auto& [b, exists]: visited) if (exists) {
			if (b.a == current) {
				ordered_ap_keys[b.b] = index;
				current = b.b;
				exists = false;
				index++;
				break;
			}
			else if (b.b == current) {
				ordered_ap_keys[b.a] = index;
				current = b.a;
				exists = false;
				index++;
				break;
			}
		}
	}

	std::cout << "DONE" << std::endl << std::flush;

	tend = time(0);
	std::cout << "It took " << difftime(tend, tstart) << " to compute and link the external hull." << std::endl;

	linked_ap = ap_graph;
	return ap_graph;
}

std::vector<Bond> Hull::link_ap() {

	std::cout << "STARTING AP LINKAGE..." << std::endl << std::flush;

	time_t tstart, tend;
	tstart = time(0);

	if (ap.empty()) determine_ap();
	std::vector<Bond> ap_graph;
	std::map<Bond, bool> ap_map;
	std::map<Bond, bool> visited;
	std::map<Key, std::vector<Key>> support; // takes an AP key and returns its corresponding supporting key(s)

	// find the non-cluster points
	for (auto const& b: ap) {
		if (!dense_cluster_map[b.a] == true) {
			ap_points[b.a] = true;
			support[b.a].push_back(b.b);
		}
		else if (!dense_cluster_map[b.b] == true) {
			ap_points[b.b] = true;
			support[b.b].push_back(b.a);
		}
		ap_map[b] = true;
	}

	for (auto const& [k, v_s]: support) {
		std::vector<Key> neighbourhood;
		std::vector<Bond> candidates;

		for (auto const& n: square(k, 0.75)) if (ap_points[n] && k != n) {
			if (eucl_dist(k, n) == 0.5) neighbourhood.push_back(n);
			else if (eucl_dist(k, n) < 0.36 && eucl_dist(k, n) > 0.35) neighbourhood.push_back(n);
		}

		for (auto const& n: neighbourhood) {

			std::vector<Key> v_s_n = support[n];
			for (auto const& s: v_s) for (auto const& s_n: v_s_n) {
				if (eucl_dist(k, n) == 0.5) {
					if (!visited[make_valid_bond(k, n)] && parallel_aligned(make_valid_bond(s_n, s), make_valid_bond(k, n))) {
							ap_graph.push_back(make_valid_bond(k, n));
							visited[make_valid_bond(k, n)] = true;
						}
				}

				else if (eucl_dist(k, n) < 0.36 && eucl_dist(k, n) > 0.35) {
					std::vector<Key> clust_k, clust_n;
					for (auto const& u: fine_neighbours_list[s]) if (dense_cluster_map[u]) clust_k.push_back(u);
					for (auto const& u: fine_neighbours_list[s_n]) if (dense_cluster_map[u]) clust_n.push_back(u);

					float dist = eucl_dist(k, n);
					if (!visited[make_valid_bond(k, n)] && common_element(clust_k, clust_n) && dist < 0.36) { // only diagonal links we want are the 0.25*sqrt(2)-length ones
						ap_graph.push_back(make_valid_bond(k, n));
						visited[make_valid_bond(k, n)] = true;
					}
				}
			}
		}
	}

	Key first(0, 0), last(0, 0), current(0, 0);
	float total_distance = 0;
	for (auto const& [k, exists]: ap_points) if (exists) {
		if (k.i == 0) {
			first = k, current = first;
			ordered_hull_keys[k] = total_distance;
		}
		else if (k.i == grid_size-1) last = k;
	}

	while (current != last) {
		for (auto& [b, exists]: visited) if (exists) {
			if (b.a == current) {
				total_distance += eucl_dist(b.b, current);
				ordered_ap_keys[b.a] = total_distance;
				ordered_ap.push_back(b.b);
				current = b.b;
				exists = false;
				break;
			}
			else if (b.b == current) {
				total_distance += eucl_dist(b.a, current);
				ordered_ap_keys[b.a] = total_distance;
				ordered_ap.push_back(b.a);
				current = b.a;
				exists = false;
				break;
			}
		}
	}

	std::cout << "DONE" << std::endl << std::flush;

	tend = time(0);
	std::cout << "It took " << difftime(tend, tstart) << " to compute and link the external hull." << std::endl;

	linked_ap = ap_graph;
	if (g.grid_size > 50) write_ap();
	return ap_graph;
}

std::vector<Bond> Hull::old_link_hull() {

	std::cout << "STARTING HULL LINKAGE..." << std::endl << std::flush;

	time_t tstart, tend;
	tstart = time(0);

	if (hull.empty()) determine_hull();
	std::vector<Bond> hull_graph;
	std::map<Bond, bool> hull_map;
	std::map<Bond, bool> visited;
	std::map<Key, std::vector<Key>> support; // takes an AP key and returns its corresponding supporting key(s)

	// find the non-cluster points
	for (auto const& b: hull) {
		if (!cluster_map[b.a] == true) {
			hull_points[b.a] = true;
			support[b.a].push_back(b.b);
		}
		else if (!cluster_map[b.b] == true) {
			hull_points[b.b] = true;
			support[b.b].push_back(b.a);
		}
		hull_map[b] = true;
	}

	// for each point of the hull, find its natural neighbours using the white scaffolding
	for (auto const& [k, v_s]: support) {
		std::vector<Key> neighbourhood;
		std::vector<Bond> candidates;

		for (auto const& n: square(k, 0.75)) if (hull_points[n] && k != n) {
			neighbourhood.push_back(n);
			//if (eucl_dist(k, n) == 0.5) neighbourhood.push_back(n);
			//else if (eucl_dist(k, n) < 0.36 && eucl_dist(k, n) > 0.35) neighbourhood.push_back(n);
		}

		for (auto const& n: neighbourhood) {

			std::vector<Key> v_s_n = support[n];
			for (auto const& s: v_s) for (auto const& s_n: v_s_n) {
				if (g.bond_parameters[make_valid_bond(s_n, s)]) {
					if (!visited[make_valid_bond(k, n)] && parallel_aligned(make_valid_bond(s_n, s), make_valid_bond(k, n))) {
						if (inbetween_key(k, n, hull_points) == false) {
							if (!common_element(v_s, v_s_n)) {
								hull_graph.push_back(make_valid_bond(k, n));
								visited[make_valid_bond(k, n)] = true;
							}
						}
					}
				}
				else if (!visited[make_valid_bond(k, n)] && perpendicular(make_valid_bond(s_n, n), make_valid_bond(s, k))) {
					if (s_n == s && k != n) {
						hull_graph.push_back(make_valid_bond(k, n));
						visited[make_valid_bond(k, n)] = true;
					}
					else {
						std::vector<Key> clust_k, clust_n;
						for (auto const& u: fine_neighbours_list[s]) if (cluster_map[u]) clust_k.push_back(u);
						for (auto const& u: fine_neighbours_list[s_n]) if (cluster_map[u]) clust_n.push_back(u);

						float dist = eucl_dist(k, n);
						if (common_element(clust_k, clust_n) && dist < 0.36) { // only diagonal links we want are the 0.25*sqrt(2)-length ones
							hull_graph.push_back(make_valid_bond(k, n));
							visited[make_valid_bond(k, n)] = true;
						}
					}
				}
			}
		}
	}

	// order out all of the hull keys using the previously generated graph
	Key first(0, 0), last(0, 0), current(0, 0);
	float total_distance = 0;
	for (auto const& [k, exists]: hull_points) if (exists) {
		if (k.i == 0) {
			first = k, current = first;
			ordered_hull_keys[k] = total_distance;
		}
		else if (k.i == grid_size-1) last = k;
	}

	while (current != last) {
		for (auto& [b, exists]: visited) if (exists) {
			if (b.a == current) {
				total_distance += eucl_dist(b.b, current);
				ordered_hull_keys[b.b] = total_distance;
				current = b.b;
				exists = false;
				break;
			}
			else if (b.b == current) {
				total_distance += eucl_dist(b.a, current);
				ordered_hull_keys[b.a] = total_distance;
				current = b.a;
				exists = false;
				break;
			}
		}
	}

	std::cout << "DONE" << std::endl << std::flush;
	tend = time(0);
	std::cout << "It took " << difftime(tend, tstart) << " to compute and link the complete hull." << std::endl;

	linked_hull = hull_graph;
	return linked_hull;

}

std::vector<Bond> Hull::link_hull() {

	std::cout << "STARTING HULL LINKAGE..." << std::endl << std::flush;

	time_t tstart, tend;
	tstart = time(0);

	if (hull.empty()) determine_hull();
	std::vector<Bond> hull_graph;
	std::map<Bond, bool> hull_map;
	std::map<Bond, bool> visited;
	std::map<Key, std::vector<Key>> support; // takes an AP key and returns its corresponding supporting key(s)

	// find the non-cluster points
	for (auto const& b: hull) {
		if (!cluster_map[b.a] == true) {
			hull_points[b.a] = true;
			support[b.a].push_back(b.b);
		}
		else if (!cluster_map[b.b] == true) {
			hull_points[b.b] = true;
			support[b.b].push_back(b.a);
		}
		hull_map[b] = true;
	}

	// for each point of the hull, find its natural neighbours using the white scaffolding
	for (auto const& [k, v_s]: support) {
		std::vector<Key> neighbourhood;
		std::vector<Bond> candidates;

		for (auto const& n: square(k, 0.75)) if (hull_points[n] && k != n) {
			// make a list of the eligible neighbours (only 0.5 vert/hor and 0.353 diagonal)
			if (eucl_dist(k, n) == 0.5) neighbourhood.push_back(n);
			else if (eucl_dist(k, n) < 0.36 && eucl_dist(k, n) > 0.35) neighbourhood.push_back(n);
		}

		for (auto const& n: neighbourhood) {

			std::vector<Key> v_s_n = support[n];
			for (auto const& s: v_s) for (auto const& s_n: v_s_n) {
				if (eucl_dist(k, n) == 0.5) {
					// for 0.5 vert/hor links, they have to be non-visited and perpendicular to their scaffolds
					if (!visited[make_valid_bond(k, n)] && parallel_aligned(make_valid_bond(s_n, s), make_valid_bond(k, n))) {
							hull_graph.push_back(make_valid_bond(k, n));
							visited[make_valid_bond(k, n)] = true;
						}
				}
					// otherwise for diagonal links, we distinguish 2 cases
				else if (eucl_dist(k, n) < 0.36 && eucl_dist(k, n) > 0.35) {
					std::vector<Key> clust_k, clust_n;
					for (auto const& u: fine_neighbours_list[s]) if (cluster_map[u]) clust_k.push_back(u);
					for (auto const& u: fine_neighbours_list[s_n]) if (cluster_map[u]) clust_n.push_back(u);

					if (!visited[make_valid_bond(k, n)] && common_element(clust_k, clust_n)) {
						hull_graph.push_back(make_valid_bond(k, n));
						visited[make_valid_bond(k, n)] = true;
					}
				}
			}
		}
	}


	// order out all of the hull keys using the previously generated graph
	Key first(0, 0), last(0, 0), current(0, 0);
	float total_distance = 0;
	for (auto const& [k, exists]: hull_points) if (exists) {
		if (k.i == 0) {
			first = k, current = first;
			ordered_hull_keys[k] = total_distance;
			ordered_hull.push_back(k);
		}
		else if (k.i == grid_size-1) last = k;
	}

	while (current != last) {
		for (auto& [b, exists]: visited) if (exists) {
			if (b.a == current) {
				total_distance += eucl_dist(b.b, current);
				ordered_hull_keys[b.b] = total_distance;
				ordered_hull.push_back(b.b);
				current = b.b;
				exists = false;
				break;
			}
			else if (b.b == current) {
				total_distance += eucl_dist(b.a, current);
				ordered_hull_keys[b.a] = total_distance;
				ordered_hull.push_back(b.a);
				current = b.a;
				exists = false;
				break;
			}
		}
	}

	std::cout << "DONE" << std::endl << std::flush;
	tend = time(0);
	std::cout << "It took " << difftime(tend, tstart) << " to compute and link the complete hull." << std::endl;

	linked_hull = hull_graph;
	if (g.grid_size > 50) write_hull();
	return linked_hull;

}

void Hull::write_hull() {
	if (linked_hull.empty()) link_hull();
	//namespace fs = std::filesystem;
	//fs::create_directory("data/hull/q="+ std::to_string(q) +"/"+ std::to_string(gs));

	std::string name = gen_name();
	std::ofstream hfile;
	int q = static_cast<int>(g.q), gs = static_cast<int>(g.grid_size);
	hfile.open("raw_grids/hull/"+name+"_gs=" + std::to_string(gs) + ".csv", std::ios_base::app);

	hfile << "#Ns= " << ordered_hull.size()-1 << std::endl;
	hfile << "# Hull of a fixed-boundary cluster" << std::endl;
	hfile << "#gs=" << grid_size << std::endl << "#q=" << q << std::endl << "#ne=" << g.evolve_steps;

	for (int i = 0; i < ordered_hull.size()-1; i++) {
		Key k1 = ordered_hull[i], k2 = ordered_hull[i+1];
		hfile << k1.j << " " << k1.i << " " << k2.j << " " << k2.i << std::endl;
	}

	hfile.close();
}

void Hull::write_ap() {
	std::string name = gen_name();
	std::ofstream apfile;
	int q = static_cast<int>(g.q), gs = static_cast<int>(g.grid_size);
	apfile.open("raw_grids/ap/"+name+"_gs=" + std::to_string(gs) + ".csv", std::ios_base::app);

	apfile << "#Ns= " << ordered_ap.size()-1 << std::endl;
	apfile << "# Accessible perimeter of a fixed-boundary cluster" << std::endl;
	apfile << "#gs=" << grid_size << std::endl << "#q=" << q << std::endl << "#ne=" << g.evolve_steps;

	for (int i = 0; i < ordered_ap.size()-1; i++)	apfile << ordered_ap[i].j << std::endl;

	apfile.close();
}

float Hull::hull_length() {
	if (linked_hull.empty()) link_hull();
	float dist = 0;
	for (auto const& b: linked_hull) {
		dist += bond_length(b);
	}
	return dist;
}

float Hull::hull_width() {
	if (linked_hull.empty()) link_hull();

	float width = 0, mean_col = 0;;

	for (auto const& k: ordered_hull) mean_col += k.j;
	mean_col = mean_col/hull.size();

	for (auto const& k: ordered_hull) {
		width += pow(mean_col - k.j, 2);
	}

	return sqrt(width/hull.size());
}

float Hull::ap_length() {
	if (linked_ap.empty()) link_ap();
	float dist = 0;
	for (auto const& b: linked_ap) dist += bond_length(b);
	return dist;
}

float Hull::ap_width() {
	if (ap.empty()) determine_ap();
	float width = 0, mean_col = 0;;

	for (auto const& b: ap) {
		Key k = cluster_map[b.a] ? b.b : b.a;
		mean_col += k.j;
	}
	mean_col = mean_col/ap.size();

	for (auto const& b: ap) {
		Key k = cluster_map[b.a] ? b.b : b.a;
		width += std::abs(mean_col - k.j);
	}

	return width/ap.size();
}

std::vector<std::map<double, double>> Hull::yardstick_ap(std::vector<double> radii) {
	std::map<double, double> hull_data, ap_data;

	time_t tstart, tend;
	tstart = time(0);

	if (hull.empty()) determine_hull();
	if (linked_hull.empty()) link_hull();
	if (ap.empty()) determine_ap();
	if (linked_ap.empty()) link_ap();

	Key first(0, 0), last(0, 0), current(0, 0);
	for (auto const& [k, exists]: hull_points) if (exists) {
		if (k.i == 0) first = k;
		else if (k.i == grid_size-1) last = k;
	}

	for (auto const& r: radii) {
		std::map<Bond, double> candidate_crossings; // maps the crossed bond and the smallest graph length of its two keys
		double total_distance = 0;
		double current_length = 0;
		current = first;

		do {
			candidate_crossings.clear();

			for (auto const& b: linked_hull) {
				bool crossed = intercrossing_link(b, current, r);
				if (crossed) {
					Key closest = ordered_hull_keys[b.a] < ordered_hull_keys[b.b] ? b.a : b.b;
					double distance = ordered_hull_keys[closest];
					if (distance > current_length) candidate_crossings[b] = distance;
				}
			}

			// treating the case of no crossings because we've reached the end
			if (candidate_crossings.empty() && eucl_dist(current, last) <= r) total_distance += eucl_dist(current, last);
			else if (candidate_crossings.empty() && eucl_dist(current, last) > r) {
				std::cout << "No candidate crossings, but distance to end is larger than " << r << " (from " << current << ") " << std::endl;
				display_yardstick_h("nocrossyardstick.png", r);
			}
			else {
				Bond crossed_bond = candidate_crossings.begin()->first;
				double min_length = candidate_crossings.begin()->second;
				Key temp_current(0, 0);

				//std::cout << "Candidates: " << candidate_crossings.size();

				for (auto const& [b, len]: candidate_crossings) {
					//std::cout << b << " :: " << len << std::endl;
					if (len <= min_length) {
						min_length = len;
						temp_current = intercrossing_point(b, current, r);
						crossed_bond = b;
					}
				}
				total_distance += r, current = temp_current;
				current_length = std::min(ordered_hull_keys[crossed_bond.a], ordered_hull_keys[crossed_bond.b]);
			}
		} while (!candidate_crossings.empty());
		hull_data[r] = total_distance;


		total_distance = 0;
		current_length = 0;
		current = first;

		do {
			candidate_crossings.clear();

			for (auto const& b: linked_ap) {
				bool crossed = intercrossing_link(b, current, r);
				if (crossed) {
					Key closest = ordered_ap_keys[b.a] < ordered_ap_keys[b.b] ? b.a : b.b;
					double distance = ordered_ap_keys[closest];
					if (distance > current_length) candidate_crossings[b] = distance;
				}
			}

			// treating the case of no crossings because we've reached the end
			if (candidate_crossings.empty() && eucl_dist(current, last) <= r) total_distance += eucl_dist(current, last);
			else if (candidate_crossings.empty() && eucl_dist(current, last) > r) {
				std::cout << "No candidate crossings, but distance to end is larger than " << r << " (from " << current << ") " << std::endl;
				display_yardstick_ap("nocrossyardstick.png", r);
			}
			else {
				Bond crossed_bond = candidate_crossings.begin()->first;
				double min_length = candidate_crossings.begin()->second;
				Key temp_current(0, 0);

				//std::cout << "Candidates: " << candidate_crossings.size();

				for (auto const& [b, len]: candidate_crossings) {
					//std::cout << b << " :: " << len << std::endl;
					if (len <= min_length) {
						min_length = len;
						temp_current = intercrossing_point(b, current, r);
						crossed_bond = b;
					}
				}
				total_distance += r, current = temp_current;
				current_length = std::min(ordered_ap_keys[crossed_bond.a], ordered_ap_keys[crossed_bond.b]);
			}
		} while (!candidate_crossings.empty());
		ap_data[r] = total_distance;

	}

	tend = time(0);
	difftime(tend, tstart);
	std::cout << "It took " << difftime(tend, tstart) << "s to compute the yardstick." << std::endl;
	return {hull_data, ap_data};
}

std::vector<std::map<double, double>> Hull::yardstick_new(std::vector<double> radii) {
	std::map<double, double> hull_data, ap_data;

	time_t tstart, tend;
	tstart = time(0);

	if (hull.empty()) determine_hull();
	if (linked_hull.empty()) link_hull();

	Key first(0, 0), last(0, 0), current(0, 0);
	for (auto const& [k, exists]: hull_points) if (exists) {
		if (k.i == 0) first = k;
		else if (k.i == grid_size-1) last = k;
	}

	for (auto const& r: radii) {
		std::map<Bond, double> candidate_crossings; // maps the crossed bond and the smallest graph length of its two keys
		double total_distance = 0;
		double current_length = 0;
		current = first;

		do {
			candidate_crossings.clear();

			for (auto const& b: linked_hull) {
				bool crossed = intercrossing_link(b, current, r);
				if (crossed) {
					Key closest = ordered_hull_keys[b.a] < ordered_hull_keys[b.b] ? b.a : b.b;
					double distance = ordered_hull_keys[closest];
					if (distance > current_length) candidate_crossings[b] = distance;
				}
			}

			// treating the case of no crossings because we've reached the end
			if (candidate_crossings.empty() && eucl_dist(current, last) <= r) total_distance += eucl_dist(current, last);
			else if (candidate_crossings.empty() && eucl_dist(current, last) > r) {
				std::cout << "No candidate crossings, but distance to end is larger than " << r << " (from " << current << ") " << std::endl;
				display_yardstick_h("nocrossyardstick.png", r);
			}
			else {
				Bond crossed_bond = candidate_crossings.begin()->first;
				double min_length = candidate_crossings.begin()->second;
				Key temp_current(0, 0);

				//std::cout << "Candidates: " << candidate_crossings.size();

				for (auto const& [b, len]: candidate_crossings) {
					//std::cout << b << " :: " << len << std::endl;
					if (len <= min_length) {
						min_length = len;
						temp_current = intercrossing_point(b, current, r);
						crossed_bond = b;
					}
				}
				total_distance += r, current = temp_current;
				current_length = std::min(ordered_hull_keys[crossed_bond.a], ordered_hull_keys[crossed_bond.b]);
			}
		} while (!candidate_crossings.empty());
		hull_data[r] = total_distance;


		total_distance = 0;
		current_length = 0;
		current = first;

		do {
			candidate_crossings.clear();

			for (auto const& b: linked_hull) {
				bool crossed = intercrossing_link(b, current, r);
				if (crossed) {
					Key closest = ordered_hull_keys[b.a] < ordered_hull_keys[b.b] ? b.a : b.b;
					double distance = ordered_hull_keys[closest];
					if (distance > current_length) candidate_crossings[b] = distance;
				}
			}

			// treating the case of no crossings because we've reached the end
			if (candidate_crossings.empty() && eucl_dist(current, last) <= r) total_distance += eucl_dist(current, last);
			else if (candidate_crossings.empty() && eucl_dist(current, last) > r) {
				std::cout << "No candidate crossings, but distance to end is larger than " << r << " (from " << current << ") " << std::endl;
				display_yardstick_ap("nocrossyardstick.png", r);
			}
			else {
				Bond crossed_bond = candidate_crossings.begin()->first;
				double min_length = candidate_crossings.begin()->second;
				Key temp_current(0, 0);

				//std::cout << "Candidates: " << candidate_crossings.size();

				for (auto const& [b, len]: candidate_crossings) {
					//std::cout << b << " :: " << len << std::endl;
					if (len >= min_length) {
						min_length = len;
						temp_current = intercrossing_point(b, current, r);
						crossed_bond = b;
					}
				}
				total_distance += r, current = temp_current;
				current_length = std::min(ordered_hull_keys[crossed_bond.a], ordered_hull_keys[crossed_bond.b]);
			}
		} while (!candidate_crossings.empty());
		ap_data[r] = total_distance;

	}

	tend = time(0);
	difftime(tend, tstart);
	std::cout << "It took " << difftime(tend, tstart) << "s to compute the yardstick." << std::endl;
	return {hull_data, ap_data};
}

std::vector<std::map<float, float>> Hull::yardstick(std::vector<float> radii) {
	std::map<float, float> hull_data, ap_data;

	time_t tstart, tend;
	tstart = time(0);

	if (hull.empty()) determine_hull();
	if (linked_hull.empty()) link_hull();

	for(auto const& r: radii) std::cout << r << " ";

	Key first(0, 0), last(0, 0), current(0, 0);
	for (auto const& [k, exists]: hull_points) if (exists) {
		if (k.i == 0) first = k;
		else if (k.i == grid_size-1) last = k;
	}

	for (auto const& r: radii) {
		std::map<Bond, float> candidate_crossings; // maps the crossed bond and the smallest graph length of its two keys
		std::map<Key, bool> visited;
		float total_distance = 0;
		for (auto const& [k, exists]: hull_points) if (exists) visited[k] = false;
		current = first;

		// we first make a loop that treats the hull (i.e. takes the last crossing on circles)
		do {
			candidate_crossings.clear();
			// find all the spots where the circle crosses a bond
			for (auto const& b: linked_hull) {
				bool crossed = crossing_link(b, current, r);
				if ((!visited[b.a] || !visited[b.b]) && crossed) {
					candidate_crossings[b] = std::min(ordered_hull_keys[b.a], ordered_hull_keys[b.b]);
				}
				//std::cout << "candidate number: " << candidate_crossings.size() << std::endl;
				//for (auto const& [b,f]: candidate_crossings) std::cout << b;
			}
			// if we get no crossings, that normally means we've reached the end of the graph
			if (candidate_crossings.empty() && eucl_dist(current, last) < r) total_distance += eucl_dist(current, last);
			else if (candidate_crossings.empty() && eucl_dist(current, last) > r) {
				std::cout << "WTFa" << r << " :: " << current << std::endl;
				display_yardstick_h("WTFyardstick.png", r);
			}
			else {
				Bond crossed_bond = candidate_crossings.begin()->first;
				float min_length = candidate_crossings.begin()->second;
				Key temp_current(0, 0);

			 	for (auto const& [b, length]: candidate_crossings) {
					if (length <= min_length) {
						min_length = length, temp_current = crossing_point(b, current, r), crossed_bond = b;
					}
				}

				if (temp_current.i < 0 || temp_current.j < 0) std::cout << "WTFFa" << r << " :: " << current << std::endl;
				total_distance += r, current = temp_current;
				for (auto const& [k, length]: ordered_hull_keys) {
					if (length < candidate_crossings[crossed_bond]) visited[k] = true;
				}
			}
		} while (!candidate_crossings.empty());
		hull_data[r] = total_distance;

		// then we take care of the accessible perimeter (i.e. first crossing on circles)
		total_distance = 0;
		visited.clear();
		for (auto const& [k, exists]: hull_points) if (exists) visited[k] = false;
		current = first;

		do {
			candidate_crossings.clear();
			// find all the spots where the circle crosses a bond
			for (auto const& b: linked_hull) {
				int crossed = crossing_link(b, current, r);
				if ((!visited[b.a] || !visited[b.b]) && crossed) {
					candidate_crossings[b] = std::min(ordered_hull_keys[b.a], ordered_hull_keys[b.b]);
				}
			}
			// if we get no crossings, that means we've reached the end of the graph
			if (candidate_crossings.empty() && eucl_dist(current, last) < r) total_distance += eucl_dist(current, last);
			else if (candidate_crossings.empty() && eucl_dist(current, last) > r) {
				std::cout << "WTFb" << r << " :: " << current << std::endl;
			}
			else {
				Bond crossed_bond = candidate_crossings.begin()->first;
				float max_length = candidate_crossings.begin()->second;
				Key temp_current(0, 0);

				for (auto const& [b, length]: candidate_crossings) {
					if (length >= max_length) {
						max_length = length, temp_current = crossing_point(b, current, r), crossed_bond = b;
					}
				}
				if (temp_current.i < 0 || temp_current.j < 0) std::cout << "WTFFb" << r << " :: " << current << std::endl;
				total_distance += r, current = temp_current;
				for (auto const& [k, length]: ordered_hull_keys) {
					if (length < candidate_crossings[crossed_bond]) visited[k] = true;
				}
			}
		} while (!candidate_crossings.empty());
		ap_data[r] = total_distance;
	}

	for (int i = 0; i < radii.size()-1; i++) if (hull_data[radii[i]] > 8*hull_data[radii[i+1]]) {
		std::cout << "WTF" << radii[i] << " :: " << current << std::endl;
		display_yardstick_h("WTFyardstick"+ std::to_string(radii[i+1]) + ".png", radii[i+1]);
	}

	std::vector<std::map<float, float>> result;
	result.push_back(hull_data), result.push_back(ap_data);

	tend = time(0);
	difftime(tend, tstart);
	std::cout << "It took " << difftime(tend, tstart) << "s to compute the yardstick." << std::endl;
	return result;
}

std::vector<Key> Hull::yardstick_rilder_test(double stick) {

	double *x1,*y1,*x2,*y2,px,py,r;
  double sx,sy,s1x,s1y,s2x,s2y,t,t1,t2,epsilon;
	double l;
  long Ns;
	long i,j,ok,g;

	Key first(0, 0), last(0, 0), current(0, 0);
	for (auto const& [k, exists]: hull_points) if (exists) {
		if (k.i == 0) first = k;
		else if (k.i == grid_size-1) last = k;
	}

	std::vector<Key> keys;
	keys.push_back(ordered_hull[0]);

	Ns = ordered_hull.size()-1;
	x1=(double*)calloc(Ns,sizeof(double));
  y1=(double*)calloc(Ns,sizeof(double));
  x2=(double*)calloc(Ns,sizeof(double));
  y2=(double*)calloc(Ns,sizeof(double));

	for (int i = 0; i < Ns; i++) {
		x1[i] = ordered_hull[i].j;
		y1[i] = ordered_hull[i].i;
		x2[i] = ordered_hull[i+1].j;
		y2[i] = ordered_hull[i+1].i;
	}
	/* Main Loop: ------------------------------------------------------------*/

/* Check if the segments of the object are in the right orientation */
  for(i=1;i<Ns;i++){
    j=i; ok=1;
    while((j<Ns)&&(ok==1)){
      if((fabs(x1[j]-x2[i-1])<1e-8*stick)&&(fabs(y1[j]-y2[i-1])<1e-8*stick)){
        px=x1[j]; x1[j]=x1[i]; x1[i]=px;
        py=y1[j]; y1[j]=y1[i]; y1[i]=py;
        px=x2[j]; x2[j]=x2[i]; x2[i]=px;
        py=y2[j]; y2[j]=y2[i]; y2[i]=py;
        ok=0;
      };
      j++;
    };
/*    if ((i%1000)==0) fprintf(stderr,"%d %d\n",i,Ns);*/
  };

  l=1;
  while(l>0){
/*  px and py are the coordenates of the points of the yardstick */
    px=x1[0]; py=y1[0];

/*  r is the distance between two consecutive points of the yardstick */
    r=0e0; i=0;

/*  (l+1) is the number of sticks used to "cover" the object */
    l=0;

    while(i<Ns){
      r=0e0;
      /* Find the segment where the yardstick touches the object*/
      while((r<stick)&&(i<Ns)){
        i++;
        r=distance(px,py,x2[i],y2[i]);
      };

      if(i<Ns){
        t1=0e0; s1x=x1[i]; s1y=y1[i];
        t2=1e0; s2x=x2[i]; s2y=y2[i];
        epsilon=1.0;

       /* Find where in the segment yardstick touch the object */
       while(epsilon>1e-8){ /* (bissection method) <- This can be improved */
          t=(t1+t2)/2.0;
          sx=s1x+(s2x-s1x)*t;
          sy=s1y+(s2y-s1y)*t;
          r=distance(px,py,sx,sy);
          if(r<stick){
            t1=t; s1x=sx; s1y=sy;
          }else{
            t2=t; s2x=sx; s2y=sy;
          };
          epsilon=fabs(r-stick);
        };

        /* at this point: r ~ stick */
        px=(s1x+s2x)/2.0; py=(s1y+s2y)/2.0;
				keys.push_back(Key(py, px));

        l++;
      }
			else {
				l += eucl_dist(Key(py, px), last);
			};
    };
    printf("%e %lf\n",stick,((double)(l+1))*stick);
    stick*=2e0;
    break;
  };
	return keys;
}

std::map<double, double> Hull::yardstick_rilder(std::vector<double> radii) {
	std::map<double, double> data;

	for (auto const& stick: radii) {
		double *x1,*y1,*x2,*y2,px,py,r;
	  double sx,sy,s1x,s1y,s2x,s2y,t,t1,t2,epsilon;
		double l;
	  long Ns;
		long i,j,ok,g;

		Key first(0, 0), last(0, 0), current(0, 0);
		for (auto const& [k, exists]: hull_points) if (exists) {
			if (k.i == 0) first = k;
			else if (k.i == grid_size-1) last = k;
		}

		std::vector<Key> keys;
		keys.push_back(ordered_hull[0]);

		Ns = ordered_hull.size()-1;
		x1=(double*)calloc(Ns,sizeof(double));
	  y1=(double*)calloc(Ns,sizeof(double));
	  x2=(double*)calloc(Ns,sizeof(double));
	  y2=(double*)calloc(Ns,sizeof(double));

		for (int i = 0; i < Ns; i++) {
			x1[i] = ordered_hull[i].j;
			y1[i] = ordered_hull[i].i;
			x2[i] = ordered_hull[i+1].j;
			y2[i] = ordered_hull[i+1].i;
		}
		/* Main Loop: ------------------------------------------------------------*/

	/* Check if the segments of the object are in the right orientation */
	  for(i=1;i<Ns;i++){
	    j=i; ok=1;
	    while((j<Ns)&&(ok==1)){
	      if((fabs(x1[j]-x2[i-1])<1e-8*stick)&&(fabs(y1[j]-y2[i-1])<1e-8*stick)){
	        px=x1[j]; x1[j]=x1[i]; x1[i]=px;
	        py=y1[j]; y1[j]=y1[i]; y1[i]=py;
	        px=x2[j]; x2[j]=x2[i]; x2[i]=px;
	        py=y2[j]; y2[j]=y2[i]; y2[i]=py;
	        ok=0;
	      };
	      j++;
	    };
	/*    if ((i%1000)==0) fprintf(stderr,"%d %d\n",i,Ns);*/
	  };

	  l=1;
	  while(l>0){
	/*  px and py are the coordenates of the points of the yardstick */
	    px=x1[0]; py=y1[0];

	/*  r is the distance between two consecutive points of the yardstick */
	    r=0e0; i=0;

	/*  (l+1) is the number of sticks used to "cover" the object */
	    l=0;

	    while(i<Ns){
	      r=0e0;
	      /* Find the segment where the yardstick touches the object*/
	      while((r<stick)&&(i<Ns)){
	        i++;
	        r=distance(px,py,x2[i],y2[i]);
	      };

	      if(i<Ns){
	        t1=0e0; s1x=x1[i]; s1y=y1[i];
	        t2=1e0; s2x=x2[i]; s2y=y2[i];
	        epsilon=1.0;

	       /* Find where in the segment yardstick touch the object */
	       while(epsilon>1e-8){ /* (bissection method) <- This can be improved */
	          t=(t1+t2)/2.0;
	          sx=s1x+(s2x-s1x)*t;
	          sy=s1y+(s2y-s1y)*t;
	          r=distance(px,py,sx,sy);
	          if(r<stick){
	            t1=t; s1x=sx; s1y=sy;
	          }else{
	            t2=t; s2x=sx; s2y=sy;
	          };
	          epsilon=fabs(r-stick);
	        };

	        /* at this point: r ~ stick */
	        px=(s1x+s2x)/2.0; py=(s1y+s2y)/2.0;
					keys.push_back(Key(py, px));

	        l++;
	      }
				else {
					l += eucl_dist(Key(py, px), last);
				};
	    };
	    data[stick] = (double)(l+1)*stick;
	    break;
	  };
	}

	return data;
}


std::vector<Key> Hull::test_yardstick_h(float r) {
	std::vector<Key> index_list;

	if (hull.empty()) determine_hull();
	if (linked_hull.empty()) link_hull();

	Key first(0, 0), last(0, 0), current(0, 0);
	for (auto const& [k, exists]: hull_points) if (exists) {
		if (k.i == 0) first = k;
		else if (k.i == grid_size-1) last = k;
	}

	std::map<Bond, float> candidate_crossings; // maps the crossed bond and the smallest graph length of its two keys
	std::map<Key, bool> visited;
	float total_distance = 0;
	int index = 0;
	for (auto const& [k, exists]: hull_points) if (exists) visited[k] = false;
	current = first;

	// we first make a loop that treats the hull (i.e. takes the last crossing on circles)
	do {
		candidate_crossings.clear();
		index_list.push_back(current);
		// find all the spots where the circle crosses a bond
		for (auto const& b: linked_hull) {
			bool crossed = crossing_link(b, current, r);
			if ((!visited[b.a] || !visited[b.b]) && crossed) {
				candidate_crossings[b] = std::min(ordered_hull_keys[b.a], ordered_hull_keys[b.b]);
			}
		}
		// if we get no crossings, that normally means we've reached the end of the graph
		if (candidate_crossings.empty() && eucl_dist(current, last) < r) total_distance += eucl_dist(current, last);
		else if (candidate_crossings.empty() && eucl_dist(current, last) > r) {
			std::cout << "WTF" << r << " :: " << current << std::endl;
			display_yardstick_h("WTFyardstick.png", r);
		}
		else {
			Bond crossed_bond = candidate_crossings.begin()->first;
			float min_length = candidate_crossings.begin()->second;
			Key temp_current(0, 0);

			for (auto const& [b, length]: candidate_crossings) {
				if (length <= min_length) {
					min_length = length, crossed_bond = b;
					temp_current = crossing_point(b, current, r);
				}
			}

			if (temp_current.i < 0 || temp_current.j < 0) std::cout << "WTFF" << r << " :: " << current << std::endl;
			total_distance += r, current = temp_current;
			for (auto const& [k, length]: ordered_hull_keys) {
				if (length < candidate_crossings[crossed_bond]) visited[k] = true;
			}
		}
	} while (!candidate_crossings.empty());
 	std::cout << "Calculated length for hull, radius " << r <<  " yardstick: " << total_distance << std::endl;
	index_list.push_back(last);
	return index_list;
}

std::vector<Key> Hull::test_yardstick_interh(float r) {
	std::vector<Key> index_list;

	if (hull.empty()) determine_hull();
	if (linked_hull.empty()) link_hull();

	Key first(0, 0), last(0, 0), current(0, 0);
	for (auto const& [k, exists]: hull_points) if (exists) {
		if (k.i == 0) first = k;
		else if (k.i == grid_size-1) last = k;
	}

	std::map<Bond, float> candidate_crossings; // maps the crossed bond and the smallest graph length of its two keys
	std::map<Key, bool> visited;
	float total_distance = 0;
	int index = 0;
	for (auto const& [k, exists]: hull_points) if (exists) visited[k] = false;
	current = first;


	// we first make a loop that treats the hull (i.e. takes the last crossing on circles)
	do {
		candidate_crossings.clear();
		index_list.push_back(current);
		// find all the spots where the circle crosses a bond
		for (auto const& b: linked_hull) {
			bool crossed = intercrossing_link(b, current, r);
			if ((!visited[b.a] || !visited[b.b]) && crossed) {
				candidate_crossings[b] = std::min(ordered_hull_keys[b.a], ordered_hull_keys[b.b]);
			}
		}
		// if we get no crossings, that normally means we've reached the end of the graph
		if (candidate_crossings.empty() && eucl_dist(current, last) < r) total_distance += eucl_dist(current, last);
		else if (candidate_crossings.empty() && eucl_dist(current, last) > r) {
			std::cout << "WTF" << r << " :: " << current << std::endl;
			display_yardstick_h("WTFyardstick.png", r);
		}
		else {
			Bond crossed_bond = candidate_crossings.begin()->first;
			float min_length = candidate_crossings.begin()->second;
			Key temp_current(0, 0);

			for (auto const& [b, length]: candidate_crossings) {
				if (length <= min_length) {
					min_length = length, crossed_bond = b;
					temp_current = intercrossing_point(b, current, r);
				}
			}

			if (temp_current.i < 0 || temp_current.j < 0) std::cout << "WTFF" << r << " :: " << current << std::endl;
			total_distance += r, current = temp_current;
			for (auto const& [k, length]: ordered_hull_keys) {
				if (length < candidate_crossings[crossed_bond]) visited[k] = true;
			}
		}
	} while (!candidate_crossings.empty());
 	std::cout << "Calculated length for hull, radius " << r <<  " yardstick: " << total_distance << std::endl;
	index_list.push_back(last);
	return index_list;
}


std::vector<Key> Hull::test_yardstick_ap(float r) {
	std::vector<Key> index_list;

	if (hull.empty()) determine_hull();
	if (linked_hull.empty()) link_hull();

	Key first(0, 0), last(0, 0), current(0, 0);
	for (auto const& [k, exists]: hull_points) if (exists) {
		if (k.i == 0) first = k;
		else if (k.i == grid_size-1) last = k;
	}

	std::map<Bond, float> candidate_crossings; // maps the crossed bond and the smallest graph length of its two keys
	std::map<Key, bool> visited;
	float total_distance = 0;
	int index = 0;
	for (auto const& [k, exists]: hull_points) if (exists) visited[k] = false;
	current = first;

	current = first;

	// we first make a loop that treats the hull (i.e. takes the last crossing on circles)
	do {
		candidate_crossings.clear();
		index_list.push_back(current);
		// find all the spots where the circle crosses a bond
		for (auto const& b: linked_hull) {
			int crossed = crossing_link(b, current, r);
			if ((!visited[b.a] || !visited[b.b]) && crossed) {
				candidate_crossings[b] = std::min(ordered_hull_keys[b.a], ordered_hull_keys[b.b]);
			}
		}
		// if we get no crossings, that means we've reached the end of the graph
		if (candidate_crossings.empty() && eucl_dist(current, last) < r) total_distance += eucl_dist(current, last);
		else if (candidate_crossings.empty() && eucl_dist(current, last) > r) {
			std::cout << "WTF" << r << " :: " << current << std::endl;
			//display_yardstick("WTFyardstick.png", r);
		}
		else {
			Bond crossed_bond = candidate_crossings.begin()->first;
			float max_length = candidate_crossings.begin()->second;
			Key temp_current(0, 0);

			for (auto const& [b, length]: candidate_crossings) {
				if (length >= max_length) {
					max_length = length, temp_current = crossing_point(b, current, r), crossed_bond = b;
				}
			}
			if (temp_current.i < 0 || temp_current.j < 0) std::cout << "WTFF" << r << " :: " << current << std::endl;
			total_distance += r, current = temp_current;
			for (auto const& [k, length]: ordered_hull_keys) {
				if (length < candidate_crossings[crossed_bond]) visited[k] = true;
			}
		}
	} while (!candidate_crossings.empty());

	std::cout << "Calculated length for AP, radius " << r <<  " yardstick: " << total_distance << std::endl;
	index_list.push_back(last);
	return index_list;
}


void Hull::display_hull(std::string name) {
	std::vector<Bond> hull = determine_hull();

	std::vector<Bond> cluster;
	for (auto const& b: g.bond_list) {
		if (g.bond_parameters[b] == 1 && cluster_map[b.a] == 1 && cluster_map[b.b] == 1) {
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

	cairo_set_line_width(cr, 12);

	for (auto const& b: cluster) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(6));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	cairo_set_line_width(cr, 2);

	for (auto const& b: hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(7));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}


	cairo_surface_write_to_png(surface, name.c_str());         // Output to PNG
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void Hull::display_ap(std::string name) {
	std::vector<Bond> hull = determine_ap();

	std::vector<Bond> cluster;
	for (auto const& b: g.bond_list) {
		if (g.bond_parameters[b] == 1 && cluster_map[b.a] == 1 && cluster_map[b.b] == 1) {
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

	cairo_set_line_width(cr, 12);

	for (auto const& b: cluster) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(6));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	cairo_set_line_width(cr, 2);

	for (auto const& b: hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(7));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}


	cairo_surface_write_to_png(surface, name.c_str());         // Output to PNG
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void Hull::display_linked_ap(std::string name, bool show_clust) {
	if (ap.empty()) determine_ap();
	std::map<Key, bool> ap_points;
	for (auto const& b: ap) {
		if (!dense_cluster_map[b.a]) ap_points[b.a] = true;
		else if (!dense_cluster_map[b.b]) ap_points[b.b] = true;
	}

	std::vector<Bond> cluster;
	for (auto const& b: g.bond_list) {
		if (dense_cluster_bond_parameters[b] == 1 && dense_cluster_map[b.a] == 1 && dense_cluster_map[b.b] == 1) {
			cluster.push_back(b);
		}
	}

	if (linked_ap.empty()) link_ap(); //link_ap();

	int WIDTH = 8096;
	int HEIGHT = 8096;

	int u = HEIGHT/grid_size;
	int v = WIDTH/grid_size;

	cairo_surface_t *surface;
	cairo_t *cr;

	surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, HEIGHT, WIDTH);
	cr = cairo_create(surface);

	cairo_set_line_width(cr, 50);

	if (show_clust) for (auto const& b: cluster) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(8));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	cairo_set_line_width(cr, 2);
/*
	for (auto const& b: ap) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(7));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	for (auto const& [k, b]: ap_points) if (b) {
			std::vector<int> rgb(color_hull(9));

			cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
			cairo_move_to(cr, (k.j + 0.5)*u, (k.i + 0.5)*v);
			cairo_rectangle (cr, (k.j + 0.5)*u, (k.i + 0.5)*v, 5, 5);
			cairo_fill(cr);
			cairo_stroke(cr);
		}

	for (auto const& b: linked_ap) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(6));
		//std::cout << Bond(k1, k2) << std::endl;

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}
*/
	if (show_clust) for (auto const& b: ap) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(4));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	if (show_clust) for (auto const& [k, b]: ap_points) if (b) {
			std::vector<int> rgb(color_hull(9));

			cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
			cairo_move_to(cr, (k.j + 0.5)*u, (k.i + 0.5)*v);
			cairo_rectangle (cr, (k.j + 0.5)*u, (k.i + 0.5)*v, 5, 5);
			cairo_fill(cr);
			cairo_stroke(cr);
		}

	for (auto const& b: linked_ap) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(6));
		//std::cout << Bond(k1, k2) << std::endl;

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	cairo_surface_write_to_png(surface, name.c_str());         // Output to PNG
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void Hull::display_linked_hull(std::string name, bool show_clust) {
	if (hull.empty()) determine_hull();
	std::map<Key, bool> hull_points;
	for (auto const& b: hull) {
		if (!cluster_map[b.a]) hull_points[b.a] = true;
		else if (!cluster_map[b.b]) hull_points[b.b] = true;
	}

	std::vector<Bond> cluster;
	for (auto const& b: g.bond_list) {
		if (g.bond_parameters[b] == 1 && cluster_map[b.a] == 1 && cluster_map[b.b] == 1) {
			cluster.push_back(b);
		}
	}

	if (linked_hull.empty()) link_hull();

	int WIDTH = 8096;
	int HEIGHT = 8096;

	int u = HEIGHT/grid_size;
	int v = WIDTH/grid_size;

	cairo_surface_t *surface;
	cairo_t *cr;

	surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, HEIGHT, WIDTH);
	cr = cairo_create(surface);

	cairo_set_line_width(cr, 50);

	if (show_clust) for (auto const& b: cluster) {
	Key k1 = b.a;
	Key k2 = b.b;
	std::vector<int> rgb(color_hull(8));

	cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
	cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
	cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
	cairo_stroke(cr);
}


	cairo_set_line_width(cr, 2);

	if (show_clust) for (auto const& b: hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(4));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	if (show_clust) for (auto const& [k, b]: hull_points) if (b) {
			std::vector<int> rgb(color_hull(9));

			cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
			cairo_move_to(cr, (k.j + 0.5)*u, (k.i + 0.5)*v);
			cairo_rectangle (cr, (k.j + 0.5)*u, (k.i + 0.5)*v, 5, 5);
			cairo_fill(cr);
			cairo_stroke(cr);
		}

	for (auto const& b: linked_hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(6));
		//std::cout << Bond(k1, k2) << std::endl;

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	cairo_surface_write_to_png(surface, name.c_str());         // Output to PNG
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void Hull::display_yardstick_h(std::string name, float radius) {
	if (hull.empty()) determine_hull();
	std::map<Key, bool> hull_points;
	for (auto const& b: hull) {
		if (!cluster_map[b.a]) hull_points[b.a] = true;
		else if (!cluster_map[b.b]) hull_points[b.b] = true;
	}

	std::vector<Bond> cluster;
	for (auto const& b: g.bond_list) {
		if (g.bond_parameters[b] == 1 && cluster_map[b.a] == 1 && cluster_map[b.b] == 1) {
			cluster.push_back(b);
		}
	}

	if (linked_hull.empty()) link_hull();
	std::vector<Key> test = test_yardstick_h(radius);

	int WIDTH = 8192;
	int HEIGHT = 8192;

	int u = HEIGHT/grid_size;
	int v = WIDTH/grid_size;

	cairo_surface_t *surface;
	cairo_t *cr;

	surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, HEIGHT, WIDTH);
	cr = cairo_create(surface);

	cairo_set_line_width(cr, 2);
/*
	for (auto const& b: cluster) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(8));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	cairo_set_line_width(cr, 6);

	for (auto const& b: hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(7));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	for (auto const& [k, b]: hull_points) if (b) {
			std::vector<int> rgb(color_hull(9));

			cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
			cairo_move_to(cr, (k.j + 0.5)*u, (k.i + 0.5)*v);
			cairo_rectangle (cr, (k.j + 0.5)*u, (k.i + 0.5)*v, 5, 5);
			cairo_fill(cr);
			cairo_stroke(cr);
		}
*/
	for (auto const& b: linked_hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(6));
		//std::cout << Bond(k1, k2) << std::endl;

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	for (int i = 0; i < test.size()-1; i++) {
		Key k1 = test[i], k2 = test[i+1];
		cairo_set_line_width(cr, 4);

		cairo_select_font_face(cr, "Purisa",
			CAIRO_FONT_SLANT_NORMAL,
			CAIRO_FONT_WEIGHT_BOLD);

		std::vector<int> rgb(color_hull(4));
		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);

		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		//cairo_rectangle (cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v, 10, 10);
		//cairo_set_font_size(cr, 100);
		//cairo_show_text(cr, std::to_string(i).c_str());
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}


	cairo_surface_write_to_png(surface, name.c_str());         // Output to PNG
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void Hull::display_yardstick_ap(std::string name, float radius) {
	if (hull.empty()) determine_hull();
	std::map<Key, bool> hull_points;
	for (auto const& b: hull) {
		if (!cluster_map[b.a]) hull_points[b.a] = true;
		else if (!cluster_map[b.b]) hull_points[b.b] = true;
	}

	std::vector<Bond> cluster;
	for (auto const& b: g.bond_list) {
		if (g.bond_parameters[b] == 1 && cluster_map[b.a] == 1 && cluster_map[b.b] == 1) {
			cluster.push_back(b);
		}
	}

	if (linked_hull.empty()) link_hull();
	std::vector<Key> test = test_yardstick_ap(radius);

	int WIDTH = 8192;
	int HEIGHT = 8192;

	int u = HEIGHT/grid_size;
	int v = WIDTH/grid_size;

	cairo_surface_t *surface;
	cairo_t *cr;

	surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, HEIGHT, WIDTH);
	cr = cairo_create(surface);

	cairo_set_line_width(cr, 2);
/*
	for (auto const& b: cluster) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(8));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	cairo_set_line_width(cr, 4);

	for (auto const& b: hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(7));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	for (auto const& [k, b]: hull_points) if (b) {
			std::vector<int> rgb(color_hull(9));

			cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
			cairo_move_to(cr, (k.j + 0.5)*u, (k.i + 0.5)*v);
			cairo_rectangle (cr, (k.j + 0.5)*u, (k.i + 0.5)*v, 5, 5);
			cairo_fill(cr);
			cairo_stroke(cr);
		}
*/
	for (auto const& b: linked_hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(6));
		//std::cout << Bond(k1, k2) << std::endl;

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	for (int i = 0; i < test.size()-1; i++) {
		Key k1 = test[i], k2 = test[i+1];
		cairo_set_line_width(cr, 4);

		cairo_select_font_face(cr, "Purisa",
      CAIRO_FONT_SLANT_NORMAL,
      CAIRO_FONT_WEIGHT_BOLD);

		std::vector<int> rgb(color_hull(4));
		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);

		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		//cairo_rectangle (cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v, 10, 10);
		//cairo_set_font_size(cr, 100);
		//cairo_show_text(cr, std::to_string(i).c_str());
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	cairo_surface_write_to_png(surface, name.c_str());         // Output to PNG
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}

void Hull::display_yardstick_rilder(std::string name, float radius) {
	if (hull.empty()) determine_hull();
	std::map<Key, bool> hull_points;
	for (auto const& b: hull) {
		if (!cluster_map[b.a]) hull_points[b.a] = true;
		else if (!cluster_map[b.b]) hull_points[b.b] = true;
	}

	std::vector<Bond> cluster;
	for (auto const& b: g.bond_list) {
		if (g.bond_parameters[b] == 1 && cluster_map[b.a] == 1 && cluster_map[b.b] == 1) {
			cluster.push_back(b);
		}
	}

	if (linked_hull.empty()) link_hull();
	std::vector<Key> test = yardstick_rilder_test(radius);

	int WIDTH = 8192;
	int HEIGHT = 8192;

	int u = HEIGHT/grid_size;
	int v = WIDTH/grid_size;

	cairo_surface_t *surface;
	cairo_t *cr;

	surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, HEIGHT, WIDTH);
	cr = cairo_create(surface);

	cairo_set_line_width(cr, 2);
/*
	for (auto const& b: cluster) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(8));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	cairo_set_line_width(cr, 6);

	for (auto const& b: hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(7));

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	for (auto const& [k, b]: hull_points) if (b) {
			std::vector<int> rgb(color_hull(9));

			cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
			cairo_move_to(cr, (k.j + 0.5)*u, (k.i + 0.5)*v);
			cairo_rectangle (cr, (k.j + 0.5)*u, (k.i + 0.5)*v, 5, 5);
			cairo_fill(cr);
			cairo_stroke(cr);
		}
*/
	for (auto const& b: linked_hull) {
		Key k1 = b.a;
		Key k2 = b.b;
		std::vector<int> rgb(color_hull(6));
		//std::cout << Bond(k1, k2) << std::endl;

		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);
		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}

	for (int i = 0; i < test.size()-1; i++) {
		Key k1 = test[i], k2 = test[i+1];
		cairo_set_line_width(cr, 4);

		cairo_select_font_face(cr, "Purisa",
			CAIRO_FONT_SLANT_NORMAL,
			CAIRO_FONT_WEIGHT_BOLD);

		std::vector<int> rgb(color_hull(4));
		cairo_set_source_rgb(cr, rgb[0], rgb[1], rgb[2]);

		cairo_move_to(cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v);
		//cairo_rectangle (cr, (k1.j + 0.5)*u, (k1.i + 0.5)*v, 10, 10);
		//cairo_set_font_size(cr, 100);
		//cairo_show_text(cr, std::to_string(i).c_str());
		cairo_line_to(cr, (k2.j + 0.5)*u, (k2.i + 0.5)*v);
		cairo_stroke(cr);
	}


	cairo_surface_write_to_png(surface, name.c_str());         // Output to PNG
	cairo_destroy(cr);
	cairo_surface_destroy(surface);
}
