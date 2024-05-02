#include <utility>
#include <memory>
#include <optional>

#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/utilities.h>

#include "convex_hash_map.h"

std::atomic_flag print_lock = ATOMIC_FLAG_INIT;

void atomic_print(const std::string& message) {
    while (print_lock.test_and_set(std::memory_order_acquire)); // lock
    std::cout << message << std::endl;
    print_lock.clear(std::memory_order_release); // unlock
}

// **************************************************************
// Parallel Convel Hull in 3D
// From the paper:
// Randomized Incremental Convex Hull is Highly Parallel
// Blelloch, Gu, Shun and Sun
// **************************************************************

using real = float;
using point_id = int;

using tri = std::array<point_id,3>;
using edge = std::array<point_id,2>;

struct point {
  point_id id; real x; real y; real z;
  const bool operator<(const point p) const {return id < p.id;}
};

struct vect {
  double x; double y; double z;
  vect(double x, double y, double z) : x(x), y(y), z(z) {}
};

struct triangle {
  tri t;
  point_id pid; // point in convex hull not in triangle
  parlay::sequence<point> conflicts;
  triangle();
  triangle(tri t, point_id pid, parlay::sequence<point> c) : t(t), pid(pid), conflicts(c) {}
};

inline auto get_normal_vect(const point& a, const point& b, const point& c) {
  return vect{(b.y - a.y) * (c.z - a.z) - (b.z - a.z) * (c.y - a.y),
              (b.z - a.z) * (c.x - a.x) - (b.x - a.x) * (c.z - a.z),
              (b.x - a.x) * (c.y - a.y) - (b.y - a.y) * (c.x - a.x)};
}

// check above to define visible points
inline auto is_above(const point& a, const vect& normal, const point& target) {
  auto dot = [=] (const vect& v1, const vect& v2) {
    return v1.x * v2.x + v1.y * v2.y + v1.z * v2.z;};
  double dist = dot(vect(a.x - target.x, a.y - target.y, a.z - target.z), normal);
  return dist >= 0;
}

// **************************************************************
// The main body
// **************************************************************

using Points = parlay::sequence<point>;

struct Convex_Hull_3d {
  using triangle_ptr = std::shared_ptr<triangle>;
  convex_hash_map<edge, triangle_ptr> map_facets;
  Points points;
  convex_hash_map<tri, bool> convex_hull;
  point_id n;

  point_id min_conflicts(triangle_ptr& t) {
    return (t->conflicts.size() == 0) ? n : t->conflicts[0].id;}

  auto get_visible_points(tri t, point_id pid, Points p) {
    auto normal = get_normal_vect(points[t[0]], points[t[1]], points[t[2]]);
    auto is_convex_above = is_above(points[t[0]], normal, points[pid]);
    auto keep = parlay::tabulate(p.size(), [&] (long i) {
      return (is_convex_above != is_above(points[t[0]], normal, p[i]));
    });
    return parlay::pack(p, keep);
  }

  // recursive routine
  // convex hull in 3d space
  void process_ridge(triangle_ptr& t1, edge r, triangle_ptr& t2) {
    if (t1->conflicts.size() == 0 && t2->conflicts.size() == 0) {
      return;
    } else if (min_conflicts(t2) == min_conflicts(t1)) {
      // H \ {t1, t2}
      convex_hull.remove(t1->t);
      convex_hull.remove(t2->t);
    } else if (min_conflicts(t2) < min_conflicts(t1)) {
      process_ridge(t2, r, t1);
    } else {
      point_id pid = min_conflicts(t1);
      tri t{r[0], r[1], pid};

      // C(t) <- v' in C(t1) U C(t2) | visible(v', t)
      auto uni = parlay::merge(t1->conflicts, t2->conflicts);
      auto keep = parlay::tabulate(uni.size(), [&] (long i) {
        return ((i != 0) && (uni[i].id != uni[i-1].id));});
      auto uni_dedup = parlay::pack(uni, keep);

      auto t_new = std::make_shared<triangle>(t,
                                              t1->pid, 
                                              get_visible_points(t, 
                                                                 t1->pid, 
                                                                 uni_dedup));

      // H <- (H \ {t1}) U {t}
      convex_hull.remove(t1->t);
      convex_hull.insert_and_set(t, true);

      auto check_edge = [&] (edge e, triangle_ptr& tp) {
        auto key = (e[0] < e[1]) ? e : edge{e[1], e[0]};
        if (map_facets.insert_and_set(key, tp)) return;

        auto tt = map_facets.get_value(e, tp).value();
        process_ridge(tp, e, tt);};
      parlay::par_do3([&] {process_ridge(t_new, r, t2);},
                      [&] {check_edge(edge{r[0], pid}, t_new);},
                      [&] {check_edge(edge{r[1], pid}, t_new);});
    }
  }

  // toplevel routine
  // The result is set of facets: result_facets
  // assum that P has more than 4 points
  Convex_Hull_3d(const Points& P) :
      map_facets(convex_hash_map<edge, triangle_ptr>(6*P.size())),
      convex_hull(convex_hash_map<tri, bool>(6*P.size())),
      n(P.size()) {
    points = P;
    tri init_tri[4] = {{0, 1, 2},
                       {1, 2, 3},
                       {0, 2, 3},
                       {0, 1, 3}};
    point_id remain[4] = {3, 0, 1, 2};

    // get initial convex hull
    convex_hull.insert_and_set(init_tri[0], true);
    convex_hull.insert_and_set(init_tri[1], true);
    convex_hull.insert_and_set(init_tri[2], true);
    convex_hull.insert_and_set(init_tri[3], true);

    auto keep = parlay::tabulate(P.size(), [&] (long i) {
      return (i > 3);});
    Points target_points = parlay::pack(P, keep);

    // get visible foreach t
    triangle_ptr t[4];
    parlay::parallel_for(0, 4, [&](auto i) {
      t[i] = std::make_shared<triangle>(init_tri[i],
                                        remain[i], 
                                        get_visible_points(init_tri[i], 
                                                           remain[i], 
                                                           target_points));
    });

    struct bound { int t1_idx; int t2_idx; edge e;};
    bound share_info[6] = {{0, 1, edge{1, 2}},
                           {0, 2, edge{0, 2}},
                           {0, 3, edge{0, 1}},
                           {1, 2, edge{2, 3}},
                           {1, 3, edge{1, 3}},
                           {2, 3, edge{0, 3}}};
    
    parlay::parallel_for(0, 6, [&](auto i) {
      process_ridge(t[share_info[i].t1_idx],
                    share_info[i].e,
                    t[share_info[i].t2_idx]);
    });    
  }
};

parlay::sequence<tri> convex_hull_3d(const Points& P) {
  Convex_Hull_3d ch3d(P);
  return ch3d.convex_hull.keys();
}
