#include <iostream>
#include <fstream>
#include <string>
#include <random>

#include <parlay/primitives.h>
#include <parlay/sequence.h>
#include <parlay/random.h>

#include "convex_hull_3d.h"

// **************************************************************
// Driver
// **************************************************************
int main(int argc, char* argv[]) {
  auto usage = "Usage: <n>";
  if (argc != 2) std::cout << usage << std::endl;
  else {
    point_id n;
    try {n = std::stoi(argv[1]);}
    catch (...) {std::cout << usage << std::endl; return 1;}
    parlay::random_generator gen(0);
    std::uniform_real_distribution<real> dis(0.0,1.0);

    // generate n random points in a unit square
    auto points = parlay::tabulate(n, [&] (point_id i) -> point {
      auto r = gen[i];
      return point{i, dis(r), dis(r), dis(r)};});

    parlay::internal::timer t("Time");
    parlay::sequence<tri> result;
    for (int i=0; i < 5; i++) {
      result = convex_hull_3d(points);
      t.next("convex_hull_3d");
    }

    std::cout << parlay::num_workers() << std::endl;

    // print all points info
    std::ofstream inFile("convex_hull.in");
    if (!inFile.is_open()) {
      std::cerr << "Failed to open output file!" << std::endl;
    } else {
      for (const auto& p : points) {
        inFile << p.x << " " << p.y << " " << p.z << std::endl;
      }
      inFile.close();
    }    

    // print convex_hull
    std::ofstream outFile("convex_hull.out");
    if (!outFile.is_open()) {
      std::cerr << "Failed to open output file!" << std::endl;
    } else {
      for (const auto& tri : result) {
        point_id id = tri[0];
        outFile << points[id].x << " " << points[id].y << " " << points[id].z;

        id = tri[1];
        outFile << " " << points[id].x << " " << points[id].y << " " << points[id].z;

        id = tri[2];
        outFile << " " << points[id].x << " " << points[id].y << " " << points[id].z << std::endl;
      }
      outFile.close();
    }    

    std::cout << "number of triangles in the mesh = " << result.size() << std::endl;
  }
}
