#include <iostream>
#include <fstream>
#include <vector>
#include <utility>
#include "./delaunay.h"

int main(int argc,char **argv){
  if(argc < 2){
    std::cerr << argv[0] << " [filename]" << std::endl;
    return 1;
  }
  std::ifstream ifs(argv[1]);
  std::vector<std::pair<double,double>> points;
  int id;
  double x,y;
  while(ifs >> id >> x >> y){
    points.emplace_back(x,y);
  }
  std::cerr << "triangulate" << std::endl;
  delaunay_triangulate(points);
}
