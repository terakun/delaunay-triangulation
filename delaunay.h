#pragma once
#include <utility>
#include <vector>
#include <complex>
#include <array>
#include <cassert>

struct vertex {
  std::complex<double> z;
  vertex(const std::complex<double> &z):z(z){}
  vertex(double x,double y):z(x,y){}
  ~vertex(){}
};

struct face;

struct halfedge {
  vertex *v;
  face *f;
  halfedge *inv;
  halfedge *next;
  halfedge *new_e;
  halfedge(halfedge * const e):v(e->v),f(nullptr),inv(nullptr),next(nullptr),new_e(nullptr) {
    if(e->inv){
      inv = e->inv;
      e->inv->inv = this;
    }
  }
  halfedge(vertex *v,face *f=nullptr,halfedge *inv=nullptr,halfedge *next=nullptr):v(v),f(f),inv(inv),next(next),new_e(nullptr) {}
  ~halfedge() {}
};

struct face {
  halfedge *e;
  face *child[3]; // 探索用
  bool is_leaf; // 探索用
  face(halfedge *e=nullptr):e(e),child{nullptr,nullptr,nullptr},is_leaf(true) {}
  ~face() {}
  std::array<std::complex<double>,3> get_points() const{
    std::array<std::complex<double>,3> zs;
    halfedge *te = e;
    for(int i=0;i<3;++i){
      zs[i] = te->v->z;
      te = te->next;
    }
    return zs;
  }
  bool contains(double x,double y)const {
    auto zs = get_points();
    std::complex<double> z(x,y);

    double cross[3] = {
      (std::conj(zs[1]-zs[0])*(z-zs[0])).imag(),
      (std::conj(zs[2]-zs[1])*(z-zs[1])).imag(),
      (std::conj(zs[0]-zs[2])*(z-zs[2])).imag()
    };

    return (cross[0] > 0)&&(cross[1] > 0)&&(cross[2] > 0);
  }
  bool contains(const vertex *v) const{
    auto zs = get_points();
    std::complex<double> z(v->z);
    std::complex<double> p = 
      ( zs[0]*std::norm(zs[1]) + zs[1]*std::norm(zs[2]) + zs[2]*std::norm(zs[0]) 
      - zs[2]*std::norm(zs[1]) - zs[1]*std::norm(zs[0]) - zs[0]*std::norm(zs[2]) ) /
      ( zs[0]*std::conj(zs[1]) + zs[1]*std::conj(zs[2]) + zs[2]*std::conj(zs[0]) 
      - zs[2]*std::conj(zs[1]) - zs[1]*std::conj(zs[0]) - zs[0]*std::conj(zs[2]) );
    double r2 = std::norm(p - zs[0]);
    return std::norm(z-p) < r2;
  }
};

void delaunay_triangulate(const std::vector<std::pair<double,double>> &points);

