#include <iostream>
#include <algorithm>
#include <cstdlib>
#include <stack>
#include <tuple>
#include <fstream>
#include "./delaunay.h"

void make_triangle(halfedge *e_ab,halfedge *e_bc,halfedge *e_ca,face *f){
  e_ab->next = e_bc;
  e_bc->next = e_ca;
  e_ca->next = e_ab;

  e_ab->f = f;
  e_bc->f = f;
  e_ca->f = f;
}

std::tuple<face*,face*,halfedge*,halfedge*> flip(face *f_cab,face *f_acd,halfedge *e_ca,halfedge *e_ac) {
  if(e_ca->inv != e_ac){
    std::cout << "something wrong in flip" << std::endl;
    std::abort();
  }

  vertex *v_b = e_ca->next->next->v;
  halfedge *e_bd = new halfedge(v_b);

  vertex *v_d = e_ac->next->next->v;
  halfedge *e_db = new halfedge(v_d);

  e_bd->inv = e_db;
  e_db->inv = e_bd;

  halfedge *e_ab = new halfedge(e_ca->next);
  e_ca->next->new_e = e_ab;
  halfedge *e_bc = new halfedge(e_ca->next->next);
  e_ca->next->next->new_e = e_bc;
  halfedge *e_cd = new halfedge(e_ac->next);
  e_ac->next->new_e = e_cd;
  halfedge *e_da = new halfedge(e_ac->next->next);
  e_ac->next->next->new_e = e_da;

  face *f_abd = new face(e_bd);
  e_ca->next->f = f_abd;
  e_ac->next->next->f = f_abd;
  make_triangle(e_ab,e_bd,e_da,f_abd);

  face *f_cdb = new face(e_db);
  e_ac->next->f = f_cdb;
  e_ca->next->next->f = f_cdb;
  make_triangle(e_cd,e_db,e_bc,f_cdb);

  f_cab->child[0] = f_abd;
  f_cab->child[1] = f_cdb;
  f_cab->is_leaf = false;

  f_acd->child[0] = f_abd;
  f_acd->child[1] = f_cdb;
  f_acd->is_leaf = false;

  return std::make_tuple(f_abd,f_cdb,e_bd,e_db);
}


face* super_triangle(const std::vector<std::pair<double,double>> &points) {
  double minx = 1e10,miny = 1e10;
  double maxx = -1e10,maxy = -1e10;
  for (auto p : points) {
    minx = std::min(minx,p.first);
    miny = std::min(miny,p.second);

    maxx = std::max(maxx,p.first);
    maxy = std::max(maxy,p.second);
  }
  double x1 = minx - (maxy-miny)/std::sqrt(3) - 2.0;
  double y1 = miny - 1.0;
  vertex *p1 = new vertex(0,x1,y1);
  double x2 = maxx + (maxy-miny)/std::sqrt(3) + 2.0;
  double y2 = miny - 1.0;
  vertex *p2 = new vertex(0,x2,y2);
  double x3 = minx + (maxx-minx)/2;
  double y3 = maxy + (maxx-minx)/2*std::sqrt(3);
  vertex *p3 = new vertex(0,x3,y3);

  halfedge *e_12 = new halfedge(p1);
  halfedge *e_23 = new halfedge(p2);
  halfedge *e_31 = new halfedge(p3);

  face *f = new face(e_12);

  make_triangle(e_12,e_23,e_31,f);
  return f;
}

void delaunay_triangulate(const std::vector<std::pair<double,double>> &points) {
  face *st = super_triangle(points);
  std::array<std::complex<double>,3> stz{
    st->e->v->z , st->e->next->v->z , st->e->next->next->v->z
  };

  int id = 0;
  for(const auto &point:points){
    if(((id+1)%10000) == 0) {
      std::cerr << id+1 << std::endl;
    }

    face *f = st;
    int d=0;
    while(!f->is_leaf) {
      face *pref = f;
      for(int i=0;i<3;++i){
        if(f->child[i]&&f->child[i]->contains(point.first,point.second)){
          f = f->child[i];
          break;
        }
      }
      if(pref == f){
        std::cerr << "something wrong" << std::endl;
        std::abort();
      }
    }

    vertex *p = new vertex(id++,point.first,point.second);

    halfedge *e_ab = new halfedge(f->e);
    f->e->new_e = e_ab;
    halfedge *e_bc = new halfedge(f->e->next);
    f->e->next->new_e = e_bc;
    halfedge *e_ca = new halfedge(f->e->next->next);
    f->e->next->next->new_e = e_ca;

    face *f_abp = new face(e_ab);
    halfedge *e_bp = new halfedge(e_bc->v);
    halfedge *e_pa = new halfedge(p);
    make_triangle(e_ab,e_bp,e_pa,f_abp);

    face *f_bcp = new face(e_bc);
    halfedge *e_cp = new halfedge(e_ca->v);
    halfedge *e_pb = new halfedge(p);
    make_triangle(e_bc,e_cp,e_pb,f_bcp);

    face *f_cap = new face(e_ca);
    halfedge *e_ap = new halfedge(e_ab->v);
    halfedge *e_pc = new halfedge(p);
    make_triangle(e_ca,e_ap,e_pc,f_cap);

    f->child[0] = f_abp;
    f->child[1] = f_bcp;
    f->child[2] = f_cap;
    f->is_leaf = false;

    e_ap->inv = e_pa;
    e_pa->inv = e_ap;

    e_bp->inv = e_pb;
    e_pb->inv = e_bp;

    e_cp->inv = e_pc;
    e_pc->inv = e_cp;

    std::vector<halfedge*> edgestack;
    edgestack.push_back(e_ab);
    edgestack.push_back(e_bc);
    edgestack.push_back(e_ca);

    while(!edgestack.empty()){
      halfedge *e_ca = edgestack.back();edgestack.pop_back();
      while(e_ca->new_e){
        e_ca = e_ca->new_e;
      }
      face *f_cab = e_ca->f;

      if(e_ca->inv == nullptr) continue;
      halfedge *e_ac = e_ca->inv;

      face *f_acd = e_ac->f;

      if(f_cab->contains(e_ac->next->next->v)){
        face *f_abd,*f_cdb;
        halfedge *e_bd,*e_db;
        std::tie(f_abd,f_cdb,e_bd,e_db) = flip(f_cab,f_acd,e_ca,e_ac);

        halfedge *e_da = e_bd->next;
        halfedge *e_ab = e_da->next;

        halfedge *e_bc = e_db->next;
        halfedge *e_cd = e_bc->next;

        edgestack.push_back(e_da);
        edgestack.push_back(e_ab);
        edgestack.push_back(e_bc);
        edgestack.push_back(e_cd);
      }
    }
  }

  std::ofstream ofs_plot("plot.txt");
  std::vector<std::vector<size_t>> graph(points.size());
  std::vector<face*> s;
  s.push_back(st);
  while(!s.empty()){
    face* f = s.back(); s.pop_back();
    if(!f->is_leaf){
      for(int i=0;i<3;++i){
        if(f->child[i]){
          s.push_back(f->child[i]);
        }
      }
    }else{
      f->is_leaf = false;
      auto zs = f->get_points();
      bool has_super = false;
      for(int j=0;j<3;++j){
        for(int k=0;k<3;++k){
          if(zs[j] == stz[k]){
            has_super = true;
            break;
          }
        }
      }
      if(has_super) continue;
      halfedge *e = f->e;
      for(int i=0;i<3;++i){
        graph[e->v->id].push_back(e->next->v->id);
        graph[e->next->v->id].push_back(e->v->id);
        e = e->next;
      }
      for(int i=0;i<3;++i) ofs_plot << zs[i].real() << " " << zs[i].imag() << std::endl;
      ofs_plot << zs[0].real() << " " << zs[0].imag() << std::endl;
      ofs_plot << std::endl;

    }
  }

  std::ofstream ofs_graph("graph.txt");
  ofs_graph << graph.size() << std::endl;
  for(auto &&e:graph){
    std::sort(std::begin(e),std::end(e));
    e.erase(std::unique(std::begin(e), std::end(e)), std::end(e));
    ofs_graph << e.size();
    for(auto v:e){
      ofs_graph << " " << v;
    }
    ofs_graph << std::endl;
  }
}

