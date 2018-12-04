#include <iostream>
#include <cstdlib>
#include <stack>
#include "./delaunay.h"

void make_triangle(halfedge *e_ab,halfedge *e_bc,halfedge *e_ca,face *f){
  e_ab->next = e_bc;
  e_bc->next = e_ca;
  e_ca->next = e_ab;

  e_ab->f = f;
  e_bc->f = f;
  e_ca->f = f;
}

void flip(face *f_abc,face *f_acd,halfedge **e_ca,halfedge **e_ac) {
  if((*e_ca)->inv != *e_ac){
    std::cout << "something wrong in flip" << std::endl;
    std::abort();
  }

  vertex *v_b = (*e_ca)->next->next->v;
  halfedge *e_bd = new halfedge(v_b);

  vertex *v_d = (*e_ac)->next->next->v;
  halfedge *e_db = new halfedge(v_d);

  e_bd->inv = e_db;
  e_db->inv = e_bd;

  halfedge *e_ab = (*e_ca)->next;
  halfedge *e_bc = (*e_ca)->next->next;
  halfedge *e_cd = (*e_ac)->next;
  halfedge *e_da = (*e_ac)->next->next;

  make_triangle(e_ab,e_bd,e_da,f_abc);
  make_triangle(e_cd,e_db,e_bc,f_acd);

  f_abc->e = e_bd;
  f_acd->e = e_db;
  delete *e_ca;
  *e_ca = e_bd;
  delete *e_ac;
  *e_ac = e_db;
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
  vertex *p1 = new vertex(x1,y1);
  double x2 = maxx + (maxy-miny)/std::sqrt(3) + 2.0;
  double y2 = miny - 1.0;
  vertex *p2 = new vertex(x2,y2);
  double x3 = minx + (maxx-minx)/2;
  double y3 = maxy + (maxx-minx)/2*std::sqrt(3);
  vertex *p3 = new vertex(x3,y3);

  halfedge *e_12 = new halfedge(p1);
  halfedge *e_23 = new halfedge(p2);
  halfedge *e_31 = new halfedge(p3);

  face *f = new face(e_12);

  make_triangle(e_12,e_23,e_31,f);
  return f;
}

void delaunay_triangulate(const std::vector<std::pair<double,double>> &points) {
  std::vector<face*> faces;
  face *st = super_triangle(points);
  std::array<std::complex<double>,3> stz{
    st->e->v->z , st->e->next->v->z , st->e->next->next->v->z
  };
  faces.push_back(st);

  int cnt = 1;
  for(const auto &point:points){
    if((cnt%100) == 0) {
      std::cerr << cnt << std::endl;
    }
    cnt++;

    face *f = nullptr;
    int idx;
    for(int i=0;i<faces.size();++i){
      if(faces[i]->contains(point.first,point.second)){
        idx = i;
        f = faces[i];
        break;
      }
    }
    if(f == nullptr){
      std::cerr << "something wrong" << std::endl;
      std::abort();
    }

    vertex *p = new vertex(point.first,point.second);

    halfedge *e_ab = f->e;
    halfedge *e_bc = e_ab->next;
    halfedge *e_ca = e_bc->next;

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

    e_ap->inv = e_pa;
    e_pa->inv = e_ap;

    e_bp->inv = e_pb;
    e_pb->inv = e_bp;

    e_cp->inv = e_pc;
    e_pc->inv = e_cp;

    delete f;
    f = nullptr;

    faces[idx] = f_abp;
    faces.push_back(f_bcp);
    faces.push_back(f_cap);

    std::stack<halfedge*> edgestack;
    edgestack.push(e_ab);
    edgestack.push(e_bc);
    edgestack.push(e_ca);

    while(!edgestack.empty()){
      halfedge *e1 = edgestack.top();edgestack.pop();
      face *f1 = e1->f;
      if(e1->inv == nullptr) continue;
      halfedge *e2 = e1->inv;
      face *f2 = e2->f;

      if(f1->contains(e2->next->next->v)){
        flip(f1,f2,&e1,&e2);

        halfedge *e_da = e1;
        halfedge *e_ab = e_da->next;
        halfedge *e_bd = e_ab->next;

        halfedge *e_ad = e2;
        halfedge *e_dc = e_ad->next;
        halfedge *e_ca = e_dc->next;

        edgestack.push(e_ab);
        edgestack.push(e_bd);
        edgestack.push(e_dc);
        edgestack.push(e_ca);
      }
    }
  }

  for(int i=0;i<faces.size();++i){
    auto zs = faces[i]->get_points();
    bool has_super = false;
    for(int j=0;j<3;++j){
      for(int k=0;k<3;++k){
        if(zs[j] == stz[k]) {
          has_super = true;
          break;
        }
      }
    }
    if(has_super) continue;
    for(int i=0;i<3;++i) std::cout << zs[i].real() << " " << zs[i].imag() << std::endl;
    std::cout << zs[0].real() << " " << zs[0].imag() << std::endl;
    std::cout << std::endl;
  }
}


