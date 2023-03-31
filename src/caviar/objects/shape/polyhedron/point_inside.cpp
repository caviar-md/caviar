
//========================================================================
//
// Copyright (C) 2019 by Morad Biagooi and Ehsan Nedaaee Oskoee.
//
// This file is part of the CAVIAR package.
//
// The CAVIAR package is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 3.0 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the CAVIAR distribution.
//
//========================================================================

#include "caviar/objects/shape/polyhedron/point_inside.h"
#include "caviar/objects/shape/polyhedron/polyhedron.h"
#include "caviar/utility/interpreter_io_headers.h"

#include <string>
#include <cmath>
#include <fstream>


namespace caviar {

namespace shape {
namespace polyhedron {

static inline int int_floor(double x) { 
    return (int)(x+100000) - 100000; 
}

Point_Inside::Point_Inside (CAVIAR *fptr) : Pointers{fptr} {
  point_is_inside_method = -1;
}

Point_Inside::~Point_Inside () { }

// this function is used for particle distribution
bool Point_Inside::is_inside(shape::polyhedron::Polyhedron & p_object, const Vector<double> &v0) {
  auto method = point_is_inside_method;
  if (method == -1){
    // since the method of ray is not perfect yet, we use these two additional
    // vector values to add to the rays.
    Vector<double> v3 = {1e-5,3e-5,25e-6};
    Vector<double> v4 = {-1e-5,2e-5,-3e-5};
 
    const auto r0 = ray_tells_point_is_inside (p_object, v0, 0);  
    const auto r1 = ray_tells_point_is_inside (p_object, v0, 1);  
    const auto r2 = ray_tells_point_is_inside (p_object, v0, 2); 
    const auto r3 = ray_tells_point_is_inside (p_object, v0+v3, 0);  
    const auto r4 = ray_tells_point_is_inside (p_object, v0+v4, 1);     
    const int res = r0 + r1 + r2 + r3 + r4;
    
    if (res == 5)
      return true;

  } else if (method == 0) {
    return ray_tells_point_is_inside (p_object, v0, 0);
  } else if (method == 1) {
    return ray_tells_point_is_inside (p_object, v0, 1);
  } else if (method == 2) {
    return ray_tells_point_is_inside (p_object, v0, 2);
  } else if (method == 123) {
    const auto r0 = ray_tells_point_is_inside (p_object, v0, 0);  
    const auto r1 = ray_tells_point_is_inside (p_object, v0, 1);  
    const auto r2 = ray_tells_point_is_inside (p_object, v0, 2); 
    const int res = r0 + r1 + r2;
    if (res == 3)
      return true;
  } else {
    error->all(FC_FILE_LINE_FUNC,"Undefined point inside method");
  }

  return false;   
}



bool Point_Inside::is_inside_grid(shape::polyhedron::Polyhedron & p_object, const Vector<double> &v, const double r) {
  if (!is_inside(p_object, v)) return false;
  if (in_contact_grid(p_object, v, r)) return false; 
  return true;
}

bool Point_Inside::is_inside_all(shape::polyhedron::Polyhedron & p_object, const Vector<double> &v, const double r) {
  if (!is_inside(p_object, v)) return false;
  if (in_contact_all(p_object, v, r)) return false; 
  return true;
}


// this is called in polyhedron::Handler->in_contact() if activated.
// the same algorithm as (in_contact_all), except all the polygons are going to be checked
// any fix should be done on both of them.
// for the optimization reasons, we won't make them one functions.
bool Point_Inside::in_contact_grid (shape::polyhedron::Polyhedron & p_object, 
    const Vector<double> &v, const Real_t radius, Vector<double> &contact_vector) {

  
  const auto & vertex = p_object.vertex;
  const auto & face = p_object.face;
  const auto & normal = p_object.normal;
   
  Vector<Real_t> compression_vector {0.0,0.0,0.0};
     
  const auto & grid = p_object.grid;
   
  const auto & xlo = p_object.xlo;
  const auto & ylo = p_object.ylo;
  const auto & zlo = p_object.zlo;
  
  const auto & dx_part = p_object.dx_part;
  const auto & dy_part = p_object.dy_part;
  const auto & dz_part = p_object.dz_part;
  
  const auto & nx_part = p_object.nx_part;
  const auto & ny_part = p_object.ny_part;
  const auto & nz_part = p_object.nz_part;

  int xindex = int_floor((v.x-xlo)/dx_part);
  int yindex = int_floor((v.y-ylo)/dy_part);
  int zindex = int_floor((v.z-zlo)/dz_part);

  if (xindex<0) return false;
  if (yindex<0) return false;
  if (zindex<0) return false;

  if (xindex>int_floor(nx_part-1)) return false;
  if (yindex>int_floor(ny_part-1)) return false;
  if (zindex>int_floor(nz_part-1)) return false;


  bool contact_flag = false;
    
  for (auto i : grid[xindex][yindex][zindex]) {

    // renaming the vertices of the 
    const auto &fv0 = vertex[face[i][0]];
    const auto &fv1 = vertex[face[i][1]];    
    const auto &fv2 = vertex[face[i][2]];
    
    // renaming the normal    
    const auto &n = normal[i];

    // distance from v to the plane containing polygon, using fv0,
    // positve value for dis_v_poly means that v is on the side which normal vector
    // points to.
    const auto dis_v_poly = (v - fv0)*n;
    const auto dis_v_poly_sq = dis_v_poly*dis_v_poly;
    
    // projection of v on the polyhedron
    const auto vpr = v - (dis_v_poly*n);

    // vectors from the projection point to vertices
    const auto vpr_fv0 = fv0 - vpr;
    const auto vpr_fv1 = fv1 - vpr;
    const auto vpr_fv2 = fv2 - vpr;

    // cross productions of vpr_fvX,
    const auto cp_vpr_fv0_fv1 = cross_product(vpr_fv0, vpr_fv1);
    const auto cp_vpr_fv1_fv2 = cross_product(vpr_fv1, vpr_fv2);  
    const auto cp_vpr_fv2_fv0 = cross_product(vpr_fv2, vpr_fv0);    
        
    // checking the direction with the normal
    const auto n_cp_vpr_fv0_fv1 = n*cp_vpr_fv0_fv1;
    const auto n_cp_vpr_fv1_fv2 = n*cp_vpr_fv1_fv2;
    const auto n_cp_vpr_fv2_fv0 = n*cp_vpr_fv2_fv0;

    // if 'true' the vpr is on the (convex) polygon.
    if ((n_cp_vpr_fv0_fv1 >0 && n_cp_vpr_fv1_fv2 >0 && n_cp_vpr_fv2_fv0 >0) ||
        (n_cp_vpr_fv0_fv1 <0 && n_cp_vpr_fv1_fv2 <0 && n_cp_vpr_fv2_fv0 <0)) {

      const auto closest_distance = std::sqrt(dis_v_poly_sq);
      const auto compression = radius - closest_distance;
      if (compression > 0) {
        compression_vector += -compression * n;
        contact_flag = true;
      }
      continue;
    }
    
   
    
    // direction vector of lines describing the edges
    const auto l01 = fv1 - fv0;
    const auto l12 = fv2 - fv1;
    const auto l20 = fv0 - fv2;    

    // projection of v on the lines containing edges. if 'd' is the direction 
    // vector of a line constructed from 'fv1' and 'fv0', d=(fv1-fv0)/|fv1-fv0|
    // then the projection of the point 'v' would be : vpr = fv0 + ((v-fv0).d)d
    const auto vpr_l01 = fv0 + l01*(((v - fv0)*l01) / (l01*l01)); 
    const auto vpr_l12 = fv1 + l12*(((v - fv1)*l12) / (l12*l12)); 
    const auto vpr_l20 = fv2 + l20*(((v - fv2)*l20) / (l20*l20));     

    // checking whether vpr_lXY is between fvX and fvY, just by checking a dot
    // product of two vectors made by them.
    const bool vpr_on_edge01 = ((vpr_l01-fv0)*(vpr_l01-fv1) <= 0? true : false);
    const bool vpr_on_edge12 = ((vpr_l12-fv1)*(vpr_l12-fv2) <= 0? true : false);
    const bool vpr_on_edge20 = ((vpr_l20-fv2)*(vpr_l20-fv0) <= 0? true : false);
        
    // the same as the formula dis_v_lXY_sq
    const auto dist_v_l01_sq = (v - vpr_l01)*(v - vpr_l01);
    const auto dist_v_l12_sq = (v - vpr_l12)*(v - vpr_l12);
    const auto dist_v_l20_sq = (v - vpr_l20)*(v - vpr_l20);
    
    int closest_edge = -1;
    double dist_closest_edge_sq = -1;
    Vector<Real_t> closest_egde_vector {0,0,0};
    
    if (vpr_on_edge01) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l01_sq) {
        closest_edge = 01;
        dist_closest_edge_sq = dist_v_l01_sq;
        closest_egde_vector = v - vpr_l01;
      }
    }
    if (vpr_on_edge12) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l12_sq) {
        closest_edge = 12;
        dist_closest_edge_sq = dist_v_l12_sq;
        closest_egde_vector = v - vpr_l12;        
      }
    }      
    if (vpr_on_edge20) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l20_sq) {
        closest_edge = 20;
        dist_closest_edge_sq = dist_v_l20_sq;
        closest_egde_vector = v - vpr_l20;        
      }
    }
    
    if (closest_edge != -1) { 
      const auto closest_distance = std::sqrt(dist_closest_edge_sq);
      const auto compression = radius - closest_distance;
      if (compression > 0) {
        compression_vector += -compression * closest_egde_vector/ closest_distance;
        contact_flag = true;
      }
      continue;      
    } 
    
    
    // vectors from 'v' to vertices
    auto v_fv0 = v - fv0;
    auto v_fv1 = v - fv1;
    auto v_fv2 = v - fv2;
            
    // square distance from 'v' to vertices of the polyhedron    
    auto dis_v_fv0_sq = v_fv0 * v_fv0;
    auto dis_v_fv1_sq = v_fv1 * v_fv1;
    auto dis_v_fv2_sq = v_fv2 * v_fv2;
    
    // finding the closest distance from 'v' to polygon
    Real_t dist_closest_vertex_sq = dis_v_fv0_sq;    
    Vector<Real_t> closest_vertex_vector = v_fv0;
    
    if (dist_closest_vertex_sq > dis_v_fv1_sq) {
      dist_closest_vertex_sq = dis_v_fv1_sq;
      closest_vertex_vector = v_fv1;
    }
    if (dist_closest_vertex_sq > dis_v_fv2_sq) {
      dist_closest_vertex_sq = dis_v_fv2_sq;
      closest_vertex_vector = v_fv2;
    }    

    {
      auto closest_distance = std::sqrt(dist_closest_vertex_sq);
      auto compression = radius - closest_distance;
      if (compression > 0) {
        compression_vector += -compression * closest_vertex_vector / closest_distance;
        contact_flag = true;
      }
    }

  }

  contact_vector = compression_vector;
  return contact_flag;
}

// this is called in polyhedron::Handler->in_contact() if activated.
// the same algorithm as (in_contact_grid), except all the polygons are going to be checked
// any fix should be done on both of them.
// for the optimization reasons, we won't make them one functions.
bool Point_Inside::in_contact_all (shape::polyhedron::Polyhedron & p_object,
    const Vector<double> &v, const Real_t radius, Vector<double> &contact_vector) {

  
  const auto & vertex = p_object.vertex;
  const auto & face = p_object.face;
  const auto & normal = p_object.normal;
   
  Vector<Real_t> compression_vector {0.0,0.0,0.0};
     
  bool contact_flag = false;
    
  for (unsigned int i = 0; i< face.size(); ++i)  {

    // renaming the vertices of the 
    const auto &fv0 = vertex[face[i][0]];
    const auto &fv1 = vertex[face[i][1]];    
    const auto &fv2 = vertex[face[i][2]];
    
    // renaming the normal    
    const auto &n = normal[i];

    // distance from v to the plane containing polygon, using fv0,
    // positve value for dis_v_poly means that v is on the side which normal vector
    // points to.
    const auto dis_v_poly = (v - fv0)*n;
    const auto dis_v_poly_sq = dis_v_poly*dis_v_poly;
    
    // projection of v on the polyhedron
    const auto vpr = v - (dis_v_poly*n);

    // vectors from the projection point to vertices
    const auto vpr_fv0 = fv0 - vpr;
    const auto vpr_fv1 = fv1 - vpr;
    const auto vpr_fv2 = fv2 - vpr;

    // cross productions of vpr_fvX,
    const auto cp_vpr_fv0_fv1 = cross_product(vpr_fv0, vpr_fv1);
    const auto cp_vpr_fv1_fv2 = cross_product(vpr_fv1, vpr_fv2);  
    const auto cp_vpr_fv2_fv0 = cross_product(vpr_fv2, vpr_fv0);    
        
    // checking the direction with the normal
    const auto n_cp_vpr_fv0_fv1 = n*cp_vpr_fv0_fv1;
    const auto n_cp_vpr_fv1_fv2 = n*cp_vpr_fv1_fv2;
    const auto n_cp_vpr_fv2_fv0 = n*cp_vpr_fv2_fv0;

    // if 'true' the vpr is on the (convex) polygon.
    if ((n_cp_vpr_fv0_fv1 >0 && n_cp_vpr_fv1_fv2 >0 && n_cp_vpr_fv2_fv0 >0) ||
        (n_cp_vpr_fv0_fv1 <0 && n_cp_vpr_fv1_fv2 <0 && n_cp_vpr_fv2_fv0 <0)) {

      const auto closest_distance = std::sqrt(dis_v_poly_sq);
      const auto compression = radius - closest_distance;
      if (compression > 0) {
        compression_vector += -compression * n;
        contact_flag = true;
      }
      continue;
    }
    
   
    
    // direction vector of lines describing the edges
    const auto l01 = fv1 - fv0;
    const auto l12 = fv2 - fv1;
    const auto l20 = fv0 - fv2;    

    // projection of v on the lines containing edges. if 'd' is the direction 
    // vector of a line constructed from 'fv1' and 'fv0', d=(fv1-fv0)/|fv1-fv0|
    // then the projection of the point 'v' would be : vpr = fv0 + ((v-fv0).d)d
    const auto vpr_l01 = fv0 + l01*(((v - fv0)*l01) / (l01*l01)); 
    const auto vpr_l12 = fv1 + l12*(((v - fv1)*l12) / (l12*l12)); 
    const auto vpr_l20 = fv2 + l20*(((v - fv2)*l20) / (l20*l20));     

    // checking whether vpr_lXY is between fvX and fvY, just by checking a dot
    // product of two vectors made by them.
    const bool vpr_on_edge01 = ((vpr_l01-fv0)*(vpr_l01-fv1) <= 0? true : false);
    const bool vpr_on_edge12 = ((vpr_l12-fv1)*(vpr_l12-fv2) <= 0? true : false);
    const bool vpr_on_edge20 = ((vpr_l20-fv2)*(vpr_l20-fv0) <= 0? true : false);
        
    // the same as the formula dis_v_lXY_sq
    const auto dist_v_l01_sq = (v - vpr_l01)*(v - vpr_l01);
    const auto dist_v_l12_sq = (v - vpr_l12)*(v - vpr_l12);
    const auto dist_v_l20_sq = (v - vpr_l20)*(v - vpr_l20);
    
    int closest_edge = -1;
    double dist_closest_edge_sq = -1;
    Vector<Real_t> closest_egde_vector {0,0,0};
    
    if (vpr_on_edge01) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l01_sq) {
        closest_edge = 01;
        dist_closest_edge_sq = dist_v_l01_sq;
        closest_egde_vector = v - vpr_l01;
      }
    }
    if (vpr_on_edge12) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l12_sq) {
        closest_edge = 12;
        dist_closest_edge_sq = dist_v_l12_sq;
        closest_egde_vector = v - vpr_l12;        
      }
    }      
    if (vpr_on_edge20) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l20_sq) {
        closest_edge = 20;
        dist_closest_edge_sq = dist_v_l20_sq;
        closest_egde_vector = v - vpr_l20;        
      }
    }
    
    if (closest_edge != -1) { 
      const auto closest_distance = std::sqrt(dist_closest_edge_sq);
      const auto compression = radius - closest_distance;
      if (compression > 0) {
        compression_vector += -compression * closest_egde_vector/ closest_distance;
        contact_flag = true;
      }
      continue;      
    } 
    
    
    // vectors from 'v' to vertices
    auto v_fv0 = v - fv0;
    auto v_fv1 = v - fv1;
    auto v_fv2 = v - fv2;
            
    // square distance from 'v' to vertices of the polyhedron    
    auto dis_v_fv0_sq = v_fv0 * v_fv0;
    auto dis_v_fv1_sq = v_fv1 * v_fv1;
    auto dis_v_fv2_sq = v_fv2 * v_fv2;
    
    // finding the closest distance from 'v' to polygon
    Real_t dist_closest_vertex_sq = dis_v_fv0_sq;    
    Vector<Real_t> closest_vertex_vector = v_fv0;
    
    if (dist_closest_vertex_sq > dis_v_fv1_sq) {
      dist_closest_vertex_sq = dis_v_fv1_sq;
      closest_vertex_vector = v_fv1;
    }
    if (dist_closest_vertex_sq > dis_v_fv2_sq) {
      dist_closest_vertex_sq = dis_v_fv2_sq;
      closest_vertex_vector = v_fv2;
    }    

    {
      auto closest_distance = std::sqrt(dist_closest_vertex_sq);
      auto compression = radius - closest_distance;
      if (compression > 0) {
        compression_vector += -compression * closest_vertex_vector / closest_distance;
        contact_flag = true;
      }
    }

  }

  contact_vector = compression_vector;
  return contact_flag;
}



bool Point_Inside::ray_tells_point_is_inside (shape::polyhedron::Polyhedron &p_object, const Vector<Real_t> &v, const int ray_axis) { 


    const auto & vertex = p_object.vertex;
    const auto & face = p_object.face;
    const auto & normal = p_object.normal;

  
    bool inside_shape = false; // inside atleast one of the shapes
  
    for (unsigned int i=0; i<face.size(); ++i) {

      Vector<Real_t> v_pr = v;
      const auto v0 = vertex[face[i][0]];
      const auto n = normal[i];

#define VEC_GREATER_THAN_VERTEX(VAR,AXIS) \
   if (VAR.AXIS > vertex[face[i][0]].AXIS && VAR.AXIS > vertex[face[i][1]].AXIS && VAR.AXIS > vertex[face[i][2]].AXIS ) continue;
#define VEC_LESS_THAN_VERTEX(VAR,AXIS) \
  if (VAR.AXIS < vertex[face[i][0]].AXIS && VAR.AXIS < vertex[face[i][1]].AXIS && VAR.AXIS < vertex[face[i][2]].AXIS ) continue;
      
      VEC_GREATER_THAN_VERTEX(v,z)
      VEC_GREATER_THAN_VERTEX(v,x)
      VEC_GREATER_THAN_VERTEX(v,y)      
      
      if (ray_axis == 0) {
         if (std::abs(n.x) < 1e-9) continue;        
        v_pr.x = v0.x + (n.z*(-v.z + v0.z) + n.y*(-v.y + v0.y))/n.x;
        VEC_LESS_THAN_VERTEX(v,y)
        VEC_LESS_THAN_VERTEX(v,z)
        VEC_LESS_THAN_VERTEX(v_pr,x)

        VEC_GREATER_THAN_VERTEX(v_pr,x)        

      } else if (ray_axis == 1) {
         if (std::abs(n.y) < 1e-9) continue;        
        v_pr.y = v0.y + (n.x*(-v.x + v0.x) + n.z*(-v.z + v0.z))/n.y;
        VEC_LESS_THAN_VERTEX(v,x)
        VEC_LESS_THAN_VERTEX(v,z)
        VEC_LESS_THAN_VERTEX(v_pr,y)    
        VEC_GREATER_THAN_VERTEX(v_pr,y)        

      } else if (ray_axis == 2) {
         if (std::abs(n.z) < 1e-9) continue;        
        v_pr.z = v0.z + (n.x*(-v.x + v0.x) + n.y*(-v.y + v0.y))/n.z;
        VEC_LESS_THAN_VERTEX(v,x)
        VEC_LESS_THAN_VERTEX(v,y)
        VEC_LESS_THAN_VERTEX(v_pr,z)     
        VEC_GREATER_THAN_VERTEX(v_pr,z)        

      } else continue;
      
#undef VEC_LESS_THAN_VERTEX
#undef VEC_GREATER_THAN_VERTEX

      const auto PA = vertex[face[i][0]] - v_pr;  
      const auto PB = vertex[face[i][1]] - v_pr;        
      const auto PC = vertex[face[i][2]] - v_pr;
      
      const auto PAB = cross_product(PA, PB);
      const auto PBC = cross_product(PB, PC);      
      const auto PCA = cross_product(PC, PA);
      
      const auto PABPBC = PAB*PBC;
      const auto PBCPCA = PBC*PCA;
      const auto PCAPAB = PCA*PAB;
      
      if ((PABPBC<0 && PBCPCA<0 && PCAPAB<0) || (PABPBC>0 && PBCPCA>0 && PCAPAB>0)){
        inside_shape = !inside_shape;        
      } 
      
    }
    
  // inside at least one of the shapes    
  if (inside_shape) return true; 
   
  return false;
    
}


// usage in particle initial distributions (random, grid or ...)
// the same algorithm as (in_contact_all), except all the polygons are going to be checked
// any fix should be done on both of them.
// for the optimization reasons, we won't make them one functions.
bool Point_Inside::in_contact_grid (shape::polyhedron::Polyhedron & p_object, 
    const Vector<double> &v, const Real_t radius) {

  
  const auto & vertex = p_object.vertex;
  const auto & face = p_object.face;
  const auto & normal = p_object.normal;
   
    
  const auto & grid = p_object.grid;
   
  const auto & xlo = p_object.xlo;
  const auto & ylo = p_object.ylo;
  const auto & zlo = p_object.zlo;
  
  const auto & dx_part = p_object.dx_part;
  const auto & dy_part = p_object.dy_part;
  const auto & dz_part = p_object.dz_part;
  
  const auto & nx_part = p_object.nx_part;
  const auto & ny_part = p_object.ny_part;
  const auto & nz_part = p_object.nz_part;

  int xindex = int_floor((v.x-xlo)/dx_part);
  int yindex = int_floor((v.y-ylo)/dy_part);
  int zindex = int_floor((v.z-zlo)/dz_part);

  if (xindex<0) return false;
  if (yindex<0) return false;
  if (zindex<0) return false;

  if (xindex>int_floor(nx_part-1)) return false;
  if (yindex>int_floor(ny_part-1)) return false;
  if (zindex>int_floor(nz_part-1)) return false;


    
  for (auto i : grid[xindex][yindex][zindex]) {

    // renaming the vertices of the 
    const auto &fv0 = vertex[face[i][0]];
    const auto &fv1 = vertex[face[i][1]];    
    const auto &fv2 = vertex[face[i][2]];
    
    // renaming the normal    
    const auto &n = normal[i];

    // distance from v to the plane containing polygon, using fv0,
    // positve value for dis_v_poly means that v is on the side which normal vector
    // points to.
    const auto dis_v_poly = (v - fv0)*n;
    const auto dis_v_poly_sq = dis_v_poly*dis_v_poly;
    
    // projection of v on the polyhedron
    const auto vpr = v - (dis_v_poly*n);

    // vectors from the projection point to vertices
    const auto vpr_fv0 = fv0 - vpr;
    const auto vpr_fv1 = fv1 - vpr;
    const auto vpr_fv2 = fv2 - vpr;

    // cross productions of vpr_fvX,
    const auto cp_vpr_fv0_fv1 = cross_product(vpr_fv0, vpr_fv1);
    const auto cp_vpr_fv1_fv2 = cross_product(vpr_fv1, vpr_fv2);  
    const auto cp_vpr_fv2_fv0 = cross_product(vpr_fv2, vpr_fv0);    
        
    // checking the direction with the normal
    const auto n_cp_vpr_fv0_fv1 = n*cp_vpr_fv0_fv1;
    const auto n_cp_vpr_fv1_fv2 = n*cp_vpr_fv1_fv2;
    const auto n_cp_vpr_fv2_fv0 = n*cp_vpr_fv2_fv0;

    // if 'true' the vpr is on the (convex) polygon.
    if ((n_cp_vpr_fv0_fv1 >0 && n_cp_vpr_fv1_fv2 >0 && n_cp_vpr_fv2_fv0 >0) ||
        (n_cp_vpr_fv0_fv1 <0 && n_cp_vpr_fv1_fv2 <0 && n_cp_vpr_fv2_fv0 <0)) {

      const auto closest_distance = std::sqrt(dis_v_poly_sq);
      const auto compression = radius - closest_distance;
      if (compression > 0) {
        return true;
      }
      continue;
    }
    
   
    
    // direction vector of lines describing the edges
    const auto l01 = fv1 - fv0;
    const auto l12 = fv2 - fv1;
    const auto l20 = fv0 - fv2;    

    // projection of v on the lines containing edges. if 'd' is the direction 
    // vector of a line constructed from 'fv1' and 'fv0', d=(fv1-fv0)/|fv1-fv0|
    // then the projection of the point 'v' would be : vpr = fv0 + ((v-fv0).d)d
    const auto vpr_l01 = fv0 + l01*(((v - fv0)*l01) / (l01*l01)); 
    const auto vpr_l12 = fv1 + l12*(((v - fv1)*l12) / (l12*l12)); 
    const auto vpr_l20 = fv2 + l20*(((v - fv2)*l20) / (l20*l20));     

    // checking whether vpr_lXY is between fvX and fvY, just by checking a dot
    // product of two vectors made by them.
    const bool vpr_on_edge01 = ((vpr_l01-fv0)*(vpr_l01-fv1) <= 0? true : false);
    const bool vpr_on_edge12 = ((vpr_l12-fv1)*(vpr_l12-fv2) <= 0? true : false);
    const bool vpr_on_edge20 = ((vpr_l20-fv2)*(vpr_l20-fv0) <= 0? true : false);
        
    // the same as the formula dis_v_lXY_sq
    const auto dist_v_l01_sq = (v - vpr_l01)*(v - vpr_l01);
    const auto dist_v_l12_sq = (v - vpr_l12)*(v - vpr_l12);
    const auto dist_v_l20_sq = (v - vpr_l20)*(v - vpr_l20);
    
    int closest_edge = -1;
    double dist_closest_edge_sq = -1;
    Vector<Real_t> closest_egde_vector {0,0,0};
    
    if (vpr_on_edge01) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l01_sq) {
        closest_edge = 01;
        dist_closest_edge_sq = dist_v_l01_sq;
        closest_egde_vector = v - vpr_l01;
      }
    }
    if (vpr_on_edge12) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l12_sq) {
        closest_edge = 12;
        dist_closest_edge_sq = dist_v_l12_sq;
        closest_egde_vector = v - vpr_l12;        
      }
    }      
    if (vpr_on_edge20) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l20_sq) {
        closest_edge = 20;
        dist_closest_edge_sq = dist_v_l20_sq;
        closest_egde_vector = v - vpr_l20;        
      }
    }
    
    if (closest_edge != -1) { 
      const auto closest_distance = std::sqrt(dist_closest_edge_sq);
      const auto compression = radius - closest_distance;
      if (compression > 0) {
        return true;
      }
      continue;      
    } 
    
    
    // vectors from 'v' to vertices
    auto v_fv0 = v - fv0;
    auto v_fv1 = v - fv1;
    auto v_fv2 = v - fv2;
            
    // square distance from 'v' to vertices of the polyhedron    
    auto dis_v_fv0_sq = v_fv0 * v_fv0;
    auto dis_v_fv1_sq = v_fv1 * v_fv1;
    auto dis_v_fv2_sq = v_fv2 * v_fv2;
    
    // finding the closest distance from 'v' to polygon
    Real_t dist_closest_vertex_sq = dis_v_fv0_sq;    
    Vector<Real_t> closest_vertex_vector = v_fv0;
    
    if (dist_closest_vertex_sq > dis_v_fv1_sq) {
      dist_closest_vertex_sq = dis_v_fv1_sq;
      closest_vertex_vector = v_fv1;
    }
    if (dist_closest_vertex_sq > dis_v_fv2_sq) {
      dist_closest_vertex_sq = dis_v_fv2_sq;
      closest_vertex_vector = v_fv2;
    }    

    {
      auto closest_distance = std::sqrt(dist_closest_vertex_sq);
      auto compression = radius - closest_distance;
      if (compression > 0) {
        return true;
      }
    }

  }

  return false;
}


  // usage in particle initial distributions (random, grid or ...)
// the same algorithm as (in_contact_grid), except all the polygons are going to be checked
// any fix should be done on both of them.
// for the optimization reasons, we won't make them one functions.
bool Point_Inside::in_contact_all (shape::polyhedron::Polyhedron & p_object,
    const Vector<double> &v, const Real_t radius) {

  
  const auto & vertex = p_object.vertex;
  const auto & face = p_object.face;
  const auto & normal = p_object.normal;
   
     
   
  for (unsigned int i = 0; i< face.size(); ++i)  {

    // renaming the vertices of the 
    const auto &fv0 = vertex[face[i][0]];
    const auto &fv1 = vertex[face[i][1]];    
    const auto &fv2 = vertex[face[i][2]];
    
    // renaming the normal    
    const auto &n = normal[i];

    // distance from v to the plane containing polygon, using fv0,
    // positve value for dis_v_poly means that v is on the side which normal vector
    // points to.
    const auto dis_v_poly = (v - fv0)*n;
    const auto dis_v_poly_sq = dis_v_poly*dis_v_poly;
    
    // projection of v on the polyhedron
    const auto vpr = v - (dis_v_poly*n);

    // vectors from the projection point to vertices
    const auto vpr_fv0 = fv0 - vpr;
    const auto vpr_fv1 = fv1 - vpr;
    const auto vpr_fv2 = fv2 - vpr;

    // cross productions of vpr_fvX,
    const auto cp_vpr_fv0_fv1 = cross_product(vpr_fv0, vpr_fv1);
    const auto cp_vpr_fv1_fv2 = cross_product(vpr_fv1, vpr_fv2);  
    const auto cp_vpr_fv2_fv0 = cross_product(vpr_fv2, vpr_fv0);    
        
    // checking the direction with the normal
    const auto n_cp_vpr_fv0_fv1 = n*cp_vpr_fv0_fv1;
    const auto n_cp_vpr_fv1_fv2 = n*cp_vpr_fv1_fv2;
    const auto n_cp_vpr_fv2_fv0 = n*cp_vpr_fv2_fv0;

    // if 'true' the vpr is on the (convex) polygon.
    if ((n_cp_vpr_fv0_fv1 >0 && n_cp_vpr_fv1_fv2 >0 && n_cp_vpr_fv2_fv0 >0) ||
        (n_cp_vpr_fv0_fv1 <0 && n_cp_vpr_fv1_fv2 <0 && n_cp_vpr_fv2_fv0 <0)) {

      const auto closest_distance = std::sqrt(dis_v_poly_sq);
      const auto compression = radius - closest_distance;
      if (compression > 0) {
        return true;
      }
      continue;
    }
    
   
    
    // direction vector of lines describing the edges
    const auto l01 = fv1 - fv0;
    const auto l12 = fv2 - fv1;
    const auto l20 = fv0 - fv2;    

    // projection of v on the lines containing edges. if 'd' is the direction 
    // vector of a line constructed from 'fv1' and 'fv0', d=(fv1-fv0)/|fv1-fv0|
    // then the projection of the point 'v' would be : vpr = fv0 + ((v-fv0).d)d
    const auto vpr_l01 = fv0 + l01*(((v - fv0)*l01) / (l01*l01)); 
    const auto vpr_l12 = fv1 + l12*(((v - fv1)*l12) / (l12*l12)); 
    const auto vpr_l20 = fv2 + l20*(((v - fv2)*l20) / (l20*l20));     

    // checking whether vpr_lXY is between fvX and fvY, just by checking a dot
    // product of two vectors made by them.
    const bool vpr_on_edge01 = ((vpr_l01-fv0)*(vpr_l01-fv1) <= 0? true : false);
    const bool vpr_on_edge12 = ((vpr_l12-fv1)*(vpr_l12-fv2) <= 0? true : false);
    const bool vpr_on_edge20 = ((vpr_l20-fv2)*(vpr_l20-fv0) <= 0? true : false);
        
    // the same as the formula dis_v_lXY_sq
    const auto dist_v_l01_sq = (v - vpr_l01)*(v - vpr_l01);
    const auto dist_v_l12_sq = (v - vpr_l12)*(v - vpr_l12);
    const auto dist_v_l20_sq = (v - vpr_l20)*(v - vpr_l20);
    
    int closest_edge = -1;
    double dist_closest_edge_sq = -1;
    Vector<Real_t> closest_egde_vector {0,0,0};
    
    if (vpr_on_edge01) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l01_sq) {
        closest_edge = 01;
        dist_closest_edge_sq = dist_v_l01_sq;
        closest_egde_vector = v - vpr_l01;
      }
    }
    if (vpr_on_edge12) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l12_sq) {
        closest_edge = 12;
        dist_closest_edge_sq = dist_v_l12_sq;
        closest_egde_vector = v - vpr_l12;        
      }
    }      
    if (vpr_on_edge20) {
      if (closest_edge == -1 || dist_closest_edge_sq > dist_v_l20_sq) {
        closest_edge = 20;
        dist_closest_edge_sq = dist_v_l20_sq;
        closest_egde_vector = v - vpr_l20;        
      }
    }
    
    if (closest_edge != -1) { 
      const auto closest_distance = std::sqrt(dist_closest_edge_sq);
      const auto compression = radius - closest_distance;
      if (compression > 0) {
        return true;
      }
      continue;      
    } 
    
    
    // vectors from 'v' to vertices
    auto v_fv0 = v - fv0;
    auto v_fv1 = v - fv1;
    auto v_fv2 = v - fv2;
            
    // square distance from 'v' to vertices of the polyhedron    
    auto dis_v_fv0_sq = v_fv0 * v_fv0;
    auto dis_v_fv1_sq = v_fv1 * v_fv1;
    auto dis_v_fv2_sq = v_fv2 * v_fv2;
    
    // finding the closest distance from 'v' to polygon
    Real_t dist_closest_vertex_sq = dis_v_fv0_sq;    
    Vector<Real_t> closest_vertex_vector = v_fv0;
    
    if (dist_closest_vertex_sq > dis_v_fv1_sq) {
      dist_closest_vertex_sq = dis_v_fv1_sq;
      closest_vertex_vector = v_fv1;
    }
    if (dist_closest_vertex_sq > dis_v_fv2_sq) {
      dist_closest_vertex_sq = dis_v_fv2_sq;
      closest_vertex_vector = v_fv2;
    }    

    {
      auto closest_distance = std::sqrt(dist_closest_vertex_sq);
      auto compression = radius - closest_distance;
      if (compression > 0) {
        return true;
      }
    }

  }

  return false;
}

} //polyhedron
} //shape

} // namespace caviar

