
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

#include "caviar/objects/shape/polyhedron/postprocess.h"
#include "caviar/objects/shape/polyhedron/polyhedron.h"

#include <string>
#include <cmath>
#include <fstream>

CAVIAR_NAMESPACE_OPEN

namespace shape
{
  namespace polyhedron
  {

    static inline int int_floor(double x)
    {
      return (int)(x + 100000) - 100000;
    }

    Postprocess::Postprocess(CAVIAR *fptr) : Pointers{fptr} {}

    Postprocess::~Postprocess() {}

    void Postprocess::lowest_highest_coord(shape::polyhedron::Polyhedron &p_object)
    {

      const auto &vertex = p_object.vertex;
      const auto &tol = p_object.grid_tol;
      auto &xlo = p_object.xlo;
      auto &xhi = p_object.xhi;
      auto &ylo = p_object.ylo;
      auto &yhi = p_object.yhi;
      auto &zlo = p_object.zlo;
      auto &zhi = p_object.zhi;

      /*
        if (!g_coordinates_init) {
          g_coordinates_init = true;
          gxlo = vertex[0].x;  gylo = vertex[0].y;
          gzlo = vertex[0].z;  gxhi = vertex[0].x;
          gyhi = vertex[0].y;  gzhi = vertex[0].z;
        }
        */
      xlo = vertex[0].x;
      ylo = vertex[0].y;
      zlo = vertex[0].z;
      xhi = vertex[0].x;
      yhi = vertex[0].y;
      zhi = vertex[0].z;

      for (auto v : vertex)
      {
        if (xlo > v.x)
          xlo = v.x;
        if (ylo > v.y)
          ylo = v.y;
        if (zlo > v.z)
          zlo = v.z;
        if (xhi < v.x)
          xhi = v.x;
        if (yhi < v.y)
          yhi = v.y;
        if (zhi < v.z)
          zhi = v.z;
      }

      xlo -= tol;
      xhi += tol;
      ylo -= tol;
      yhi += tol;
      zlo -= tol;
      zhi += tol;

      {
        //    if (gxlo > xlo) gxlo = xlo;    if (gylo > ylo) gylo = ylo;    if (gzlo > zhi) gzlo = zlo;
        //    if (gxhi < xhi) gxhi = xhi;    if (gyhi < yhi) gyhi = yhi;    if (gzhi < zlo) gzhi = zhi;
      }

      std::cout << " min. and max. of Preprocess shape "
                << " ? "
                << " coordinates:"
                << " xlo: " << xlo << " xhi: " << xhi
                << " ylo: " << ylo << " yhi: " << yhi
                << " zlo: " << zlo << " zhi: " << zhi << std::endl;
    }

    void Postprocess::make_grid(shape::polyhedron::Polyhedron &p_object)
    {
      const auto &vertex = p_object.vertex;
      const auto &face = p_object.face;

      auto &grid = p_object.grid;

      const auto &xlo = p_object.xlo, &xhi = p_object.xhi;
      const auto &ylo = p_object.ylo, &yhi = p_object.yhi;
      const auto &zlo = p_object.zlo, &zhi = p_object.zhi;
      auto &dx_part = p_object.dx_part;
      auto &dy_part = p_object.dy_part;
      auto &dz_part = p_object.dz_part;

      const auto &nx_part = p_object.nx_part;
      const auto &ny_part = p_object.ny_part;
      const auto &nz_part = p_object.nz_part;

      // tol==0 makes a bug. I think tol has to be larger than biggest particle
      // radius. or more maybe because of the pointy edges.
      const auto &tol = p_object.grid_tol;

      grid.resize(nx_part);
      for (unsigned int i = 0; i < nx_part; ++i)
        grid[i].resize(ny_part);

      for (unsigned int i = 0; i < nx_part; ++i)
        for (unsigned int j = 0; j < ny_part; ++j)
          grid[i][j].resize(nz_part);

      dx_part = (xhi - xlo) / nx_part;
      dy_part = (yhi - ylo) / ny_part;
      dz_part = (zhi - zlo) / nz_part;

      for (unsigned int i = 0; i < face.size(); ++i)
      {
        Real_t fxlo = vertex[face[i][0]].x,
               fylo = vertex[face[i][0]].y,
               fzlo = vertex[face[i][0]].z;

        Real_t fxhi = fxlo,
               fyhi = fylo,
               fzhi = fzlo;

        for (unsigned int j = 1; j < face[i].size(); ++j)
        {
          auto vx = vertex[face[i][j]].x;
          auto vy = vertex[face[i][j]].y;
          auto vz = vertex[face[i][j]].z;

          if (fxlo > vx)
            fxlo = vx;
          if (fxhi < vx)
            fxhi = vx;
          if (fylo > vy)
            fylo = vy;
          if (fyhi < vy)
            fyhi = vy;
          if (fzlo > vz)
            fzlo = vz;
          if (fzhi < vz)
            fzhi = vz;
        }

        fxlo -= tol;
        fylo -= tol;
        fzlo -= tol;
        fxhi += tol;
        fyhi += tol;
        fzhi += tol;

        int xindex_lo = int_floor((fxlo - xlo) / dx_part);
        int xindex_hi = int_floor((fxhi - xlo) / dx_part);
        int yindex_lo = int_floor((fylo - ylo) / dy_part);
        int yindex_hi = int_floor((fyhi - ylo) / dy_part);
        int zindex_lo = int_floor((fzlo - zlo) / dz_part);
        int zindex_hi = int_floor((fzhi - zlo) / dz_part);

        if (xindex_lo < 0)
          xindex_lo = 0;
        if (yindex_lo < 0)
          yindex_lo = 0;
        if (zindex_lo < 0)
          zindex_lo = 0;

        if (xindex_hi > int_floor(nx_part - 1))
          xindex_hi = nx_part - 1;
        if (yindex_hi > int_floor(ny_part - 1))
          yindex_hi = ny_part - 1;
        if (zindex_hi > int_floor(nz_part - 1))
          zindex_hi = nz_part - 1;

        for (int jx = xindex_lo; jx <= xindex_hi; ++jx)
          for (int jy = yindex_lo; jy <= yindex_hi; ++jy)
            for (int jz = zindex_lo; jz <= zindex_hi; ++jz)
            {

              grid[jx][jy][jz].push_back(i);
            }
      }

      /*
        for (unsigned int i=0;i<face.size();++i) {
          Real_t fxlo = vertex[face[i][0]].x - toll,
                 fxhi = vertex[face[i][0]].x + toll,
                 fylo = vertex[face[i][0]].y - toll,
                 fyhi = vertex[face[i][0]].y + toll,
                 fzlo = vertex[face[i][0]].z - toll,
                 fzhi = vertex[face[i][0]].z + toll;

          for (unsigned int j=1;j<face[i].size();++j) {
            if (fxlo > vertex[face[i][j]].x - toll)
              fxlo = vertex[face[i][j]].x - toll;
            if (fxhi < vertex[face[i][j]].x + toll)
              fxhi = vertex[face[i][j]].x + toll;
            if (fylo > vertex[face[i][j]].y - toll)
              fylo = vertex[face[i][j]].y - toll;
            if (fyhi < vertex[face[i][j]].y + toll)
              fyhi = vertex[face[i][j]].y + toll;
            if (fzlo > vertex[face[i][j]].z - toll)
              fzlo = vertex[face[i][j]].z - toll;
            if (fzhi < vertex[face[i][j]].z + toll)
              fzhi = vertex[face[i][j]].z + toll;
          }

          int xindex_lo = int((fxlo-xlo)/dx_part);
          int xindex_hi = int((fxhi-xlo)/dx_part);
          int yindex_lo = int((fylo-ylo)/dy_part);
          int yindex_hi = int((fyhi-ylo)/dy_part);
          int zindex_lo = int((fzlo-zlo)/dz_part);
          int zindex_hi = int((fzhi-zlo)/dz_part);


          if (xindex_lo<0) xindex_lo = 0;
          if (yindex_lo<0) yindex_lo = 0;
          if (zindex_lo<0) zindex_lo = 0;

          if (xindex_hi>int(nx_part-1)) xindex_hi = nx_part-1;
          if (yindex_hi>int(ny_part-1)) yindex_hi = ny_part-1;
          if (zindex_hi>int(nz_part-1)) zindex_hi = nz_part-1;


          for (int jx=xindex_lo; jx<=xindex_hi; ++jx)
            for (int jy=yindex_lo; jy<=yindex_hi; ++jy)
              for (int jz=zindex_lo; jz<=zindex_hi; ++jz) {

                grid[jx][jy][jz].push_back (i);
              }
        }
      */
    }
  } // polyhedron
} // shape

CAVIAR_NAMESPACE_CLOSE
