
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
//
//

#ifndef CAVIAR_PYTHON_UTILS_DEF_H
#define CAVIAR_PYTHON_UTILS_DEF_H


#include <boost/python.hpp>
#include <boost/python/list.hpp>

/* // raw pointers setters and getters. Not sure if it works.
#define FC_PYDEF_SET_PTR(CLS,VAR,TYP) \
  void CLS :: set_ ## VAR (std::shared_ptr<TYP> x_ptr) { \
    VAR = x_ptr.get(); \
  }

#define FC_PYDEF_GET_PTR(CLS,VAR,TYP) \
  std::shared_ptr<TYP> CLS :: get_ ## VAR () { \
    return std::shared_ptr<TYP> (VAR); \
  }
*/

// Shared pointer
#define FC_PYDEF_SET_PTR(CLS,VAR,TYP) \
  void CLS :: set_ ## VAR (std::shared_ptr<TYP> x_ptr) { \
    VAR = x_ptr; \
  }

#define FC_PYDEF_GET_PTR(CLS,VAR,TYP) \
  std::shared_ptr<TYP> CLS :: get_ ## VAR () { \
    return VAR; \
  }

#define FC_PYDEF_SETGET_PTR(CLS,x,y) \
  FC_PYDEF_SET_PTR(CLS,x,y) \
  FC_PYDEF_GET_PTR(CLS,x,y)


// T is the type of the std::vector, e.g double
#define FC_PYDEF_GET_STDVEC(CLS,VAR,TYP) \
  boost::python::list CLS :: get_ ## VAR () {   \
    boost::python::list list; \
    for (unsigned int i = 0; i < VAR.size(); ++i) { \
      list.append(VAR[i]); \
    } \
    return list; \
  }


// T is the type of the std::vector, e.g double
#define FC_PYDEF_SET_STDVEC(CLS,VAR,TYP) \
  void CLS :: set_ ## VAR (boost::python::list& ns) { \
    VAR.resize(len(ns)); \
    for (unsigned int i = 0; i < len(ns); ++i) { \
      VAR[i] = boost::python::extract<TYP>(ns[i]); \
    } \
  }



#define FC_PYDEF_GET_STDVEC2D(CLS,VAR,TYP) \
  boost::python::list CLS :: get_ ## VAR() {   \
    boost::python::list list; \
    for (unsigned int i = 0; i < VAR.size(); ++i) { \
      boost::python::list listtmp; \
      for (unsigned int j = 0; j < VAR[i].size(); ++j) { \
        listtmp.append(VAR[i][j]); \
      } \
      list.append(listtmp); \
    } \
    return list; \
  }


#define FC_PYDEF_SET_STDVEC2D(CLS,VAR,TYP) \
  void CLS :: set_ ## VAR (const std::vector< std::vector< TYP > > &v) {   \
    VAR.resize(v.size()); \
    for (unsigned int i = 0; i < v.size(); ++i) { \
      VAR[i].resize(v[i].size()); \
      for (unsigned int j = 0; j < v[i].size(); ++j) { \
        VAR[i][j] = v[i][j]; \
      } \
    } \
  }

// T is the type of the caviar::Vector, e.g double
#define FC_PYDEF_GET_CAVVEC(CLS,VAR,TYP) \
  boost::python::list CLS :: get_ ## VAR () {   \
    boost::python::list list; \
    list.append(VAR.x); \
    list.append(VAR.y); \
    list.append(VAR.z); \
    return list; \
  }


// T is the type of the caviar::Vector, e.g double
#define FC_PYDEF_SET_CAVVEC(CLS,VAR,TYP) \
  void CLS :: set_ ## VAR (boost::python::list& ns) { \
    VAR.x = boost::python::extract<TYP>(ns[0]); \
    VAR.y = boost::python::extract<TYP>(ns[1]); \
    VAR.z = boost::python::extract<TYP>(ns[2]); \
  }


#define FC_PYDEF_SETGET_CAVVEC(CLS,VAR,TYP) \
  FC_PYDEF_SET_CAVVEC(CLS,VAR,TYP) \
  FC_PYDEF_GET_CAVVEC(CLS,VAR,TYP)


#define FC_PYDEF_SETGET_STDVEC(CLS,VAR,TYP) \
  FC_PYDEF_SET_STDVEC(CLS,VAR,TYP) \
  FC_PYDEF_GET_STDVEC(CLS,VAR,TYP)

//
#define FC_PYDEF_SETGET_STDVEC2D(CLS,VAR,TYP) \
  FC_PYDEF_SET_STDVEC2D(CLS,VAR,TYP) \
  FC_PYDEF_GET_STDVEC2D(CLS,VAR,TYP)



#endif
