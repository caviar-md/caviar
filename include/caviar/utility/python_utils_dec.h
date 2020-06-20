
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

#ifndef CAVIAR_PYTHON_UTILS_DEC_H
#define CAVIAR_PYTHON_UTILS_DEC_H

#ifndef BOOST_BIND_GLOBAL_PLACEHOLDERS
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#endif

namespace boost {
namespace python {
class list;
}
}

#define FC_PYDEC_SET_PTR(x,y) \
  void set_ ## x (std::shared_ptr<y> x_ptr);

#define FC_PYDEC_GET_PTR(x,y) \
  std::shared_ptr<y> get_ ## x () ;

#define FC_PYDEC_SETGET_PTR(x,y) \
  FC_PYDEC_SET_PTR(x,y) \
  FC_PYDEC_GET_PTR(x,y)


// T is the type of the std::vector, e.g double
#define FC_PYDEC_GET_STDVEC(x,T) \
  boost::python::list get_ ## x () ;


// T is the type of the std::vector, e.g double
#define FC_PYDEC_SET_STDVEC(x,T) \
  void set_ ## x (boost::python::list& ns) ;


#define FC_PYDEC_GET_STDVEC2D(x,T) \
  boost::python::list get_ ## x() ;


#define FC_PYDEC_SET_STDVEC2D(x,T) \
  void set_ ## x (const std::vector< std::vector< T > > &v) ;

// T is the type of the caviar::Vector, e.g double
#define FC_PYDEC_GET_CAVVEC(x,T) \
  boost::python::list get_ ## x () ;


// T is the type of the caviar::Vector, e.g double
#define FC_PYDEC_SET_CAVVEC(x,T) \
  void set_ ## x (boost::python::list& ns) ;

#define FC_PYDEC_SETGET_CAVVEC(x,y) \
  FC_PYDEC_SET_CAVVEC(x,y) \
  FC_PYDEC_GET_CAVVEC(x,y)


#define FC_PYDEC_SETGET_STDVEC(x,y) \
  FC_PYDEC_SET_STDVEC(x,y) \
  FC_PYDEC_GET_STDVEC(x,y)

//
#define FC_PYDEC_SETGET_STDVEC2D(x,y) \
  FC_PYDEC_SET_STDVEC2D(x,y) \
  FC_PYDEC_GET_STDVEC2D(x,y)



#endif
