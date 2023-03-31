
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

#ifndef CAVIAR_OBJECTS_SHAPE_MIXED_H
#define CAVIAR_OBJECTS_SHAPE_MIXED_H

#include "caviar/objects/shape.h"

namespace caviar {

namespace shape {

/**
 * This class creates mixes of different shapes in order to define new shapes.
 * 
 */
class Mixed : public caviar::Shape {
  public:
    Mixed (class CAVIAR *);
    ~Mixed ();
//    bool read (caviar::interpreter::Parser *, class Object_container *);    
    bool read(class caviar::interpreter::Parser *);
    void satisfy_Mixed ();
    
    bool inside_check;
    
    bool is_inside (const Vector<double> &);
    bool is_inside (const Vector<double> &, const double r);
    bool in_contact(const Vector<double> &v, const double r, Vector<double> & contact_vector);    
    
    //bool is_all (const Vector<double> &); //checks 'is_inside()' if 'inside_check==true'

    std::vector<caviar::Shape*> shapes;
    std::vector<int> operators; // 1:and_inside, -1:and_outside, 2:or_inside, -2:or_outside
    bool shape_add;

};
} //shape

} // namespace caviar
#endif
