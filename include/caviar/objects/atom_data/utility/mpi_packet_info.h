
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

#ifndef CAVIAR_OBJECTS_ATOMDATA_UTILITY_MPI_PACKET_INFO_H
#define CAVIAR_OBJECTS_ATOMDATA_UTILITY_MPI_PACKET_INFO_H

#include "caviar/utility/objects_common_headers.h"

CAVIAR_NAMESPACE_OPEN

namespace atom_data
{
  // Angle contain data for rigid atomic angles which may be used in
  // constraint algorithms or soft atomic angles in spring_angle force_fields
  struct MPI_packet_info
  {
    int id_s, id_e; 
    int type_s,type_e; 
    int pos_s, pos_e; 
    int vel_s, vel_e; 
    int acc_s, acc_e; 
    int msd_s, msd_e; 
    int pos_o_s, pos_o_e; 
    int vel_o_s, vel_o_e;
    int acc_o_s, acc_o_e;        
  } ;
}
CAVIAR_NAMESPACE_CLOSE
#endif