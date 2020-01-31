
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

#include "colors.inc"

camera {
	//right x*ImageWidth/ImageHeight
	location <0,100,100>
	look_at  <50,50,50>
}

light_source { <-12, -4, -3> color White}
light_source { <-2, -40, -30> color White}
