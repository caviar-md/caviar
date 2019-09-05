
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
	location <-10,-10,-10>
	look_at  <0,0,0>
}

light_source { <-12, -4, -3> color White}
light_source { <-2, -40, -30> color White}
#include "colors.inc"

camera {
	location <-10,-10,-10>
	look_at  <0,0,0>
}

light_source { <-12, -4, -3> color White}
light_source { <-2, -40, -30> color White}

union {
	sphere {<40.382,-22.6001,20> 1 }
	sphere {<253.123,25.8939,20> 2 }
	sphere {<167.487,15.7405,20> 3}
	sphere {<88.6407,5.4224,20> 4}
	sphere {<45.4901,-0.895497,20>  5}
	texture {
		pigment { color Yellow }
 }
}
