/*
  Teem: Tools to process and visualize scientific data and images              
  Copyright (C) 2011, 2010, 2009  University of Chicago
  Copyright (C) 2008, 2007, 2006, 2005  Gordon Kindlmann
  Copyright (C) 2004, 2003, 2002, 2001, 2000, 1999, 1998  University of Utah

  This library is free software; you can redistribute it and/or
  modify it under the terms of the GNU Lesser General Public License
  (LGPL) as published by the Free Software Foundation; either
  version 2.1 of the License, or (at your option) any later version.
  The terms of redistributing and/or modifying this software also
  include exceptions to the LGPL that facilitate static linking.

  This library is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  Lesser General Public License for more details.

  You should have received a copy of the GNU Lesser General Public License
  along with this library; if not, write to Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
*/

#include "air.h"
#include <teemEndian.h>

/*
******** airMyEndian()
**
** determine at run-time if we are little (1234) or big (4321) endian
*/
int
airMyEndian(void) {
  int tmpI, ret;
  char endian;
  
  tmpI = 1;
  endian = !(*((char*)(&tmpI)));
  if (endian) {
    ret = airEndianBig;
  }
  else {
    ret = airEndianLittle;
  }    
  return ret;
}

const char *
_airEndianStr[] = {
  "(unknown endian)",
  "little",
  "big"
};

const char *
_airEndianDesc[] = {
  "unknown endianness",
  "Intel and compatible",
  "Everyone besides Intel and compatible"
};

const int
_airEndianVal[] = {
  airEndianUnknown,
  airEndianLittle,
  airEndianBig,
};

const airEnum
_airEndian = {
  "endian",
  2,
  _airEndianStr, _airEndianVal,
  _airEndianDesc,
  NULL, NULL,
  AIR_FALSE
};

const airEnum *const
airEndian = &_airEndian;

