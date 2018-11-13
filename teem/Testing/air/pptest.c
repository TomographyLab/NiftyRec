/*
  Teem: Tools to process and visualize scientific data and images              
  Copyright (c) 2011, 2010, 2009  University of Chicago
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


#include "teem/air.h"

/*
** the purpose of this test is currently limited to generating
** warnings if the format specifications for size_t and ptrdiff_t are
** incorrect.
**
** other tests relating to printing & parsing stirngs will go here later
*/

int
main(int argc, const char *argv[]) {
  size_t sz;
  ptrdiff_t pd;
  char buff[AIR_STRLEN_SMALL],
    stmp1[AIR_STRLEN_SMALL], stmp2[AIR_STRLEN_SMALL];

  AIR_UNUSED(argc);
  
  sz = 424242;
  pd = -424242;
  sprintf(buff, "sz: %s, pd: %s", 
          airSprintSize_t(stmp1, sz), airSprintPtrdiff_t(stmp2, pd));
  printf("%s %s\n", buff, airSprintSize_t(stmp1, strlen(buff)));

  exit(0);
}



