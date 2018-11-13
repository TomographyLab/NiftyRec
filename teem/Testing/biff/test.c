/*
  Teem: Tools to process and visualize scientific data and images              
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


#include "teem/biff.h"

int
main() {
  char *tmp;

  biffAdd("axis", "the first error axis");
  biffAdd("axis", "the second error axis");
  biffAdd("axis", "the third error axis");
  biffAdd("chard", "the first error chard");
  biffAdd("chard", "the second error chard");
  biffAdd("chard", "the third error chard");
  biffAdd("bingo", "zero-eth bingo message");
  biffMove("bingo", NULL, "chard");
  biffAdd("bingo", "the first error bingo");
  biffAdd("bingo", "the second bll boo boo boo error bingo");
  biffAdd("bingo", "the third error bingo");
  printf("%s\n", (tmp = biffGet("bingo")));
  free(tmp);
  biffDone("bingo");
  printf("%s\n", (tmp = biffGet("chard")));
  free(tmp);
  biffDone("chard");
  printf("%s\n", (tmp = biffGet("axis")));
  free(tmp);
  biffDone("axis");
  
  biffAdd("harold", "the first error harold");
  biffAdd("harold", "the second error harold");
  biffAdd("harold", "the third error harold");
  printf("%s\n", (tmp = biffGet("harold")));
  free(tmp);

  biffAdd("axis", "the first error axis");
  biffAdd("axis", "the second error axis");
  biffAdd("axis", "the third error axis");
  biffAdd("axis", "the fourth error axis");
  biffAdd("axis", "the fifth error axis");
  printf("%s\n", (tmp = biffGet("axis")));
  free(tmp);
  biffDone("axis");

  biffAddf("test", "%s: this is a test of biffAddf %d %f", "me", 1, 2.0);
  printf("%s\n", (tmp = biffGet("test")));
  free(tmp);
  biffDone("test");

  biffAddf("test2", "%s: this is a test of biffAddf %d %f", "me", 1, 2.0);
  biffMovef("test3", "test2", "%s: testing biffMove %d.", "me", 69);
  printf("%s\n", (tmp = biffGet("test3")));
  free(tmp);
  biffDone("test3");

  exit(0);
}



