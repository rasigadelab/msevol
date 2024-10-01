#pragma once
#include "../msmodel/universe.h"
//using namespace hgraph;

/* 5 elements of each type, unassigned */
void setup_00(Universe& u);

/* Gene -x2-> Chrom -x1-> Cell -x100-> Patch */
void setup_01(Universe& u);

/* Chrom0 -x1-> (Cell x5) -x100-> (Patch x5)*/
void setup_02(Universe& u);