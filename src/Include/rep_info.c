#ifndef REP_INFO_C
#define REP_INFO_C

#include "./rep_info.h"

rep_info *rep;

int get_index_of_pbc_replica(){
  // finds index corresponding to label=0
  int ridx_lab0;
  for(ridx_lab0=0; 0!=rep->label[ridx_lab0]; ++ridx_lab0){}

  return ridx_lab0;
}

#endif
