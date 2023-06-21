#ifndef REP_INFO_H
#define REP_INFO_H

typedef struct rep_info_t {    
	int replicas_total_number;
  int is_evenodd;
	int defect_boundary;
	int defect_coordinates[3];
	double *cr_vec;
	int *label;
} rep_info;
extern rep_info *rep;

#endif
