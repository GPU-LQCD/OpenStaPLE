#ifndef FIELD_CORR_H_
#define FIELD_CORR_H_

#include "./geometry.h"
#include "./su3_utilities.h"

void calc_field_corr_single_orientation(
										 __restrict su3_soa * const u,
										 __restrict su3_soa * const field_corr,
										 su3_soa * field_corr_aux,
										 su3_soa * loc_plaq,
										 dcomplex_soa * const trace_local,
										 d_complex * const corr,
										 __restrict single_su3 * const closed_corr,
										 const int mu, const int nu, const int ro);

void random_gauge_transformation(__restrict su3_soa * const u,
																 __restrict su3_soa * const u_new,
																 __restrict su3_soa * const m_soa);
#endif
