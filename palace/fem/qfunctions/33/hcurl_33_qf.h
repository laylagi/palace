// Copyright Amazon.com, Inc. or its affiliates. All Rights Reserved.
// SPDX-License-Identifier: Apache-2.0

#ifndef PALACE_LIBCEED_HCURL_33_QF_H
#define PALACE_LIBCEED_HCURL_33_QF_H

#include "../coeff/coeff_3_qf.h"
#include "utils_33_qf.h"

#include <stdio.h>


CEED_QFUNCTION(f_apply_hcurl_33)(void *__restrict__ ctx, CeedInt Q,
                                 const CeedScalar *const *in, CeedScalar *const *out)
{
  const CeedScalar *qdata = in[0], *u = in[1];
  CeedScalar *v = out[0];

  const CeedInt stride = 2 + 9; // attr, w * |J|, (adjJt / |J|) colwise

  // printf("Q: %d\nAs int:", Q);
  // CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  // {
  //   for (CeedInt j = 0; j < stride; j++)
  //     printf("qdata[%d][%d] %d ", i, j, (CeedInt)qdata[i * stride + j]);
  //   printf("\n");
  // }

  CeedPragmaSIMD for (CeedInt i = 0; i < Q; i++)
  {
    const CeedScalar *qdata_i = qdata + i * stride;
    const CeedScalar* adjJt_loc = qdata_i + 2;
    const CeedScalar u_loc[3] = {u[i + Q * 0], u[i + Q * 1], u[i + Q * 2]};
    CeedScalar coeff[9], v_loc[3];

    // printf("hello: %d\n", i);
    // printf("qdata_i[0] %d\n", (CeedInt)qdata_i[0]);
    // printf("qdata_i[1] %f\n", qdata_i[1]);
    CoeffUnpack3((const CeedIntScalar *)ctx, (CeedInt)qdata_i[0], coeff);
    // return 1;
    MultAtBCx33(adjJt_loc, coeff, adjJt_loc, u_loc, v_loc);

    v[i + Q * 0] = qdata_i[1] * v_loc[0];
    v[i + Q * 1] = qdata_i[1] * v_loc[1];
    v[i + Q * 2] = qdata_i[1] * v_loc[2];
  }
  return 0;
}

#endif  // PALACE_LIBCEED_HCURL_33_QF_H
