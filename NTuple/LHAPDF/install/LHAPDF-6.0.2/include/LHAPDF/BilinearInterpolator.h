// -*- C++ -*-
//
// This file is part of LHAPDF
// Copyright (C) 2012-2013 The LHAPDF collaboration (see AUTHORS for details)
//
#pragma once

#include "LHAPDF/Interpolator.h"

namespace LHAPDF {


  /// Implementation of bilinear interpolation
  class BilinearInterpolator : public Interpolator {
  public:
    double _interpolateXQ2(const KnotArray1F& subgrid, double x, size_t ix, double q2, size_t iq2) const;
  };


}
