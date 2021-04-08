//
//  OTORModel.hpp
//  GlowCurveAnalsys
//
//  Created by jeremy hepker on 7/14/19.
//  Copyright © 2019 Jeremy Hepker. All rights reserved.
//

#ifndef OTORModel_hpp
#define OTORModel_hpp

#include <stdio.h>
#include <cmath>
#include <vector>
#include <boost/math/special_functions/lambert_w.hpp>
#include <boost/math/special_functions/expint.hpp>

using boost::math::lambert_w0;
using boost::math::expint;

void OTORModel(std::vector<double>& x, std::vector<double>& peak, double Tm, double Im, double E);

#endif /* OTORModel_hpp */
