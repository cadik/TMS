// SPDX-License-Identifier: Apache-2.0
// 
// Copyright 2015 Conrad Sanderson (http://conradsanderson.id.au)
// Copyright 2015 National ICT Australia (NICTA)
// 
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
// http://www.apache.org/licenses/LICENSE-2.0
// 
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
// ------------------------------------------------------------------------


#include <armadillo>
#include "catch.hpp"

using namespace arma;


TEST_CASE("mat_neg_1")
  {
  mat A = 
    "\
     0.061198   0.201990   0.019678  -0.493936  -0.126745   0.051408;\
     0.437242   0.058956  -0.149362  -0.045465   0.296153   0.035437;\
    -0.492474  -0.031309   0.314156   0.419733   0.068317  -0.454499;\
     0.336352   0.411541   0.458476  -0.393139  -0.135040   0.373833;\
     0.239585  -0.428913  -0.406953  -0.291020  -0.353768   0.258704;\
    ";
  
  mat B = -A;
  
  REQUIRE( B(0,0) == Approx(-0.061198) );
  REQUIRE( B(1,0) == Approx(-0.437242) );
  REQUIRE( B(2,0) == Approx(+0.492474) );
  REQUIRE( B(3,0) == Approx(-0.336352) );
  REQUIRE( B(4,0) == Approx(-0.239585) );
  
  REQUIRE( B(0,1) == Approx(-0.201990) );
  REQUIRE( B(1,1) == Approx(-0.058956) );
  REQUIRE( B(2,1) == Approx(+0.031309) );
  REQUIRE( B(3,1) == Approx(-0.411541) );
  REQUIRE( B(4,1) == Approx(+0.428913) );
  
  REQUIRE( B(0,5) == Approx(-0.051408) );
  REQUIRE( B(1,5) == Approx(-0.035437) );
  REQUIRE( B(2,5) == Approx(+0.454499) );
  REQUIRE( B(3,5) == Approx(-0.373833) );
  REQUIRE( B(4,5) == Approx(-0.258704) );
  
  REQUIRE( accu(B + A) == Approx(0.0).margin(0.001) );
  
  // REQUIRE_THROWS(  );
  }
