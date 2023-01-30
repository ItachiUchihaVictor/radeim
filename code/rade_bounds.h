/**
 *
 * Copyright 2014 Matteo Riondato <matteo@cs.brown.edu>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * 
 */

#ifndef RADE_BOUNDS_H
#define RADE_BOUNDS_H

long double get_max_deviation_bound(const long double rademacher_bound);
long double get_rademacher_bound_standard();
long double get_rademacher_bound_optimize();
long double get_rademacher_bound_optimize_mathematica();
long double get_rademacher_bound_refined();
long double get_rademacher_bound_new();
#endif
