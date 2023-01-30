/**
 *
 * copyright 2014 matteo riondato <matteo@cs.brown.edu>
 *
 * licensed under the apache license, version 2.0 (the "license");
 * you may not use this file except in compliance with the license.
 * you may obtain a copy of the license at
 *
 *      http://www.apache.org/licenses/license-2.0
 *
 * unless required by applicable law or agreed to in writing, software
 * distributed under the license is distributed on an "as is" basis,
 * without warranties or conditions of any kind, either express or implied.
 * see the license for the specific language governing permissions and
 * limitations under the license.
 * 
 */

#ifndef RADE_SAMPLE_H
#define RADE_SAMPLE_H

#include <random>
#include <cstdio>

extern mt19937_64 generator;
unsigned long sample(FILE __restrict__ *ds_FILE, unsigned int to_sample);
void print_sample();
#endif
