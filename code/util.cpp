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

#ifdef DEBUG
#include <cassert>
#include <cfloat>
#endif
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <map>
#include <set>
#include <sstream>

#include "util.h"

using namespace std;

/**
 * Convert an itemset to a string of items separated by whitespaces.
 *
 */
string __attribute__((pure)) itemset_to_str(const set<unsigned int>& itemset)  {
	stringstream stream;
	for (unsigned int item : itemset) {
		stream << item << " ";
	}
	stream.unget(); // remove last space
	return stream.str();
}

/**
 * Return the logarithm of (n choose k)
 *
 * This can probably be optimized a lot in terms of numerical stability/accuracy
 */
long double __attribute__((pure)) log_of_binomial(unsigned int n, unsigned int k) {
	if (n < k) {
		return -1.0; // special.
	}
	if (n == k) {
		return 0.0;
	}

	unsigned int first_index;
	unsigned int second_index;
	if (k < n - k) {
		first_index = n - k + 1;
		second_index = k;
	} else {
		first_index = k + 1;
		second_index = n - k;
	}
	long double to_return = 0.0;
	for (unsigned int i = first_index; i <= n; ++i) {
		to_return += logl(i);
	}
	// Starting from 2 because log(1) = 0
	for (unsigned int i = 2; i <= second_index; ++i) {
		to_return -= logl(i);
	}
	return to_return;
}

/**
 * Compare sizes of two sets (transactions).
 *
 * Return true if i is larger than j.
 *
 */
bool __attribute__((pure)) set_is_larger(set<unsigned int> i, set<unsigned int> j) {
	return (i.size() > j.size());
}
