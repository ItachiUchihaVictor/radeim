/**
 * Useful functions
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

#ifndef UTIL_H
#define UTIL_H

#ifdef DEBUG
#include <cassert>
#include <cfloat>
#endif
#include <cmath>
#include <cstring>
#include <map>
#include <set>

using namespace std;

static const string COMMA = ",";
static const string SPACE = " ";
static const string VERBOSE_HEADER = "INFO: ";
static const string LOG_HEADER = "LOG: ";
static const string ERROR_HEADER = "ERROR: ";

// Allowing for very long transactions, just in case
static const short TRANSACTION_LINE_MAXLEN = 10000;

void *malloc_check(size_t size, const char *err_msg);
string itemset_to_str(const set<unsigned int>& itemset);
set<unsigned int> transaction_to_set(const char *line, bool is_utility_line=false);
string sort_transaction(const char *line);
bool get_next_line(char *line, int size, FILE *in_file);
long double log_of_binomial(unsigned int n, unsigned int k);
bool set_is_larger(set<unsigned int> i, set<unsigned int> j);

inline long double logsum(const long double log_a, const long double log_b) {
#ifdef DEBUG
	long double to_return = max(log_a, log_b) + log1pl(expl(fminl(log_a, log_b) - max(log_a, log_b)));
	assert(to_return > max(log_a, log_b) || to_return < nexttoward(max(log_a,log_b), LDBL_MAX));
	return to_return;
#else
	return max(log_a, log_b) + log1pl(expl(min(log_a, log_b) - max(log_a, log_b)));
#endif
}

/**
 * Call and check malloc. Print message to stderr and exit if fails.
 *
 */
inline void *malloc_check(size_t size, const char *err_msg) {
	void *ptr = malloc(size);
	if (ptr == NULL) {
		cerr << "CRITICAL: malloc() failed in " << err_msg << endl;
		exit(EXIT_FAILURE);
	}
	return ptr;
}

/**
 * Return the logarithm of "a-b" from log(a) and log(b) using the identities
 * log(a-b) = log(a) + log(1 + b/a) and b/A = exp(log(b) - log(a))
 * Assumes that a > b
 */
inline long double logdiff(const long double log_a, const long double log_b) {
	return log_a + log1pl(-expl(log_b - log_a));
}


inline unsigned int get_max_item_supp(const map<unsigned int, unsigned int>& items) {
	unsigned int max_item_supp = 0u;
	for (auto item_freq_pair : items) {
		if (item_freq_pair.second > max_item_supp) {
			max_item_supp = item_freq_pair.second;
		}
	}
	return max_item_supp;
}

/**
 * Sort the items in a transaction.
 *
 */
inline __attribute__((pure)) string sort_transaction(const char *line) {
	return itemset_to_str(transaction_to_set(line));
}
#endif

/**
 * Convert a transaction to a set of items.
 *
 */
inline set<unsigned int> __attribute__((pure)) transaction_to_set(const char *line, bool is_utility_line) {
	// We copy the line because we modify it.
	char *copy = (char *) malloc_check(strlen(line) + 1, "transaction_to_set");
	char *orig_copy = copy;
	strncpy(copy, line, strlen(line) + 1);
	set<unsigned int> items;
	if (__builtin_expect(!!(is_utility_line), 0)) {
		// If this is an utility dataset, consider only the first part of the transaction
		char *colons_location = strchr(copy, ':');
		if (colons_location != NULL) {
			*colons_location = '\0';
		}
	}
	char *token;
	while ((token = strsep(&copy, " \t")) != NULL) {
		if (*token != '\0') {
			items.insert(strtol(token, NULL, 10));
		}
	}
	free(orig_copy);
	return items;
}

/**
 * Get next line from in_file.
 *
 * The line is copied in line. Only at most size-1 characters are copied.
 *
 * Metadata / comments (lines that starts with '#') and empty lines are skipped.
 *
 */
inline bool get_next_line(char *line, int size, FILE __restrict__ *in_file) {
	char c;
	// Skip metadata and empty lines
	while (true) {
		c = fgetc(in_file);
		if (c == '#') {
			fgets(line, size, in_file);
		} else if (c == '\n') {
			// do nothing
		} else {
			break;
		}
	}
	// Check for end of file
	if (feof(in_file)) {
		return false;
	}

	// Read line
	ungetc(c, in_file);
	fgets(line, size, in_file);

	return true;
}
