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

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <iostream>
#include <map>

#include "rade.h"
#include "rade_bounds_ci.h"
#include "util.h"

using namespace std;

/**
 * Compute an upper bound to the logarithm of the number of closed itemsets
 * using the lengths of the transactions in the sample
 */
long double get_log_CI_bound_lengths(const map<unsigned int, unsigned int>& items) {
	long double bound = log(items.size());
	for (auto element : all_sampled_transactions) {
		unsigned int size = 0u;
		if (items.size() < sample_items_sample_as_set.size()) {
			for (unsigned int item : element.first) {
				if (items.count(item) > 0) {
					++size;
				}
			}
		} else {
			size = (element.first).size();
		}
		if (size > 1u) {
			long double log_B = M_LN2 * size;
			long double log_to_remove = log1p(size); // "1p" for the empty set.
			log_B = logdiff(log_B, log_to_remove);
			bound = logsum(bound, log_B);
		}
	}
	return bound;
}

/**  
 * Compute an upper bound to the logarithm of the number of closed itemsets
 * using the lengths of the transactions in the sample and the support of the
 * items in the sample.
 */
long double get_log_CI_bound_refined(const map<unsigned int, unsigned int>& items) {
	// All the items are potential CIs, so consider them.
	long double bound = log(items.size());
	// For the following map:
	// key: item 'k', value: a map with key an integer 'ell' and with value the number
	// of sampled transactions containing the item 'k' and 'ell' items "larger" than
	// k. An item 'a' is "larger" than an item 'b' if the support of 'a' is
	// greater than the support of 'b' or, if supp(a)==supp(b), if a > b. This
	// implies a total order in the items.
	// FIXME This is really not a good name...
	map<unsigned int,map<unsigned int, unsigned int> > transaction_lengths;
	// Populate the above map
	for (auto element : all_sampled_transactions) {
		set<unsigned int> transaction_set;
		if (items.size() < sample_items_sample_as_set.size()) {
			// Filter out items that are not in items
			for (unsigned int item : element.first) {
				if (items.count(item) > 0) {
					transaction_set.insert(item);
				}
			}
			// Pass to the next transaction if no items from items are in
			// this transaction.
			if (transaction_set.size() == 0) {
				continue;
			}
		} else {
			transaction_set = element.first;
		}
		
		//
		vector<unsigned int> transaction(transaction_set.begin(),
				transaction_set.end());
		// Sort the items in the transaction by reverse order according to
		// frequency and, in case of equal frequency, based on their value
		// (larger before smaller).
		sort(transaction.begin(), transaction.end(), [&items](unsigned
					int a, unsigned int b)->bool { return (items.at(a) >
						items.at(b) || (items.at(a) ==
							items.at(b) && a > b));});
		// In the following loop we start from 1 because there are no potential
		// CIs to be counted for the "largest" item in the transaction: they have
		// all been accounted for when dealing with "smaller" items. Therefore,
		// it is useless to keep track of the number of transactions where an
		// item is the largest.
		for (unsigned int i = 1u; i < transaction.size(); ++i) { 
			unsigned int item = transaction[i];
			if (transaction_lengths.count(item) == 0) {
				map<unsigned int, unsigned int> item_map;
				transaction_lengths[item] = item_map;
			}
			if (transaction_lengths[item].count(i) == 0) {
				transaction_lengths[item][i] = 0u;
			}
			transaction_lengths[item][i]++;
		}
	} // end for(auto element : all_sampled_transactions)

	// Now, for each item 'a', and each "length" 'ell' in the map of 'a', we first
	// compute the number of transactions with length greater than 'ell', then
	// we "virtually" impose a total order between the transactions of length at
	// least ell (this is to avoid multiple counting) and, for each of these
	// transactions, and compute for each transaction 't',
	// the number of potential CIs containing 'a' that can be built using 't'
	// and the transactions "longer" than t.
	// This number is the minimum between 2^{number of "transactions" longer
	// than t} and 2^{ell} (ell is the length of t).
	for (auto element : transaction_lengths) {
		unsigned int trans_of_larger_length = 0u;
		for (auto rit = (element.second).rbegin(); rit !=
				(element.second).rend(); ++rit) {
			for (unsigned int i = 0u; i < rit->second; ++i) {
				bound = logsum(bound, min(trans_of_larger_length + i, rit->first) * M_LN2);
			}
			trans_of_larger_length += rit->second;
		}
	}
	return bound;
}

/**
 * Compute an upper bound to the logarithm of the number of closed itemsets
 * using the number of items in the sample and the maximum length of a
 * transaction in the sample
 */
long double get_log_CI_bound_items(const map<unsigned int, unsigned int>& items) {
	long double bound = log(items.size());
	unsigned int max_length = 0u;
	if (items.size() < sample_items_sample_as_set.size()) {
		for (auto element : all_sampled_transactions) {
			unsigned int items_in_transaction = 0u;
			for (unsigned int item : element.first) {
				if (items.count(item) > 0) {
					++items_in_transaction;
				}
			}
			if (items_in_transaction > max_length) {
				max_length = items_in_transaction;
			}
		}
	} else {
		max_length = max_transaction_length;
	}
	for (unsigned int i = 2u; i <= max_length; ++i) {
		long double log_B = log_of_binomial(items.size(), i);
		if (log_B >= 0.0) {
			bound = logsum(bound, log_B);
		}
	}
	return bound;
}

/**
 * Compute upper bound to empirical VC-dimension. This is basically a copy of
 * get_d_index_bound from stats.h .
 *
 * The algorithm to compute the upper bound to the d-index is Algorithm 1 from
 * M. Riondato and E. Upfal, Efficient Discovery of Association Rules and
 * Frequent Itemsets through Sampling with Tight Performance Guarantees, ACM
 * Transactions on Knowledge Discovery from Data (in press).
 */
unsigned int get_eVC_bound(const map<unsigned int, unsigned int>& items) {
	unsigned int eVC_bound = 0u;
	
	// This is to store the long transactions 
	vector<set<unsigned int>> long_transactions_d; 

	// Read first transaction
	auto iterator = all_sampled_transactions.begin();
	if (items.size() < sample_items_sample_as_set.size()) {
		set<unsigned int> first_transaction_set;
		do {
			for (unsigned int item : iterator->first) {
				if (items.count(item) > 0) {
					first_transaction_set.insert(item);
				}
			}
			++iterator;
		} while (first_transaction_set.size() == 0 && iterator != all_sampled_transactions.end());
		if (first_transaction_set.size() > 0) {
			long_transactions_d.push_back(first_transaction_set);
			eVC_bound = 1u;
		}
	} else {
		long_transactions_d.push_back(iterator->first);
		eVC_bound = 1u;
		++iterator;
	}

	// Iterate over all transactions to compute bound to d-index
	while (iterator != all_sampled_transactions.end())  {
		set<unsigned int> curr_transaction;
		if (items.size() < sample_items_sample_as_set.size()) {
			for (unsigned int item : iterator->first) {
				if (items.count(item) > 0) {
					curr_transaction.insert(item);
				}
			}
			if (curr_transaction.size() == 0) {
				++iterator;
				continue;
			}
		} else {
			curr_transaction = iterator->first;
		}
		++iterator; 

		unsigned int curr_transaction_len = curr_transaction.size();

		// Only process if longer than candidate d-index
		if (curr_transaction_len > eVC_bound) {
			// We now check whether this transaction is a subset of
			// one of the long transactions we stored
			bool process = true;
			for (set<unsigned int> long_transaction :
					long_transactions_d) {
				if (includes(long_transaction.begin(),
						long_transaction.end(),
						curr_transaction.begin(),
						curr_transaction.end()))  {
					process = false;
					break;
				}
			}
			// If it was a subset of a long transaction, skip it and
			// go to the next one
			if (!process) {
				continue; // iterate the main while loop
			}
			// The transaction is longer than d-index and not a
			// subset of a longer transaction. Add it to the list of
			// long transactions
			long_transactions_d.push_back(curr_transaction);
			// Sort the long transactions by size in reverse order
			sort(long_transactions_d.begin(), long_transactions_d.end(),
					set_is_larger);
			// Compute new candidate d-index
			if (long_transactions_d[eVC_bound].size() > eVC_bound) {
				++eVC_bound;
			}
			// Trim the vector of long transactions
			long_transactions_d.resize(eVC_bound);
		}
	}

	if (verbose) {
		cerr << VERBOSE_HEADER << "eVC_bound: " << eVC_bound << endl;
	}
	return eVC_bound;
}

/**
 * Compute an upper bound to the logarithm of the number of closed itemsets
 * using the size of the sample
 */
long double get_log_CI_bound_size(const map<unsigned int, unsigned int>& items) { 
	unsigned int local_sample_size = 0u;
	if (items.size() < sample_items_sample_as_set.size()) {
		for (auto element : all_sampled_transactions) {
			for (unsigned int item : element.first) {
				if (items.count(item) > 0) {
					++local_sample_size;
					break;
				}
			}
		}
	} else {
		local_sample_size = sample_size;
	}
	return local_sample_size * M_LN2;
}

/**
 * Compute an upper bound to the logarithm of the number of closed itemsets
 * using the empirical VC-dimension. See Thm. 3.7 in Anthony and Bartlett,
 * Neural Network Learning: Theoretical Foundations, Cambridge University Press,
 * 1999)
 */
long double get_log_CI_bound_eVC(const map<unsigned int, unsigned int>& items) {
	unsigned int eVC_bound = get_eVC_bound(items);
	return (eVC_bound <= sample_size) ?  eVC_bound * (1 + logl(sample_size) -
			logl(eVC_bound)) : eVC_bound * logl(sample_size + 1);
}

/**
 * Return an upper bound to the logarithm of the number of closed itemsets built
 * on the items in 'items', using the selected methods to compute the bound.
 */
long double get_log_CI_bound(const map<unsigned int, unsigned int>& items) {
	long double log_CI_bound_lengths = (use_lengths_CI_method) ?
		get_log_CI_bound_lengths(items) : LDBL_MAX;
	long double log_CI_bound_items = (use_items_CI_method) ?
		get_log_CI_bound_items(items) : LDBL_MAX;
	long double log_CI_bound_eVC = (use_eVC_CI_method) ? get_log_CI_bound_eVC(items) :
		DBL_MAX;
	long double log_CI_bound_size = (use_size_CI_method) ?
		get_log_CI_bound_size(items) : LDBL_MAX;
	long double log_CI_bound_refined = (use_refined_CI_method) ? 
		get_log_CI_bound_refined(items) : LDBL_MAX;
	long double log_CI_bound = min(min(min(min(log_CI_bound_lengths,
							log_CI_bound_items),
						log_CI_bound_eVC),
					log_CI_bound_size), log_CI_bound_refined);
	if (verbose) {
		bool printed_something = false;
		cerr << VERBOSE_HEADER;
		if (use_lengths_CI_method) {
			cerr << "log_CI_bound_lengths: " << log_CI_bound_lengths;
			printed_something = true;
		}
		if (use_items_CI_method) {
			if (printed_something) {
				cerr << ", ";
			}
			cerr << "log_CI_bound_items: " << log_CI_bound_items;
			printed_something = true;
		} 
		if (use_eVC_CI_method) {
			if (printed_something) {
				cerr << ", ";
			} 
			cerr << "log_CI_bound_eVC: " << log_CI_bound_eVC; 
			printed_something = true;
		}
		if (use_size_CI_method) {
			if (printed_something) {
				cerr << ", ";
			} 
			cerr << "log_CI_bound_size: " << log_CI_bound_size;
			printed_something = true;
		}
		if (use_refined_CI_method) {
			if (printed_something) {
				cerr << ", ";
			} 
			cerr << "log_CI_bound_refined: " << log_CI_bound_refined;
		}
		if (use_lengths_CI_method + use_items_CI_method + use_eVC_CI_method +
				use_size_CI_method + use_refined_CI_method > 1) {
			cerr << ", minimum: ";
			if (log_CI_bound == log_CI_bound_lengths) {
				cerr << "lengths";
			} else if (log_CI_bound == log_CI_bound_items) {
				cerr << "items";
			} else if (log_CI_bound == log_CI_bound_eVC) {
				cerr << "eVC";
			} else if (log_CI_bound == log_CI_bound_size) {
				cerr << "size";
			} else {
				cerr << "refined";
			}
			 cerr << " (" << log_CI_bound << ")";
		}
		cerr << endl;
	}
	return log_CI_bound;
}

