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
#include <cmath>
#include <iostream>
#include <map>
#include <vector>

#include <nlopt.h>

#include "rade.h"
#include "rade_bounds.h"
#include "rade_bounds_ci.h"
#include "util.h"

using namespace std;

/**
 *
 * Get the quantity that bounds the maximum deviation
 */
long double  __attribute__((pure)) get_max_deviation_bound(const long double rademacher_bound) {
	static double fixed_part = sqrt(2.0 * log(2.0 / delta));
	return 2.0 * rademacher_bound + fixed_part / sqrt(sample_size);
}

/**
 * Compute a bound to the empirical Rademacher average using a bound on the
 * shatter coefficient and the standard proof
 */
long double get_rademacher_bound_standard() {
	unsigned int max_item_supp = get_max_item_supp(sample_items_sample_as_bag);
	if (verbose) {
		cerr << VERBOSE_HEADER <<
			"rademacher_bound_standard max_item_supp: "
			<< max_item_supp << endl;
	}
	// Compute bounds to the number of Closed Itemsets
	long double log_CI_bound = get_log_CI_bound(sample_items_sample_as_set);
	return sqrtl(max_item_supp * 2.0 * log_CI_bound) / sample_size;
}

/*
 * Compute the optimization function and the gradient
 */
double opt_fun(unsigned, const double *x, double *grad, void *supps_quants_pointer) {
	vector<pair<long double, long double> > *supps_logquants = (vector<pair<long double, long double> > *) supps_quants_pointer;

	long double log_sum_square_none_none = -1.0;
	long double log_sum_square_none_lin = -1.0;
	auto iter = (*supps_logquants).begin();
	long double s_square = powl(x[0], 2.0); //s^2
	long double exponent = s_square * iter->first + iter->second;
	log_sum_square_none_none = exponent;
	log_sum_square_none_lin = exponent + iter->second + logl(iter->first);
	++iter;
	while (iter != (*supps_logquants).end()) {
		exponent = s_square * iter->first;
		log_sum_square_none_none = logsum(log_sum_square_none_none, exponent + iter->second);
		log_sum_square_none_lin = logsum(log_sum_square_none_lin, exponent + iter->second + logl(iter->first));
		++iter;
	}
	// XXX NLOpt only uses doubles, not long doubles :(
	if (grad != NULL) {
		grad[0] = (double) (expl(log_sum_square_none_lin + logl(2.0) - log_sum_square_none_none) - (log_sum_square_none_none / powl(x[0], 2.0)));
	}
	//if (verbose) {
	//	cerr << "x: " << x[0] << ", val: " << (double) (log_sum_square_none_none / x[0]) << " grad: " << ((grad) ? grad[0] : -1) << endl;
	//}
	return (double) (log_sum_square_none_none / x[0]);
}

/*
 * Compute a bound to the empirical Rademacher average using a new refined
 * procedure that uses optimization.
 */
long double get_rademacher_bound_optimize() {
	vector<pair<long double, long double> > supps_logquants;
	supps_logquants.clear();
	long double factor = 2.0 * powl(sample_size, 2.0);
	// All the items are potential CIs, so add their contribution.
	for (auto element : sample_items_sample_as_set) {
		pair<long double, long double> to_add(((long double) element.second) / factor, 0.0);
		supps_logquants.push_back(to_add);
	}
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
		// Filter out items that are not in sample_items
		for (unsigned int item : element.first) {
			if (sample_items_sample_as_set.count(item) > 0) {
				transaction_set.insert(item);
			}
		}
		// Pass to the next transaction if no items from sample_items are in
		// this transaction.
		if (transaction_set.size() == 0) {
			continue;
		}

		//
		vector<unsigned int> transaction(transaction_set.begin(),
				transaction_set.end());
		// Sort the items in the transaction by reverse order according to
		// frequency and, in case of equal frequency, based on their value
		// (larger before smaller).
		sort(transaction.begin(), transaction.end(),
				[](unsigned int a, unsigned int
					b)->bool { return (sample_items_sample_as_set.at(a) >
						sample_items_sample_as_set.at(b) ||
						(sample_items_sample_as_set.at(a) ==
							sample_items_sample_as_set.at(b) && a > b));});
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
	}
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
					pair<long double, long double> to_add(((long double)
						sample_items_sample_as_set.at(element.first)) / factor,
							min(trans_of_larger_length + i, rit->first) *
							M_LN2);
					supps_logquants.push_back(to_add);
			}
			trans_of_larger_length += rit->second;
		}
	}

	// Setup and solve the optimization problem
	// FIXME We are using a derivative-free algorithm despite knowing how to compute
	// the gradient because there seem to be some problems in NLopt
	// implementation or in its interpretation of our gradient or in our
	// computation of the gradient (see opt_fun).
	// We would even know how to compute the Hessian (second derivative), but
	// we're not using that either.
	// Note that we are also using a "local" optimization algorithm because the
	// function is convex.
	nlopt_opt opt_prob = nlopt_create(NLOPT_LN_COBYLA, 1);
	nlopt_set_min_objective(opt_prob, opt_fun, &supps_logquants);
	double lb[1] = {1.0};
	nlopt_set_lower_bounds(opt_prob, lb);
	// XXX It is weird that we have to define a const variable and then use a
	// pointer to it. NLopt bug? It's not what's documented.
	const double xtol_abs = 1e-4;
	nlopt_set_xtol_abs(opt_prob, &xtol_abs);
	nlopt_set_ftol_abs(opt_prob, 1e-6);
	// Set initialization point
	double x[1];
#ifdef DEBUG
	x[0] = opt_param;
#else
	x[0] = (double) opt_minimizer;
#endif
	double min_value;

	nlopt_result result = nlopt_optimize(opt_prob, x, &min_value);
#ifdef DEBUG
	cerr << "mathematica minimizer: " << opt_minimizer << ", nlopt minimizer: " << x[0] << endl;
#endif
	if (result >= 0) {
		opt_minimizer = x[0];
		if (verbose) {
			cerr << VERBOSE_HEADER << "get_rademacher_bound_optimize: found minimum at "
				<< x[0] << ": " << min_value << ", return code: " << result <<
				endl;
		}
	} else  {
		cerr << "Optimization failed! return code: " << result << endl;
		exit(1);
	}
	nlopt_destroy(opt_prob);

	return min_value;
}

#ifdef DEBUG
/*
 * Compute a bound to the empirical Rademacher average using a new refined
 * procedure that uses optimization with Mathematica
 */
long double get_rademacher_bound_optimize_mathematica() {
	long double factor = 2.0 * powl(sample_size, 2.0);
	stringstream minimize;
	minimize << fixed << setprecision(15) << "Minimize[{Log[";
	// All the items are potential CIs, so add their contribution.
	for (auto element : sample_items_sample_as_set) {
		minimize << "Exp[(x^2)*" << ((long double) element.second) / factor << "]+";
	}
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
		// Filter out items that are not in sample_items
		for (unsigned int item : element.first) {
			if (sample_items_sample_as_set.count(item) > 0) {
				transaction_set.insert(item);
			}
		}
		// Pass to the next transaction if no items from sample_items are in
		// this transaction.
		if (transaction_set.size() == 0) {
			continue;
		}

		//
		vector<unsigned int> transaction(transaction_set.begin(),
				transaction_set.end());
		// Sort the items in the transaction by reverse order according to
		// frequency and, in case of equal frequency, based on their value
		// (larger before smaller).
		sort(transaction.begin(), transaction.end(),
				[](unsigned int a, unsigned int
					b)->bool { return (sample_items_sample_as_set.at(a) >
						sample_items_sample_as_set.at(b) ||
						(sample_items_sample_as_set.at(a) ==
							sample_items_sample_as_set.at(b) && a > b));});
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
	}
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
				minimize << "Exp[(x^2)*" << ((long double)
						sample_items_sample_as_set.at(element.first)) / factor
					<< "+" << min(trans_of_larger_length + i, rit->first) *
					M_LN2 <<
				"]+";
			}
			trans_of_larger_length += rit->second;
		}
	}
	minimize << "0]/x, x>0}, x]";
	MLPutFunction(mathlink, "EvaluatePacket", 1);
	MLPutFunction(mathlink, "ToString", 1);
	MLPutFunction(mathlink, "ToExpression", 1);
	MLPutString(mathlink, minimize.str().c_str());
	MLEndPacket(mathlink);

	int pkt;
	while ((pkt = MLNextPacket(mathlink)) && pkt != RETURNPKT) {
		MLNewPacket(mathlink);
	}
	long double min_value = LDBL_MAX;
	if(!pkt) {
		// TODO handle_error(mathlink); /* including calling MLClearError(link) */
	} else {
		const char *result;
		if (! MLGetString(mathlink, &result)) {
			// TODO
		} else {
			char *second_token;
			min_value = strtold(result + 1, &second_token);
			opt_minimizer = strtold(second_token + 8, NULL);
			MLReleaseString(mathlink, result);
		}
	}
	return min_value;
}
#endif

/**
 * Compute a bound to the empirical Rademacher average using a new refined procedure.
 */
long double get_rademacher_bound_refined() {
	long double bound = 0.0;
	long double factor =  powl(opt_param, 2.0) / (2.0 * pow(sample_size, 2.0));
	// All the items are potential CIs, so add their contribution.
	for (auto element : sample_items_sample_as_set) {
		long double to_add = factor * element.second;
		if (bound == 0.0) {
			bound = to_add;
		} else {
			bound = logsum(bound, to_add);
		}
	}
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
		// Filter out items that are not in sample_items
		for (unsigned int item : element.first) {
			if (sample_items_sample_as_set.count(item) > 0) {
				transaction_set.insert(item);
			}
		}
		// Pass to the next transaction if no items from sample_items are in
		// this transaction.
		if (transaction_set.size() == 0) {
			continue;
		}

		//
		vector<unsigned int> transaction(transaction_set.begin(),
				transaction_set.end());
		// Sort the items in the transaction by reverse order according to
		// frequency and, in case of equal frequency, based on their value
		// (larger before smaller).
		sort(transaction.begin(), transaction.end(),
				[](unsigned int a, unsigned int
					b)->bool { return (sample_items_sample_as_set.at(a) >
						sample_items_sample_as_set.at(b) ||
						(sample_items_sample_as_set.at(a) ==
							sample_items_sample_as_set.at(b) && a > b));});
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
	}
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
				bound = logsum(bound, factor * sample_items_sample_as_set.at(element.first) +
						min(trans_of_larger_length + i, rit->first) * M_LN2);
			}
			trans_of_larger_length += rit->second;
		}
	}
	return bound / opt_param;
}

/**
 * Compute a bound to the empirical Rademacher average using a new procedure.
 */
long double get_rademacher_bound_new() {
	long double bound = 0.0;

	unsigned int max_item_supp = get_max_item_supp(sample_items_sample_as_set);
	vector<unsigned int> items;
	for (auto element : sample_items_sample_as_set) {
		items.push_back(element.first);
	}
	// Sort the items in reverse order according to frequency and, in case of
	// equal frequency, based on their value (larger before smaller).
	sort(items.begin(), items.end(), [](unsigned
				int a, unsigned int b)->bool { return (sample_items_sample_as_set.at(a) >
					sample_items_sample_as_set.at(b) || (sample_items_sample_as_set.at(a) ==
						sample_items_sample_as_set.at(b) && a > b));});

	// Compute the index of the top group.
	unsigned int group_index = 0u;
	while (pow(2.0, group_index) < max_item_supp) {
		++group_index;
	}
	unsigned int group_max_supp = max_item_supp;
	unsigned int group_min_supp = (unsigned int) pow(2.0, group_index - 1u) + 1u;

	unsigned int item_index = 0u;
	map<unsigned int, unsigned int> curr_items;
	double factor = pow(opt_param,2.0) / (2.0 * pow(sample_size, 2.0));

	do {
		// Early termination if we're done for sure.
		if (item_index == items.size()) {
			break;
		}

		// Get the items with support at least group_min_supp
		unsigned int old_items_size = curr_items.size();
		unsigned int group_item_max_supp = 0u;
		while (item_index < items.size() && sample_items_sample_as_set.at(items[item_index]) >= group_min_supp) {
			curr_items[items[item_index]] = sample_items_sample_as_set.at(items[item_index]);
			if (curr_items[items[item_index]] > group_item_max_supp) {
				group_item_max_supp = curr_items[items[item_index]];
			}
			++item_index;
		}

		// Skip group if empty
		if (old_items_size == curr_items.size()) {
			// Update boundaries for the next group
			group_max_supp = group_min_supp - 1u;
			group_min_supp = (group_max_supp / 2u) + 1u;
			--group_index;
			continue;
		}

		long double curr_log_CI_bound = get_log_CI_bound(curr_items);
		long double log_to_sum = curr_log_CI_bound + (factor * group_item_max_supp);
		if (bound == 0.0) {
			bound = log_to_sum;
		} else {
			bound = logsum(bound, log_to_sum);
		}

		// Update boundaries for the next group
		group_max_supp = group_min_supp - 1u;
		group_min_supp = (group_max_supp / 2u) + 1u;

		--group_index;
	} while (group_index > 0u);

	return bound / opt_param;
}
