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

#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <map>
#include <random>
#include <set>
#include <vector>

#include "rade.h"
#include "rade_sample.h"
#include "stats.h"
#include "util.h"

using namespace std;
using namespace std::chrono;

// The key is a transaction, the value the number of times it was sampled
map<set<unsigned int>, unsigned int> all_sampled_transactions;

unsigned long sample(FILE __restrict__ *ds_FILE, unsigned int to_sample) {
	static uniform_int_distribution<int> distribution(0, dataset_size - 1);

	auto sampling_start_time = steady_clock::now();
	// Add transactions to sample
	vector<unsigned int> sampled_indices(to_sample, 0);
	for (unsigned int i = 0; i < to_sample; ++i) {
		sampled_indices[i] = distribution(generator);
	}
	sort(sampled_indices.begin(), sampled_indices.end());

	unsigned int transaction_index = -1;
	unsigned int sampled_index = 0;
	char line[TRANSACTION_LINE_MAXLEN];
	while (sampled_index < to_sample) {
		// Read line
		do {
			get_next_line(line, sizeof(line), ds_FILE);
			transaction_index++;
		} while (transaction_index < sampled_indices[sampled_index]);

		unsigned int base_transaction_id = sampled_indices[sampled_index];

		// Update data structures
		set<unsigned int> transaction = transaction_to_set(line);
		// Only update the structures if this is the first time
		// we see this transaction
		if (all_sampled_transactions.count(transaction) == 0) {
			all_sampled_transactions[transaction] = 0u;
			for (unsigned int item : transaction) {
				if (sample_items_sample_as_set.count(item) == 0) {
					sample_items_sample_as_set[item] = 1;
				} else {
					sample_items_sample_as_set[item]++;
				}
			}
			if (transaction.size() > max_transaction_length) {
				max_transaction_length = transaction.size();
			}
		}

		// Save the sampled line (possibly multiple times)
		do {
			all_sampled_transactions[transaction]++;
			sampled_index++;
			for (unsigned int item : transaction) {
				if (sample_items_sample_as_bag.count(item) == 0) {
					sample_items_sample_as_bag[item] = 1u;
				} else {
					sample_items_sample_as_bag[item]++;
				}
			}
		} while (sampled_index < to_sample &&
				sampled_indices[sampled_index] ==
				base_transaction_id);
	}
	auto sampling_end_time = steady_clock::now();
	return duration_cast<milliseconds>(sampling_end_time -
			sampling_start_time).count();
}

void __attribute__((cold)) print_sample() {
	for (auto element: all_sampled_transactions) {
		string transaction = itemset_to_str(element.first);
		for (unsigned int i = 0; i < element.second; ++i) {
			cout << transaction << endl;
		}
	}
}
