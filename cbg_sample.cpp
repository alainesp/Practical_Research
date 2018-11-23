///////////////////////////////////////////////////////////////////////////////
// Sample code using Cuckoo Breeding Ground (CBG) hashtable
///////////////////////////////////////////////////////////////////////////////
//
// Written by Alain Espinosa <alainesp at gmail.com> in 2018 and placed
// under the MIT license (see LICENSE file for a full definition).
//
///////////////////////////////////////////////////////////////////////////////

// Only need to include and copy this file
#include "cbg.hpp"

///////////////////////////////////////////////////////////////////////////////
// The most common options for sets are:
///////////////////////////////////////////////////////////////////////////////
using set_negative_fast		= cbg::Set_SoA<2, uint64_t>;// Set, faster negative queries, fastest (recommended load_factor < 60%)
using set_negative_fat		= cbg::Set_SoA<4, uint64_t>;// Set, faster negative queries, little memory waste (when high load_factor, can reach 99%)
using set_negative_balanced	= cbg::Set_SoA<3, uint64_t>;// Set, faster negative queries, balanced (recommended 60% < load_factor < 95%)

using set_positive_fast		= cbg::Set_AoS<2, uint64_t>;// Set, faster positive queries, fastest (recommended load_factor < 60%)
using set_positive_fat		= cbg::Set_AoS<4, uint64_t>;// Set, faster positive queries, little memory waste (when high load_factor, can reach 99%)
using set_positive_balanced	= cbg::Set_AoS<3, uint64_t>;// Set, faster positive queries, balanced (recommended 60% < load_factor < 95%)
// Similar for maps, one example:
// ... other maps ...
using map_positive_fast = cbg::Map_AoS<2, std::string, uint64_t>;// Map, faster positive queries, fast (recommended load_factor < 60%)
// ... other maps ...
///////////////////////////////////////////////////////////////////////////////


// Other include for benchmarking
#include <chrono>
void benchmark()
{
	const size_t MAX_BUCKETS = 100'000;
	const size_t MAX_ELEMS = MAX_BUCKETS / 100 * 70;// load_factor 70%
	const size_t MAX_REPEAT = 200;

	// Seed with a real random value, if available
	std::random_device good_random;
	uint32_t seed = good_random();
	std::mt19937_64 r(seed);
	std::mt19937_64 r_positive(seed);
	auto elapsed_insert = std::chrono::nanoseconds::zero();
	auto elapsed_positive = std::chrono::nanoseconds::zero();

	// CBG hash table
	set_negative_fast cuckoo_table(MAX_BUCKETS);

	// Fill the table
	for (size_t repeat = 0; repeat < MAX_REPEAT; repeat++)
	{
		cuckoo_table.clear();

		// Insert
		auto start = std::chrono::high_resolution_clock::now();
		for (size_t i = 0; i < MAX_ELEMS; i++)
			cuckoo_table.insert(r());
		auto end = std::chrono::high_resolution_clock::now();
		elapsed_insert += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

		// Positive lookup
		uint32_t num_found = 0;
		start = std::chrono::high_resolution_clock::now();
		for (size_t i = 0; i < MAX_ELEMS; i++)
		{
			uint64_t val = r_positive();
			if (cuckoo_table.count(val))
				num_found++;
		}

		end = std::chrono::high_resolution_clock::now();
		elapsed_positive += std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);

		if (num_found != cuckoo_table.size())
			printf("Error in positive lookup %u of %u\n", num_found, uint32_t(cuckoo_table.size()));
	}

	printf(" Avg  Insertion time per elem: %u ns\n", uint32_t(elapsed_insert.count() / cuckoo_table.size() / MAX_REPEAT));
	printf("Positive lookup time per elem: %u ns\n", uint32_t(elapsed_positive.count() / cuckoo_table.size() / MAX_REPEAT));

	// Negative lookup
	auto start = std::chrono::high_resolution_clock::now();
	uint32_t num_found = 0;

	for (size_t repeat = 0; repeat < MAX_REPEAT; repeat++)
	{
		for (size_t i = 0; i < MAX_ELEMS; i++)
			if (cuckoo_table.count(r()))
				num_found++;

		if (num_found)
			printf("This is probably an error on negative lookup\n");
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end - start);
	printf("Negative lookup time per elem: %u ns\n", uint32_t(elapsed.count() / cuckoo_table.size() / MAX_REPEAT));
}

void main()
{
	benchmark();

	// Wait for one keystroke
	printf("\nPress any key to exit...");
	char c;
	scanf_s("%c", &c, 1);
}