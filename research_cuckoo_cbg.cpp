///////////////////////////////////////////////////////////////////////////////
// Research code on Cuckoo Breeding Ground hash table
///////////////////////////////////////////////////////////////////////////////
//
// Written by Alain Espinosa <alainesp at gmail.com> in 2018 and placed
// under the MIT license (see LICENSE file for a full definition).
//
///////////////////////////////////////////////////////////////////////////////

#include <cstdint>
#include <tuple>
#include <chrono>
#include <random>
#include <memory>
#include <unordered_set>
#include <cassert>

///////////////////////////////////////////////////////////////////////////////
// Windowed Cucko hashing implemented as only one table. Contains 
// implementations for d=2 and l<=127
//
// This is for research only, see 'cbg.hpp' for production code.
///////////////////////////////////////////////////////////////////////////////
template<uint32_t NUM_ELEMS_BUCKET> class CBG_Set
{
private:
	uint64_t* data1;
	uint32_t* cache1;
	uint32_t num_buckets;
	uint32_t num_elems;

	/////////////////////////////////////////////////////////////////////
	// Utilities
	/////////////////////////////////////////////////////////////////////
	__forceinline uint64_t hash_elem(uint64_t elem) const noexcept
	{
		return elem;
	}
	__forceinline bool cmp_elems(uint64_t l, uint64_t r) const noexcept
	{
		return l == r;
	}
	// Given a value "word", produces an integer in [0,p) without division.
	// The function is as fair as possible in the sense that if you iterate
	// through all possible values of "word", then you will generate all
	// possible outputs as uniformly as possible.
	static __forceinline uint32_t fastrange32(uint32_t word, uint32_t p)
	{
		return (uint32_t)(((uint64_t)word * (uint64_t)p) >> 32);
	}
	/////////////////////////////////////////////////////////////////////
	// Cache coded utilities
	/////////////////////////////////////////////////////////////////////
	//
	// Unlucky  Bucket is   Element        
	// bucket   Reversed    Distance    Labels
	// -------  ----------  --------    -----
	//   |          |       |      |    |   |
	//   b1         b0       x3  x2     x1 x0
	// 0x0'00'00
	__forceinline uint32_t Get_Label(uint32_t pos) const noexcept
	{
		return cache1[pos] & 0xff;
	}
	__forceinline bool Is_Empty(uint32_t pos) const noexcept
	{
		return Get_Label(pos) == 0;
	}
	__forceinline void Set_Empty(uint32_t pos) noexcept
	{
		cache1[pos] &= 0x3'00'00;
	}
	__forceinline void Copy_Elem(uint32_t dest, uint32_t source) noexcept
	{
		cache1[dest] = (cache1[dest] & 0x3'00'00) | (cache1[source] & 0xffff);
	}
	__forceinline void UpdateFlag(uint32_t pos, uint32_t distance_to_base, bool is_reverse_item) noexcept
	{
		cache1[pos] = (cache1[pos] & 0x3'00'ff) | (is_reverse_item ? 0x0'80'00 : 0) | (distance_to_base << 8);
	}
	__forceinline void UpdateFlag(uint32_t pos, uint32_t distance_to_base, bool is_reverse_item, uint32_t label) noexcept
	{
		cache1[pos] = (cache1[pos] & 0x3'00'00) | (is_reverse_item ? 0x0'80'00 : 0) | (distance_to_base << 8) | label;
	}
	__forceinline bool Is_Item_In_Reverse_Bucket(uint32_t pos) const noexcept
	{
		return cache1[pos] & 0x0'80'00;
	}
	__forceinline uint32_t GetFlagDistance(uint32_t pos) const noexcept
	{
		return (cache1[pos] >> 8) & 0x7f;
	}
	__forceinline bool Is_Unlucky_Bucket(uint32_t pos) const noexcept
	{
		return cache1[pos] & 0x2'00'00;
	}
	__forceinline void Set_Unlucky_Bucket(uint32_t pos) noexcept
	{
		cache1[pos] |= 0x2'00'00;
	}
	__forceinline bool Is_Reversed_Window(uint32_t pos) const noexcept
	{
		return cache1[pos] & 0x1'00'00;
	}
	__forceinline void Set_Reversed(uint32_t pos) noexcept
	{
		cache1[pos] |= 0x1'00'00;
	}

	__forceinline void MoveElem(uint32_t dest, uint32_t orig)
	{
		data1[dest] = data1[orig];
	}
	/////////////////////////////////////////////////////////////////////
	// Insertion algorithm utilities
	/////////////////////////////////////////////////////////////////////
	std::pair<uint32_t, uint32_t> Calculate_Minimum(uint32_t bucket_pos) const noexcept
	{
		uint32_t minimum = Get_Label(bucket_pos);
		uint32_t pos = bucket_pos;

		for (uint32_t i = 1; minimum && i < NUM_ELEMS_BUCKET; i++)
		{
			uint32_t label_value = Get_Label(bucket_pos + i);
			if (minimum > label_value)
			{
				minimum = label_value;
				pos = bucket_pos + i;
			}
		}

		return std::make_pair(minimum, pos);
	}
	__forceinline uint32_t Belong_to_Bucket(uint32_t elem_pos) const noexcept
	{
		if (Is_Empty(elem_pos))
			return UINT32_MAX;

		return elem_pos + (Is_Item_In_Reverse_Bucket(elem_pos) ? NUM_ELEMS_BUCKET - 1 : 0) - GetFlagDistance(elem_pos);
	}
	uint32_t Count_Empty(uint32_t pos) const noexcept
	{
		uint32_t count = Is_Empty(pos) ? 1 : 0;

		for (uint32_t i = 1; i < NUM_ELEMS_BUCKET; i++)
			if (Is_Empty(pos + i))
				count++;

		return count;
	}
	uint32_t Count_Elems_In_Bucket_Non_Reversed(uint32_t bucket_pos) const noexcept
	{
		uint32_t count = 0;

		for (uint32_t i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			uint32_t pos = bucket_pos + i;

			if (!Is_Item_In_Reverse_Bucket(pos) && i == GetFlagDistance(pos))// Faster Belong_to_Bucket(elem_pos)
				count++;
		}

		return count;
	}
	std::pair<uint32_t, uint32_t> Count_Elems_In_Bucket_Outside_Range(uint32_t bucket_pos, uint32_t range_init) const noexcept
	{
		uint32_t count_outsize_range = 0;
		uint32_t count = 0;

		for (uint32_t i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			uint32_t pos = bucket_pos + i;

			if (!Is_Item_In_Reverse_Bucket(pos) && i == GetFlagDistance(pos))// Faster Belong_to_Bucket(elem_pos)
			{
				count++;
				if ((pos - range_init) >= NUM_ELEMS_BUCKET)// Faster check for: !(range_init <= pos < (range_init+NUM_ELEMS_BUCKET)))
					count_outsize_range++;
			}
		}

		return std::make_pair(count, count_outsize_range);
	}
	// TODO: Do this reversing maintaining elems near the entry bin
	void Reverse_Bucket(uint32_t bucket_pos) noexcept
	{
		Set_Reversed(bucket_pos);

		uint32_t j = NUM_ELEMS_BUCKET - 1;
		for (uint32_t i = NUM_ELEMS_BUCKET - 1; i < NUM_ELEMS_BUCKET; i--)
			if (Belong_to_Bucket(bucket_pos + i) == bucket_pos)// Elems belong to our bucket
			{
				for (; j < NUM_ELEMS_BUCKET && !Is_Empty(bucket_pos - j); j--)
				{
				}// Find empty space
				if (j < NUM_ELEMS_BUCKET)
				{
					UpdateFlag(bucket_pos - j, NUM_ELEMS_BUCKET - 1 - j, true, Get_Label(bucket_pos + i));
					Set_Empty(bucket_pos + i);
					MoveElem(bucket_pos - j, bucket_pos + i);
				}
				else
				{
					assert(i == 0);
					UpdateFlag(bucket_pos, NUM_ELEMS_BUCKET - 1, true, Get_Label(bucket_pos));
				}
			}
	}
	uint32_t Find_Empty_Pos_Hopscotch(uint32_t bucket_pos, uint32_t bucket_init) noexcept
	{
		constexpr bool USE_REVERSE_WIN = true;
		uint32_t empty_pos = UINT32_MAX;

		//////////////////////////////////////////////////////////////////
		// Then try to reverse the bucket
		//////////////////////////////////////////////////////////////////
		if (USE_REVERSE_WIN) {
			if (!Is_Reversed_Window(bucket_pos) && bucket_pos >= NUM_ELEMS_BUCKET)
			{
				uint32_t count_empty = Count_Empty(bucket_pos + 1 - NUM_ELEMS_BUCKET);
				if (count_empty)
				{
					uint32_t count_elems = Count_Elems_In_Bucket_Non_Reversed(bucket_pos);

					// If we can reverse
					if (count_empty > count_elems || (count_empty == count_elems && Belong_to_Bucket(bucket_pos) == bucket_pos))
					{
						if(count_elems)// Some elems
							Reverse_Bucket(bucket_pos);
						else// No elem
							Set_Reversed(bucket_pos);

							uint32_t min1, pos1;
							bucket_init = bucket_pos + (Is_Reversed_Window(bucket_pos) ? (1 - NUM_ELEMS_BUCKET) : 0);
							std::tie(min1, pos1) = Calculate_Minimum(bucket_init);
							assert(min1 == 0);
							return pos1;
					}
				}
			}

			//////////////////////////////////////////////////////////////////
			// Then try to reverse elems
			//////////////////////////////////////////////////////////////////
			if (bucket_init >= 2 * NUM_ELEMS_BUCKET)
				for (uint32_t i = 0; i < NUM_ELEMS_BUCKET; i++)
				{
					uint32_t pos_elem = bucket_init + i;
					if (!Is_Item_In_Reverse_Bucket(pos_elem))
					{
						uint32_t bucket_elem = pos_elem - GetFlagDistance(pos_elem);

						if (bucket_elem != bucket_pos)
						{
							uint32_t count_empty = Count_Empty(bucket_elem + 1 - NUM_ELEMS_BUCKET);
							if (count_empty)
							{
								size_t count_elems, count_outside;
								std::tie(count_elems, count_outside) = Count_Elems_In_Bucket_Outside_Range(bucket_elem, bucket_init);

								assert(count_elems > count_outside);
								assert(count_elems >= 1);

								// TODO: Check this when only one element
								if (count_outside < count_empty && (count_empty >= count_elems || (count_empty+1 == count_elems && Belong_to_Bucket(bucket_elem) == bucket_elem)))
								{
									Reverse_Bucket(bucket_elem);

									uint32_t min1, pos1;
									std::tie(min1, pos1) = Calculate_Minimum(bucket_init);
									assert(min1 == 0);
									return pos1;
								}
							}
						}
					}
				}
		}
		//////////////////////////////////////////////////////////////////
		// Then try to hopscotch for an empty space
		//////////////////////////////////////////////////////////////////
		uint32_t max_dist_to_move = NUM_ELEMS_BUCKET - 1;
		for (uint32_t i = 0; i <= max_dist_to_move && (bucket_init + i) < num_buckets; i++)
		{
			if (Is_Empty(bucket_init + i))
			{
				// Find element to move
				uint32_t pos_blank = bucket_init + i;
				while ((pos_blank - bucket_init) >= NUM_ELEMS_BUCKET)
				{
					uint32_t pos_swap = pos_blank + 1 - NUM_ELEMS_BUCKET;

					for (; (pos_blank - pos_swap) > (NUM_ELEMS_BUCKET - 1 - GetFlagDistance(pos_swap)); pos_swap++)
					{
					}// TODO: Use a list with the options to not recalculate again

					 // Swap elements
					data1[pos_blank] = data1[pos_swap];
					Copy_Elem(pos_blank, pos_swap);
					UpdateFlag(pos_blank, GetFlagDistance(pos_swap) + (pos_blank - pos_swap), Is_Item_In_Reverse_Bucket(pos_swap));

					pos_blank = pos_swap;
				}

				return pos_blank;
			}
			uint32_t current_max_move = i + NUM_ELEMS_BUCKET - 1 - GetFlagDistance(bucket_init + i);
			if (current_max_move > max_dist_to_move)
				max_dist_to_move = current_max_move;
		}

		return empty_pos;
	}
	bool Is_Primary(uint32_t bucket_pos)
	{
		uint64_t hash = hash_elem(data1[bucket_pos]);

		// Calculate positions given hash
		uint32_t bucket1_pos = fastrange32((uint32_t)hash, num_buckets);

		return bucket1_pos == Belong_to_Bucket(bucket_pos);
	}

public:
	CBG_Set(uint32_t expected_num_elems) noexcept : num_elems(0), num_buckets(expected_num_elems)
	{
		data1 = (uint64_t*)malloc(num_buckets * sizeof(uint64_t));
		cache1 = (uint32_t*)malloc(num_buckets * sizeof(uint32_t));
		memset(cache1, 0, num_buckets * sizeof(uint32_t));

		for (uint32_t i = 0; i < (NUM_ELEMS_BUCKET-1); i++)
			Set_Reversed(num_buckets - 1 - i);
	}
	~CBG_Set()
	{
		if (data1) free(data1);
		if (cache1) free(cache1);

		data1 = nullptr;
		cache1 = nullptr;
		num_elems = 0;
		num_buckets = 0;
	}
	uint32_t capacity() const noexcept
	{
		return num_buckets;
	}
	uint32_t size() const noexcept
	{
		return num_elems;
	}
	void clear() noexcept
	{
		num_elems = 0;
		memset(cache1, 0, num_buckets * sizeof(uint32_t));

		for (uint32_t i = 0; i < (NUM_ELEMS_BUCKET - 1); i++)
			Set_Reversed(num_buckets - 1 - i);
	}
	double load_factor() const noexcept
	{
		return size() * 100. / capacity();
	}

	// Statistics
	std::tuple<uint32_t, uint32_t, uint32_t, uint32_t> Calculate_Statistics()
	{
		uint32_t sum_labels = 0;
		uint32_t lucky_buckets = 0;
		uint32_t num_elem_primary = 0;
		uint32_t num_reverse_bucket = 0;

		for (uint32_t i = 0; i < num_buckets; i++)
		{
			sum_labels += Get_Label(i);
			if (!Is_Unlucky_Bucket(i))
				lucky_buckets++;

			if (!Is_Empty(i) && Is_Primary(i))
				num_elem_primary++;

			if (Is_Reversed_Window(i))
				num_reverse_bucket++;
		}

		return std::make_tuple(sum_labels, lucky_buckets, num_elem_primary, num_reverse_bucket);
	}
	void Fill_Max_Distance(uint32_t* distances)
	{
		for (uint32_t i = 0; i < num_buckets; i++)
			if (!Is_Empty(i))
				distances[Is_Item_In_Reverse_Bucket(i) ? (NUM_ELEMS_BUCKET - 1 - GetFlagDistance(i)) : GetFlagDistance(i)]++;
	}
	
	bool insert(uint64_t elem, uint32_t L_MAX = 7) noexcept
	{
		while (true)
		{
			uint64_t hash = hash_elem(elem);

			// Calculate positions given hash
			uint32_t bucket1_pos = fastrange32((uint32_t)hash, num_buckets);
			uint32_t bucket2_pos = fastrange32(hash >> 32, num_buckets);

			bool is_reversed_bucket1 = Is_Reversed_Window(bucket1_pos);
			bool is_reversed_bucket2 = Is_Reversed_Window(bucket2_pos);
			uint32_t bucket1_init = bucket1_pos + (is_reversed_bucket1 ? (1 - NUM_ELEMS_BUCKET) : 0);
			uint32_t bucket2_init = bucket2_pos + (is_reversed_bucket2 ? (1 - NUM_ELEMS_BUCKET) : 0);

			// Find minimun label
			uint32_t min1 = Get_Label(bucket1_init);
			uint32_t min2 = Get_Label(bucket2_init);
			uint32_t pos1 = bucket1_init;
			uint32_t pos2 = bucket2_init;
			for (uint32_t i = 1; i < NUM_ELEMS_BUCKET /*&& (min1 || min2)*/; i++)
			{
				uint32_t current_pos1 = bucket1_init + i;
				uint32_t current_pos2 = bucket2_init + i;
				uint32_t label_value1 = Get_Label(current_pos1);
				uint32_t label_value2 = Get_Label(current_pos2);
				if (min1 > label_value1)
				{
					min1 = label_value1;
					pos1 = current_pos1;
				}
				if (min2 > label_value2)
				{
					min2 = label_value2;
					pos2 = current_pos2;
				}
			}

			//////////////////////////////////////////////////////////////////
			// No secondary added, no unlucky bucket added
			//////////////////////////////////////////////////////////////////
			// First bucket had free space
			if (min1 == 0)
			{
				UpdateFlag(pos1, pos1 - bucket1_init, is_reversed_bucket1, std::min(min2 + 1, L_MAX));
				// Put elem
				data1[pos1] = elem;
				num_elems++;
				return true;
			}

			uint32_t empty_pos = Find_Empty_Pos_Hopscotch(bucket1_pos, bucket1_init);
			if (empty_pos != UINT32_MAX)
			{
				is_reversed_bucket1 = Is_Reversed_Window(bucket1_pos);
				bucket1_init = bucket1_pos + (is_reversed_bucket1 ? (1 - NUM_ELEMS_BUCKET) : 0);
				UpdateFlag(empty_pos, empty_pos - bucket1_init, is_reversed_bucket1, std::min(min2 + 1, L_MAX));

				// Put elem
				data1[empty_pos] = elem;
				num_elems++;
				return true;
			}

			///////////////////////////////////////////////////////////////////
			// Secondary added, Unlucky bucket added
			//////////////////////////////////////////////////////////////////
			if (min2 == 0)
			{
				Set_Unlucky_Bucket(bucket1_pos);
				UpdateFlag(pos2, pos2 - bucket2_init, is_reversed_bucket2, std::min(min1 + 1, L_MAX));
				// Put elem
				data1[pos2] = elem;
				num_elems++;
				return true;
			}

			// TODO: Consider this
			//if (num_elems * 10 > 9 * num_buckets)// > 90%
			{
				empty_pos = Find_Empty_Pos_Hopscotch(bucket2_pos, bucket2_init);

				if (empty_pos != UINT32_MAX)
				{
					Set_Unlucky_Bucket(bucket1_pos);
					is_reversed_bucket2 = Is_Reversed_Window(bucket2_pos);
					bucket2_init = bucket2_pos + (is_reversed_bucket2 ? (1 - NUM_ELEMS_BUCKET) : 0);
					UpdateFlag(empty_pos, empty_pos - bucket2_init, is_reversed_bucket2, std::min(min1 + 1, L_MAX));

					// Put elem
					data1[empty_pos] = elem;
					num_elems++;
					return true;
				}
			}

			// Terminating condition
			if (std::min(min1, min2) >= L_MAX)
				return false;

			if (min1 <= min2)// Selected pos in first bucket
			{
				UpdateFlag(pos1, pos1 - bucket1_init, is_reversed_bucket1, std::min(min2 + 1, L_MAX));
				// Put elem
				uint64_t victim = data1[pos1];
				data1[pos1] = elem;
				elem = victim;
			}
			else
			{
				Set_Unlucky_Bucket(bucket1_pos);
				UpdateFlag(pos2, pos2 - bucket2_init, is_reversed_bucket2, std::min(min1 + 1, L_MAX));
				// Put elem
				uint64_t victim = data1[pos2];
				data1[pos2] = elem;
				elem = victim;
			}
		}
	}
	///////////////////////////////////////////////////////////////////////////////
	// Check if an element exist
	///////////////////////////////////////////////////////////////////////////////
	uint32_t count(uint64_t elem) const noexcept
	{
		uint64_t hash = hash_elem(elem);

		// Check first bucket
		uint32_t pos = fastrange32((uint32_t)hash, num_buckets);

		if (cmp_elems(data1[pos], elem) && !Is_Empty(pos))
			return 1;

		bool is_unlucky = Is_Unlucky_Bucket(pos);
		uint32_t reverse_sum = Is_Reversed_Window(pos) ? -1 : 1;
		pos += reverse_sum;

		for (uint32_t i = 0; i < (NUM_ELEMS_BUCKET-1); i++, pos += reverse_sum)
			if (cmp_elems(data1[pos], elem) && !Is_Empty(pos))
				return 1;
		
		// Check second bucket
		if (is_unlucky)
		{
			pos = fastrange32(hash >> 32, num_buckets);

			if (cmp_elems(data1[pos], elem) && !Is_Empty(pos))
				return 1;

			reverse_sum = Is_Reversed_Window(pos) ? -1 : 1;
			pos += reverse_sum;

			for (uint32_t i = 0; i < (NUM_ELEMS_BUCKET - 1); i++, pos += reverse_sum)
				if (cmp_elems(data1[pos], elem) && !Is_Empty(pos))
					return 1;
		}

		return 0;
	}

	uint32_t cache_lines_insert(uint32_t* cache_lines, uint32_t new_cache_line, uint32_t cache_line_counter) const noexcept
	{
		for (uint32_t i = 0; i < cache_line_counter; i++)
			if (cache_lines[i] == new_cache_line)
				return cache_line_counter;
		
		cache_lines[cache_line_counter] = new_cache_line;
		return cache_line_counter + 1;
	}
	size_t count_cache_lines(uint64_t elem, uint32_t elem_size) const noexcept
	{
		constexpr uint32_t CACHE_LINE_SIZE = 64;

		uint64_t hash = hash_elem(elem);

		// Check first bucket
		uint32_t pos = fastrange32((uint32_t)hash, num_buckets);

		uint32_t cache_lines[16];
		uint32_t cache_line_counter = 1;
		cache_lines[0] = pos * elem_size / CACHE_LINE_SIZE;

		cache_line_counter = cache_lines_insert(cache_lines, (pos * elem_size + elem_size - 1) / CACHE_LINE_SIZE, cache_line_counter);
		if (cmp_elems(data1[pos], elem) && !Is_Empty(pos))
			return cache_line_counter;

		bool is_unlucky = Is_Unlucky_Bucket(pos);
		uint32_t reverse_sum = Is_Reversed_Window(pos) ? -1 : 1;
		pos += reverse_sum;

		for (uint32_t i = 0; i < (NUM_ELEMS_BUCKET - 1); i++, pos += reverse_sum)
		{
			cache_line_counter = cache_lines_insert(cache_lines, pos * elem_size / CACHE_LINE_SIZE, cache_line_counter);
			cache_line_counter = cache_lines_insert(cache_lines, (pos * elem_size + elem_size - 1) / CACHE_LINE_SIZE, cache_line_counter);
			if (cmp_elems(data1[pos], elem) && !Is_Empty(pos))
				return cache_line_counter;
		}

		// Check second bucket
		if (is_unlucky)
		{
			pos = fastrange32(hash >> 32, num_buckets);

			cache_line_counter = cache_lines_insert(cache_lines, pos * elem_size / CACHE_LINE_SIZE, cache_line_counter);
			cache_line_counter = cache_lines_insert(cache_lines, (pos * elem_size + elem_size - 1) / CACHE_LINE_SIZE, cache_line_counter);
			if (cmp_elems(data1[pos], elem) && !Is_Empty(pos))
				return cache_line_counter;

			reverse_sum = Is_Reversed_Window(pos) ? -1 : 1;
			pos += reverse_sum;

			for (uint32_t i = 0; i < (NUM_ELEMS_BUCKET - 1); i++, pos += reverse_sum)
			{
				cache_line_counter = cache_lines_insert(cache_lines, pos * elem_size / CACHE_LINE_SIZE, cache_line_counter);
				cache_line_counter = cache_lines_insert(cache_lines, (pos * elem_size + elem_size - 1) / CACHE_LINE_SIZE, cache_line_counter);
				if (cmp_elems(data1[pos], elem) && !Is_Empty(pos))
					return cache_line_counter;
			}
		}

		printf("Error in lookup\n");
		return 0;
	}
};

///////////////////////////////////////////////////////////////////////////////
// Experiments
///////////////////////////////////////////////////////////////////////////////
void test_hashset()
{
	const uint32_t MAX_BUCKETS = 100'000;
	const uint32_t MAX_ELEMS = 100'000;

	CBG_Set<4> cuckoo_table(MAX_BUCKETS);

	// Fill the table
	std::mt19937_64 r;

	for (uint32_t i = 0; i < MAX_ELEMS; i++)
		if (!cuckoo_table.insert(r()))
		{
			printf("Error inserting at: %f%%\n", cuckoo_table.load_factor());
			break;
		}

	// Negative lookup
	std::mt19937_64 nr(45667);
	uint32_t num_found = 0;

	for (uint32_t i = 0; i < MAX_ELEMS; i++)
		if (cuckoo_table.count(nr()))
			num_found++;

	if (num_found)
		printf("This is probably an error on negative lookup\n");

	// Positive lookup
	num_found = 0;
	std::mt19937_64 r_positive;

	for (uint32_t i = 0; i <= cuckoo_table.size(); i++)
	{
		if (cuckoo_table.count(r_positive()))
			num_found++;
	}

	if (num_found != cuckoo_table.size())
		printf("Error in positive lookup\n");
}

#define CUCKOO_CAPACITY 100'000
#define MAX_REPETITIONS 1'000
#define MAX_L 15
struct TestData {
	double min;
	double max;
	double sum;
	uint32_t sum_labels;
	uint32_t elems_1st_bucket;
	uint32_t lucky_buckets;
	uint32_t reverse_buckets;
};
void test_lmax()
{
	TestData data[MAX_L];

	// Initialize data
	for (int i = 0; i < MAX_L; i++)
	{
		data[i].sum = 0;
		data[i].max = 0;
		data[i].min = 100;
		data[i].sum_labels = 0;
		data[i].elems_1st_bucket = 0;
		data[i].lucky_buckets = 0;
		data[i].reverse_buckets = 0;
	}

	// Seed with a real random value, if available
	std::random_device good_random;
	// Change first template parameter from 4 to other k to see other schemes
	CBG_Set<4> cuckoo_table(CUCKOO_CAPACITY);

	auto start = std::chrono::high_resolution_clock::now();

	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS; repeat++)
	{
		uint32_t l_max = 0;
		std::mt19937_64 r(good_random());
		cuckoo_table.clear();

		// Fill the table
		while (true)
		{
			if (!cuckoo_table.insert(r(), l_max + 1))
			{
				double rate = cuckoo_table.load_factor();
				data[l_max].sum += rate;

				uint32_t sum_labels, lucky_buckets, elems_1st_bucket, reverse_bucket;
				std::tie(sum_labels, lucky_buckets, elems_1st_bucket, reverse_bucket) = cuckoo_table.Calculate_Statistics();
				data[l_max].elems_1st_bucket += elems_1st_bucket;
				data[l_max].lucky_buckets += lucky_buckets;
				data[l_max].sum_labels += sum_labels;
				data[l_max].reverse_buckets += reverse_bucket;

				if (data[l_max].min > rate)
					data[l_max].min = rate;
				if (data[l_max].max < rate)
					data[l_max].max = rate;

				l_max++;
				if (l_max >= MAX_L)
					break;
			}
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("Time: %ums\n\n", (uint32_t)elapsed.count());

	printf("--------------------------------------------------------------------------\n");
	printf("L_max          Table_Use           Labels     Primaries   Lucky    Reverse\n");
	printf("     (min_diff-  avg   +max_diff)               Elem     Buckets   Buckets\n");
	printf("--------------------------------------------------------------------------\n");
	for (uint32_t l_max = 0; l_max < MAX_L; l_max++)
	{
		double avg_use = data[l_max].sum / MAX_REPETITIONS;
		printf("%2u     %5.2f%% - %5.2f%% + %5.2f%%    %5.2f        %3.0f%%      %3.0f%%       %.1f%%\n"
			, l_max + 1
			, avg_use - data[l_max].min, avg_use, data[l_max].max - avg_use
			, data[l_max].sum_labels*1. / MAX_REPETITIONS / CUCKOO_CAPACITY
			, data[l_max].elems_1st_bucket*100. / MAX_REPETITIONS / (CUCKOO_CAPACITY*avg_use/100)
			, data[l_max].lucky_buckets*100. / MAX_REPETITIONS / CUCKOO_CAPACITY
			, data[l_max].reverse_buckets*100. / MAX_REPETITIONS / CUCKOO_CAPACITY);
	}
}

void test_table_use()
{
	const int TO_SHOW = 18;

	uint32_t firsts[TO_SHOW];
	memset(firsts, 0, sizeof(firsts));
	uint32_t num_lucky_buckets[TO_SHOW];
	memset(num_lucky_buckets, 0, sizeof(num_lucky_buckets));
	uint32_t num_reverse_buckets[TO_SHOW];
	memset(num_reverse_buckets, 0, sizeof(num_reverse_buckets));
	uint32_t elem_dist[TO_SHOW][8];
	memset(elem_dist, 0, sizeof(elem_dist));

	// Seed with a real random value, if available
	std::random_device good_random;
	CBG_Set<4> cuckoo_table(100'000);

	auto start = std::chrono::high_resolution_clock::now();

	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS; repeat++)
	{
		std::mt19937_64 r(good_random());
		cuckoo_table.clear();

		// Initial filling
		for (int i = 0; i < 48'000; i++)
			if (!cuckoo_table.insert(r()))
				break;

		// Fill the table
		int index = 0;
		bool continue_cycle = true;

		while (continue_cycle)
		{
			for (int i = 0; i < 3'000; i++)
				if (!cuckoo_table.insert(r()))
				{
					continue_cycle = false;
					break;
				}

			uint32_t lucky_buckets, elems_1st_bucket, reverse_bucket;
			std::tie(std::ignore, lucky_buckets, elems_1st_bucket, reverse_bucket) = cuckoo_table.Calculate_Statistics();
			firsts[index] += elems_1st_bucket;
			num_lucky_buckets[index] += lucky_buckets;
			num_reverse_buckets[index] += reverse_bucket;

			cuckoo_table.Fill_Max_Distance(elem_dist[index]);

			index++;
			if (index == TO_SHOW - 1)
				break;
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("Time: %ums\n\n", (uint32_t)elapsed.count());

	printf("-----------------------------------------------------------------------------------------------\n");
	printf("Table_Use   First_Bucket   Lucky_Buckets   Rev_Buckets   Elems[0]  Elems[1]  Elems[2]  Elems[3]\n");
	printf("-----------------------------------------------------------------------------------------------\n");

	for (uint32_t i = 0; i < TO_SHOW - 1; i++)
	{
		uint32_t percent = 51 + i * 3;
		printf("  %3u%%          %.1f%%           %.1f%%         %.1f%%    "
			, percent, firsts[i]* 100. / MAX_REPETITIONS / (CUCKOO_CAPACITY*percent/100.)
			, num_lucky_buckets[i] * 100. / MAX_REPETITIONS / CUCKOO_CAPACITY
			, num_reverse_buckets[i] * 100. / MAX_REPETITIONS / CUCKOO_CAPACITY);

		for (uint32_t j = 0; j < 4; j++)
			printf("    %5.2f%%", elem_dist[i][j] * 100. / MAX_REPETITIONS / (CUCKOO_CAPACITY*percent/100.));
		printf("\n");
	}
	printf("\n");
}

#define MAX_REPETITIONS_ERRORS 1'000
void test_error()
{
	double rate[MAX_REPETITIONS_ERRORS];
	double avg_rate = 0;

	// Seed with a real random value, if available
	std::random_device good_random;
	// Change first template parameter from 4 to other k to see other schemes
	CBG_Set<2> cuckoo_table(CUCKOO_CAPACITY);

	auto start = std::chrono::high_resolution_clock::now();

	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS_ERRORS; repeat++)
	{
		std::mt19937_64 r(good_random());
		cuckoo_table.clear();

		// Fill the table
		while (true)
		{
			if (!cuckoo_table.insert(r(), 7))
			{
				rate[repeat] = cuckoo_table.load_factor();
				avg_rate += rate[repeat];
				break;
			}
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("Time: %ums\n\n", (uint32_t)elapsed.count());

	avg_rate /= MAX_REPETITIONS_ERRORS;
	double points[] = { 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 200 };
	uint32_t counters[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	// Calculate differences
	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS_ERRORS; repeat++)
	{
		double diff = avg_rate - rate[repeat];

		for (int i = 0; i < (sizeof(points) / sizeof(points[0])); i++)
			if (diff < points[i])
			{
				counters[i]++;
				break;
			}
	}

	// Show the histogram
	for (int i = 0; i < (sizeof(points) / sizeof(points[0])); i++)
		printf("%03.1f  %u\n", points[i], counters[i]);
}

void test_offline()
{
	constexpr uint32_t MAX_WINDOW_SIZE = 8;
	// Seed with a real random value, if available
	std::random_device good_random;

	std::unique_ptr<uint32_t[]> table = std::make_unique<uint32_t[]>(CUCKOO_CAPACITY);
	uint32_t windows[MAX_WINDOW_SIZE];
	memset(windows, 0, sizeof(windows));

	auto start = std::chrono::high_resolution_clock::now();

	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS; repeat++)
	{
		std::mt19937 r(good_random());
		memset(table.get(), 0, CUCKOO_CAPACITY * sizeof(uint32_t));

		// Fill the table
		for (uint32_t i = 0; i < CUCKOO_CAPACITY; i++)
			table[r() % CUCKOO_CAPACITY]++;

		for (uint32_t i = 0; i < CUCKOO_CAPACITY; i++)
			if(table[i] < MAX_WINDOW_SIZE)
				windows[table[i]]++;
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("Time: %ums\n\n", (uint32_t)elapsed.count());

	printf("-----------------------------\n");
	printf("#_Elems  Bins   Lucky_Buckets\n");
	printf("-----------------------------\n");
	double lucky_buckets = 0;
	for (uint32_t i = 0; i < MAX_WINDOW_SIZE; i++)
	{
		lucky_buckets += windows[i] * 100. / MAX_REPETITIONS / CUCKOO_CAPACITY;
		printf("%u       %5.2f%%     %5.2f%%\n", i, windows[i] * 100. / MAX_REPETITIONS / CUCKOO_CAPACITY, lucky_buckets);
	}
}

void test_insertion_time()
{
	uint64_t times[100];
	memset(times, 0, sizeof(times));

	// Seed with a real random value, if available
	std::random_device good_random;
	CBG_Set<2> cuckoo_table(100'000);

	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS; repeat++)
	{
		std::mt19937_64 r(good_random());
		cuckoo_table.clear();

		// Fill the table
		int per_cent = 0;
		bool continue_cycle = true;
		while (continue_cycle)
		{
			auto start = std::chrono::high_resolution_clock::now();

			for (int i = 0; i < 1000; i++)
				if (!cuckoo_table.insert(r()))
				{
					continue_cycle = false;
					break;
				}

			auto end = std::chrono::high_resolution_clock::now();
			auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
			times[per_cent] += elapsed.count();
			per_cent++;
		}
	}

	for (uint32_t i = 0; i < 100; i++)
		printf("%02u%% time: %I64u ns\n", i, times[i] / MAX_REPETITIONS);
}

void test_cache_lines()
{
	constexpr uint32_t NUM_ELEMENT_SIZE = 8;
	uint32_t total_cache_lines[4][NUM_ELEMENT_SIZE];
	uint32_t element_size[NUM_ELEMENT_SIZE] = {64, 32, 21, 16, 12, 10, 9, 8};
	uint32_t num_elements[4] = { 50'000, 75'000, 90'000, 99'000 };

	memset(total_cache_lines, 0, sizeof(total_cache_lines));

	// Seed with a real random value, if available
	std::random_device good_random;
	CBG_Set<6> cuckoo_table(100'000);

	auto start = std::chrono::high_resolution_clock::now();

	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS; repeat++)
	{
		uint32_t seed = good_random();
		std::mt19937_64 r(seed);
		cuckoo_table.clear();

		// Initial filling
		for (uint32_t table_load = 0; table_load < 4; table_load++)
		{
			while (cuckoo_table.size() < num_elements[table_load])
				if (!cuckoo_table.insert(r()))
					printf("Error Inserting/n");

			std::mt19937_64 r_positive(seed);
			for (uint32_t i = 0; i < num_elements[table_load]; i++)
			{
				uint64_t elem = r_positive();

				for (uint32_t j = 0; j < NUM_ELEMENT_SIZE; j++)
					total_cache_lines[table_load][j] += (uint32_t)cuckoo_table.count_cache_lines(elem, element_size[j]);
			}
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("Time: %ums\n\n", (uint32_t)elapsed.count());

	printf("--------------------------------\n");
	printf("Table_Use   Elem_Size   Positive\n");
	printf("--------------------------------\n");

	for (uint32_t i = 0; i < 4; i++)
	{
		for (uint32_t j = 0; j < NUM_ELEMENT_SIZE; j++)
		{
			uint32_t percent = 51 + i * 3;
			printf("  %u%%          %2u        %.2f\n"
				, num_elements[i] / 1000
				, element_size[j]
				, total_cache_lines[i][j] * 1. / MAX_REPETITIONS / num_elements[i]);
		}
		printf("--------------------------------\n");
	}
	printf("\n");
}


void main()
{
	//test_hashset();
	//test_lmax();
	//test_table_use();
	//test_error();
	//test_offline();
	//test_insertion_time();
	test_cache_lines();

	// Wait for one keystroke
	printf("\nPress any key to exit...");
	char c;
	scanf("%c", &c);
}