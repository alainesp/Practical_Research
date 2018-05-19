///////////////////////////////////////////////////////////////////////////////
// Research code on Cuckoo Hashing Insertion (LSA_max)
///////////////////////////////////////////////////////////////////////////////
//
// Written by Alain Espinosa <alainesp at gmail.com> in 2018 and placed
// under the MIT license (see LICENSE file for a full definition).
//
///////////////////////////////////////////////////////////////////////////////

#include <stdint.h>
#include <algorithm>
#include <numeric>
#include <random>

///////////////////////////////////////////////////////////////////////////////
// Cucko hashing insertion code (Random Walk and LSA_max) implemented as only
// one table. Contains implementations for d=[2,3] and k=<any>
//
// This is for research only, see 'Cuckoo_4x2_1t_lookup' below
// for a proof-of-concept production code.
///////////////////////////////////////////////////////////////////////////////
template<uint32_t NUM_ELEMS_BUCKET, typename T, bool IS_BLOCKED> class Cuckoo_1t
{
private:
	T* data;
	uint8_t* labels;
	uint32_t num_buckets;
	uint32_t num_elems;
	std::mt19937_64 r;

public:
	// Statistics
	uint32_t same_bucket_count;
	uint32_t num_moves;
	uint32_t sum_labels;

	Cuckoo_1t(uint32_t expected_num_elems, uint32_t random_seed) noexcept : r(random_seed), num_elems(0), same_bucket_count(0), num_moves(0), sum_labels(0)
	{
		if (IS_BLOCKED)
		{
			num_buckets = (expected_num_elems + NUM_ELEMS_BUCKET - 1) / NUM_ELEMS_BUCKET;
			data = (T*)malloc(num_buckets * NUM_ELEMS_BUCKET * sizeof(T));

			labels = (uint8_t*)malloc(num_buckets * NUM_ELEMS_BUCKET * sizeof(uint8_t));
			memset(labels, 0, num_buckets * NUM_ELEMS_BUCKET * sizeof(uint8_t));
		}
		else
		{
			num_buckets = expected_num_elems;
			data = (T*)malloc((num_buckets + NUM_ELEMS_BUCKET - 1) * sizeof(T));

			labels = (uint8_t*)malloc((num_buckets + NUM_ELEMS_BUCKET - 1) * sizeof(uint8_t));
			memset(labels, 0, (num_buckets + NUM_ELEMS_BUCKET - 1) * sizeof(uint8_t));
		}
	}
	~Cuckoo_1t() noexcept
	{
		free(data);
		free(labels);
	}
	uint32_t GetCapacity() const noexcept
	{
		return IS_BLOCKED ? NUM_ELEMS_BUCKET * num_buckets : (num_buckets + NUM_ELEMS_BUCKET - 1);
	}
	uint32_t GetNumElems() const noexcept
	{
		return num_elems;
	}
	double GetTableUse() const noexcept
	{
		return 100. * num_elems / GetCapacity();
	}
	void Clear()
	{
		num_elems = 0;
		same_bucket_count = 0;
		num_moves = 0;
		sum_labels = 0;

		memset(labels, 0, (IS_BLOCKED ? (num_buckets * NUM_ELEMS_BUCKET) : (num_buckets + NUM_ELEMS_BUCKET - 1)) * sizeof(uint8_t));
	}
	inline uint64_t HashElem(T elem) const noexcept
	{
		return elem;
	}

	///////////////////////////////////////////////////////////////////////////////
	// Random Walk insertion
	///////////////////////////////////////////////////////////////////////////////
	bool Insert_RW_k2(T elem) noexcept
	{
		const uint32_t MAX_KICKS = 500;

		for (uint32_t level = 0, last_pos = UINT32_MAX; level < MAX_KICKS; level++)
		{
			uint64_t hash = HashElem(elem);

			// Calculate positions given hash
			uint64_t h1 = (uint32_t)hash;
			uint64_t h2 = (uint32_t)(hash >> 32);
			uint32_t pos1 = (h1 % num_buckets) * (IS_BLOCKED ? NUM_ELEMS_BUCKET : 1);
			uint32_t pos2 = (h2 % num_buckets) * (IS_BLOCKED ? NUM_ELEMS_BUCKET : 1);

			// Try to insert the element into a not full bucket. Note buckets fills evenly
			for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
			{
				if (labels[pos1 + i] == 0)
				{
					labels[pos1 + i] = 1;
					data[pos1 + i] = elem;
					num_elems++;
					return true;
				}
				if (labels[pos2 + i] == 0)
				{
					labels[pos2 + i] = 1;
					data[pos2 + i] = elem;
					num_elems++;
					return true;
				}
			}

			// Break infinity cicle
			if (pos1 == pos2 && last_pos==pos1 && NUM_ELEMS_BUCKET == 1)
				return false;

			// The two buckets are full
			uint32_t victim_pos;
			do
			{
				int random_chooice = r() % (2 * NUM_ELEMS_BUCKET);
				victim_pos = (random_chooice >= NUM_ELEMS_BUCKET ? pos2 : pos1) + (random_chooice % NUM_ELEMS_BUCKET);

			} while (victim_pos == last_pos);

			T victim = data[victim_pos];
			data[victim_pos] = elem;
			elem = victim;
			last_pos = victim_pos;
			num_moves++;// Not needed for algorithm to work, only statistics
		}

		return false;
	}
	bool Insert_RW_k3(T elem) noexcept
	{
		const uint32_t MAX_KICKS = 500;

		for (uint32_t level = 0, last_pos = UINT32_MAX; level < MAX_KICKS; level++)
		{
			uint64_t hash = HashElem(elem);

			// Calculate positions given hash
			uint64_t h1 = (uint32_t)hash;
			uint64_t h2 = (uint32_t)(hash >> 32);
			uint32_t pos1 = (h1 % num_buckets) * (IS_BLOCKED ? NUM_ELEMS_BUCKET : 1);
			uint32_t pos2 = ((h1 + h2) % num_buckets)  * (IS_BLOCKED ? NUM_ELEMS_BUCKET : 1);
			uint32_t pos3 = ((h1 + 2 * h2) % num_buckets) * (IS_BLOCKED ? NUM_ELEMS_BUCKET : 1);

			// Try to insert the element into a not full bucket. Note buckets fills evenly
			for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
			{
				if (labels[pos1 + i] == 0)
				{
					labels[pos1 + i] = 1;
					data[pos1 + i] = elem;
					num_elems++;
					return true;
				}
				if (labels[pos2 + i] == 0)
				{
					labels[pos2 + i] = 1;
					data[pos2 + i] = elem;
					num_elems++;
					return true;
				}
				if (labels[pos3 + i] == 0)
				{
					labels[pos3 + i] = 1;
					data[pos3 + i] = elem;
					num_elems++;
					return true;
				}
			}

			// Break infinity cicle
			if (pos1 == pos2 && pos1 == pos3 && last_pos == pos1 && NUM_ELEMS_BUCKET == 1)
				return false;

			// The three buckets are full
			uint32_t victim_pos;
			do
			{
				int random_chooice = r() % (3 * NUM_ELEMS_BUCKET);
				victim_pos = (random_chooice >= NUM_ELEMS_BUCKET ? (random_chooice >= 2 * NUM_ELEMS_BUCKET ? pos3 : pos2) : pos1) + (random_chooice % NUM_ELEMS_BUCKET);

			} while (victim_pos == last_pos);

			T victim = data[victim_pos];
			data[victim_pos] = elem;
			elem = victim;
			last_pos = victim_pos;
			num_moves++;// Not needed for algorithm to work, only statistics
		}

		return false;
	}

	///////////////////////////////////////////////////////////////////////////////
	// LSA_max insertion
	///////////////////////////////////////////////////////////////////////////////
	uint32_t Calculate_Load(uint32_t bucket_pos)
	{
		uint32_t load = 0;

		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
			load += labels[bucket_pos + i];

		return load;
	}
	uint32_t TB_LL_Global(uint32_t bucket1_pos, uint32_t bucket2_pos, uint32_t* min_label, uint32_t* min_label2)
	{
		uint32_t items_pos[2 * NUM_ELEMS_BUCKET];
		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			items_pos[i] = bucket1_pos + i;
			items_pos[i + NUM_ELEMS_BUCKET] = bucket2_pos + i;
		}

		// Find minimun label
		*min_label = UINT32_MAX;
		*min_label2 = UINT32_MAX;
		uint32_t pos = 0;

		for (int i = 0; i < 2 * NUM_ELEMS_BUCKET; i++)
		{
			if (labels[items_pos[i]] < *min_label)
			{
				*min_label2 = *min_label;
				*min_label = labels[items_pos[i]];
				pos = items_pos[i];
			}
			else if (labels[items_pos[i]] < *min_label2)
				*min_label2 = labels[items_pos[i]];
		}

		// If there are elements here
		if (*min_label)
		{
			uint32_t ll = UINT32_MAX;

			for (int i = 0; i < 2 * NUM_ELEMS_BUCKET; i++)
			{
				if (labels[items_pos[i]] == *min_label)
				{
					// Calculate alternate
					uint64_t hash = HashElem(data[items_pos[i]]);
					// Calculate positions given hash
					uint32_t h1 = (uint32_t)hash;
					uint32_t h2 = (uint32_t)(hash >> 32);
					uint32_t pos1 = (h1 % num_buckets) * NUM_ELEMS_BUCKET;
					uint32_t pos2 = (h2 % num_buckets) * NUM_ELEMS_BUCKET;
					uint32_t bucket_pos = (items_pos[i] / NUM_ELEMS_BUCKET) * NUM_ELEMS_BUCKET;
					uint32_t new_ll = Calculate_Load(bucket_pos == pos1 ? pos2 : pos1);

					if (new_ll < ll)
					{
						ll = new_ll;
						pos = items_pos[i];
					}
				}
			}
		}

		return pos;
	}
	uint32_t TB_Left_Global(uint32_t bucket1_pos, uint32_t bucket2_pos, uint32_t* min_label, uint32_t* min_label2)
	{
		uint32_t items_pos[2 * NUM_ELEMS_BUCKET];
		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			items_pos[i] = bucket1_pos + i;
			items_pos[i + NUM_ELEMS_BUCKET] = bucket2_pos + i;
		}

		// Find minimun label
		*min_label = UINT32_MAX;
		*min_label2 = UINT32_MAX;
		uint32_t pos = 0;

		for (int i = 0; i < 2 * NUM_ELEMS_BUCKET; i++)
		{
			if (labels[items_pos[i]] < *min_label)
			{
				*min_label2 = *min_label;
				*min_label = labels[items_pos[i]];
				pos = items_pos[i];
			}
			else if (labels[items_pos[i]] < *min_label2)
				*min_label2 = labels[items_pos[i]];
		}

		return pos;
	}
	uint32_t TB_Random_Global(uint32_t bucket1_pos, uint32_t bucket2_pos, uint32_t* min_label, uint32_t* min_label2)
	{
		uint32_t items_pos[2 * NUM_ELEMS_BUCKET];
		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			items_pos[i] = bucket1_pos + i;
			items_pos[i + NUM_ELEMS_BUCKET] = bucket2_pos + i;
		}

		// Find minimun label
		*min_label = UINT32_MAX;
		*min_label2 = UINT32_MAX;
		uint32_t pos = 0;

		for (int i = 0; i < 2 * NUM_ELEMS_BUCKET; i++)
		{
			if (labels[items_pos[i]] < *min_label)
			{
				*min_label2 = *min_label;
				*min_label = labels[items_pos[i]];
				pos = items_pos[i];
			}
			else if (labels[items_pos[i]] < *min_label2)
				*min_label2 = labels[items_pos[i]];
		}

		uint32_t count_min = 0;
		for (int i = 0; i < 2 * NUM_ELEMS_BUCKET; i++)
			if (labels[items_pos[i]] == *min_label)
				count_min++;

		uint32_t random_pos = r() % count_min;
		count_min = 0;
		for (int i = 0; i < 2 * NUM_ELEMS_BUCKET; i++)
			if (labels[items_pos[i]] == *min_label)
			{
				if (count_min == random_pos)
				{
					pos = items_pos[i];
					break;
				}
				count_min++;
			}

		return pos;
	}

	uint32_t TB_LL_LL(uint32_t bucket1_pos, uint32_t bucket2_pos, uint32_t* min_label, uint32_t* min_label2)
	{
		uint32_t items_pos[2 * NUM_ELEMS_BUCKET];
		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			items_pos[i] = bucket1_pos + i;
			items_pos[i + NUM_ELEMS_BUCKET] = bucket2_pos + i;
		}

		// Find minimun label
		*min_label = UINT32_MAX;
		*min_label2 = UINT32_MAX;
		uint32_t pos = 0;

		for (int i = 0; i < 2 * NUM_ELEMS_BUCKET; i++)
		{
			if (labels[items_pos[i]] < *min_label)
			{
				*min_label2 = *min_label;
				*min_label = labels[items_pos[i]];
				pos = items_pos[i];
			}
			else if (labels[items_pos[i]] < *min_label2)
				*min_label2 = labels[items_pos[i]];
		}

		uint32_t ll = UINT32_MAX;
		uint32_t bucket_displacement = Calculate_Load(bucket1_pos) <= Calculate_Load(bucket2_pos) ? 0 : NUM_ELEMS_BUCKET;

		if (*min_label)
			for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
				if (labels[items_pos[bucket_displacement + i]] == *min_label)
				{
					// Calculate alternate
					uint64_t hash = HashElem(data[items_pos[bucket_displacement + i]]);
					// Calculate positions given hash
					uint32_t h1 = (uint32_t)hash;
					uint32_t h2 = (uint32_t)(hash >> 32);
					uint32_t pos1 = (h1 % num_buckets) * NUM_ELEMS_BUCKET;
					uint32_t pos2 = (h2 % num_buckets) * NUM_ELEMS_BUCKET;
					uint32_t bucket_pos = (items_pos[bucket_displacement + i] / NUM_ELEMS_BUCKET) * NUM_ELEMS_BUCKET;
					uint32_t new_ll = Calculate_Load(bucket_pos == pos1 ? pos2 : pos1);

					if (new_ll < ll)
					{
						ll = new_ll;
						pos = items_pos[bucket_displacement + i];
					}
				}
				else// There are no elements
					for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
						if (labels[items_pos[bucket_displacement + i]] == 0)
						{
							pos = items_pos[bucket_displacement + i];
							break;
						}

		return pos;
	}
	uint32_t TB_LL_Left(uint32_t bucket1_pos, uint32_t bucket2_pos, uint32_t* min_label, uint32_t* min_label2)
	{
		uint32_t items_pos[2 * NUM_ELEMS_BUCKET];
		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			items_pos[i] = bucket1_pos + i;
			items_pos[i + NUM_ELEMS_BUCKET] = bucket2_pos + i;
		}

		// Find minimun label
		*min_label = UINT32_MAX;
		*min_label2 = UINT32_MAX;
		uint32_t pos = 0;

		for (int i = 0; i < 2 * NUM_ELEMS_BUCKET; i++)
		{
			if (labels[items_pos[i]] < *min_label)
			{
				*min_label2 = *min_label;
				*min_label = labels[items_pos[i]];
				pos = items_pos[i];
			}
			else if (labels[items_pos[i]] < *min_label2)
				*min_label2 = labels[items_pos[i]];
		}

		uint32_t bucket_displacement = Calculate_Load(bucket1_pos) <= Calculate_Load(bucket2_pos) ? 0 : NUM_ELEMS_BUCKET;

		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
			if (labels[items_pos[bucket_displacement + i]] == *min_label)
			{
				pos = items_pos[bucket_displacement + i];
				break;
			}

		return pos;
	}
	uint32_t TB_LL_Left_k3(uint32_t bucket1_pos, uint32_t bucket2_pos, uint32_t bucket3_pos, uint32_t* min_label, uint32_t* min_label2)
	{
		uint32_t items_pos[3 * NUM_ELEMS_BUCKET];
		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			items_pos[i] = bucket1_pos + i;
			items_pos[i + NUM_ELEMS_BUCKET] = bucket2_pos + i;
			items_pos[i + 2 * NUM_ELEMS_BUCKET] = bucket3_pos + i;
		}

		// Find minimun label
		*min_label = UINT32_MAX;
		*min_label2 = UINT32_MAX;
		uint32_t pos = 0;

		for (int i = 0; i < 3 * NUM_ELEMS_BUCKET; i++)
		{
			if (labels[items_pos[i]] < *min_label)
			{
				*min_label2 = *min_label;
				*min_label = labels[items_pos[i]];
				pos = items_pos[i];
			}
			else if (labels[items_pos[i]] < *min_label2)
				*min_label2 = labels[items_pos[i]];
		}

		uint32_t b1_load = Calculate_Load(bucket1_pos);
		uint32_t b2_load = Calculate_Load(bucket2_pos);
		uint32_t b3_load = Calculate_Load(bucket3_pos);

		uint32_t bucket_displacement = 0;

		if (b3_load < b1_load && b3_load <= b2_load)
			bucket_displacement = 2 * NUM_ELEMS_BUCKET;
		else if (b2_load < b1_load && b2_load <= b3_load)
			bucket_displacement = NUM_ELEMS_BUCKET;

		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
			if (labels[items_pos[bucket_displacement + i]] == *min_label)
			{
				pos = items_pos[bucket_displacement + i];
				break;
			}

		return pos;
	}
	uint32_t TB_LL_Random(uint32_t bucket1_pos, uint32_t bucket2_pos, uint32_t* min_label, uint32_t* min_label2)
	{
		uint32_t items_pos[2 * NUM_ELEMS_BUCKET];
		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			items_pos[i] = bucket1_pos + i;
			items_pos[i + NUM_ELEMS_BUCKET] = bucket2_pos + i;
		}

		// Find minimun label
		*min_label = UINT32_MAX;
		*min_label2 = UINT32_MAX;
		uint32_t pos = 0;

		for (int i = 0; i < 2 * NUM_ELEMS_BUCKET; i++)
		{
			if (labels[items_pos[i]] < *min_label)
			{
				*min_label2 = *min_label;
				*min_label = labels[items_pos[i]];
				pos = items_pos[i];
			}
			else if (labels[items_pos[i]] < *min_label2)
				*min_label2 = labels[items_pos[i]];
		}

		uint32_t bucket_displacement = Calculate_Load(bucket1_pos) <= Calculate_Load(bucket2_pos) ? 0 : NUM_ELEMS_BUCKET;

		uint32_t count_min = 0;
		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
			if (labels[items_pos[bucket_displacement + i]] == *min_label)
				count_min++;

		uint32_t random_pos = r() % count_min;
		count_min = 0;
		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
			if (labels[items_pos[bucket_displacement + i]] == *min_label)
			{
				if (count_min == random_pos)
				{
					pos = items_pos[bucket_displacement + i];
					break;
				}
				count_min++;
			}

		return pos;
	}

	bool Insert_LSA_max_k2(T elem, uint32_t L_max) noexcept
	{
		uint32_t last_bucket_pos = UINT32_MAX;

		while (true)
		{
			uint64_t hash = HashElem(elem);

			// Calculate positions given hash
			uint32_t h1 = (uint32_t)hash;
			uint32_t h2 = (uint32_t)(hash >> 32);
			uint32_t pos1 = (h1 % num_buckets) * (IS_BLOCKED ? NUM_ELEMS_BUCKET : 1);
			uint32_t pos2 = (h2 % num_buckets) * (IS_BLOCKED ? NUM_ELEMS_BUCKET : 1);

			// Find minimun label
			uint32_t min_label, min_label2;
			uint32_t pos = TB_LL_Left(pos1, pos2, &min_label, &min_label2);

			// Terminating condition
			if (min_label >= L_max)
				return false;

			sum_labels += min_label2 + 1 - labels[pos];// Not needed for algorithm to work, only statistics
			labels[pos] = min_label2 + 1;
			// Put elem
			if (min_label)
			{
				T victim = data[pos];
				data[pos] = elem;
				elem = victim;

				num_moves++;// Not needed for algorithm to work, only statistics

				// Count when an item reinsert into a bucket where it was kicked-out.
				// Not needed for algorithm to work, only statistics
				uint32_t bucket_pos = IS_BLOCKED ? (pos / NUM_ELEMS_BUCKET * NUM_ELEMS_BUCKET) : ((pos - pos1) < NUM_ELEMS_BUCKET ? pos1 : pos2);
				if (last_bucket_pos == bucket_pos)
					same_bucket_count++;
				last_bucket_pos = bucket_pos;
			}
			else// bin is empty
			{
				data[pos] = elem;
				num_elems++;
				return true;
			}
		}
	}
	bool Insert_LSA_max_k3(T elem, uint32_t L_max) noexcept
	{
		uint32_t last_bucket_pos = UINT32_MAX;

		while (true)
		{
			uint64_t hash = HashElem(elem);

			// Calculate positions given hash
			uint32_t h1 = (uint32_t)hash;
			uint32_t h2 = (uint32_t)(hash >> 32);
			uint32_t pos1 = (h1 % num_buckets) * (IS_BLOCKED ? NUM_ELEMS_BUCKET : 1);
			uint32_t pos2 = ((h1 + h2) % num_buckets) * (IS_BLOCKED ? NUM_ELEMS_BUCKET : 1);
			uint32_t pos3 = ((h1 + 2 * h2) % num_buckets) * (IS_BLOCKED ? NUM_ELEMS_BUCKET : 1);

			// Find minimun label
			uint32_t min_label, min_label2;
			uint32_t pos = TB_LL_Left_k3(pos1, pos2, pos3, &min_label, &min_label2);

			// Terminating condition
			if (min_label >= L_max)
				return false;

			sum_labels += min_label2 + 1 - labels[pos];// Not needed for algorithm to work, only statistics
			labels[pos] = min_label2 + 1;
			// Put elem
			if (min_label)
			{
				T victim = data[pos];
				data[pos] = elem;
				elem = victim;

				num_moves++;// Not needed for algorithm to work, only statistics

				// Count when an item reinsert into a bucket where it was kicked-out.
				// Not needed for algorithm to work, only statistics
				uint32_t bucket_pos = IS_BLOCKED ? (pos / NUM_ELEMS_BUCKET * NUM_ELEMS_BUCKET) : ((pos - pos1) < NUM_ELEMS_BUCKET ? pos1 : pos2);
				if (last_bucket_pos == bucket_pos)
					same_bucket_count++;
				last_bucket_pos = bucket_pos;
			}
			else// bin is empty
			{
				data[pos] = elem;
				num_elems++;
				return true;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	// Check if an element exist. Used for testing
	///////////////////////////////////////////////////////////////////////////////
	bool Exist_k2(T elem) const noexcept
	{
		uint64_t hash = HashElem(elem);

		// Calculate positions given hash
		uint32_t h1 = (uint32_t)hash;
		uint32_t h2 = (uint32_t)(hash >> 32);
		uint32_t pos1 = (h1 % num_buckets) * NUM_ELEMS_BUCKET;
		uint32_t pos2 = (h2 % num_buckets) * NUM_ELEMS_BUCKET;

		for (int i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			if (labels[pos1 + i] && data[pos1 + i] == elem)
				return true;
			if (labels[pos2 + i] && data[pos2 + i] == elem)
				return true;
		}

		return false;
	}
};

///////////////////////////////////////////////////////////////////////////////
// Proof-of-concept Cuckoo hashing insertion LSA_max implemented as only
// one table, with lookup table for scheme (2,4)
///////////////////////////////////////////////////////////////////////////////

// Lookup table for Cuckoo (2,4)-scheme used on insertion algorithm LSA_max
static uint16_t lt_2x4[] = {
	1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
	,257,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67,67
	,257,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3
	,257,323,259,135,261,135,135,135,265,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135,135
	,257,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5
	,257,323,259,71,261,71,71,71,265,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71,71
	,257,323,259,7,261,7,7,7,265,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7
	,257,323,259,391,261,327,263,208,265,331,267,208,269,208,208,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216,216
	,257,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9,9
	,257,323,259,75,261,75,75,75,265,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75,75
	,257,323,259,11,261,11,11,11,265,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11,11
	,257,323,259,391,261,327,263,144,265,331,267,144,269,144,144,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148,148
	,257,323,259,13,261,13,13,13,265,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13
	,257,323,259,391,261,327,263,80,265,331,267,80,269,80,80,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82,82
	,257,323,259,391,261,327,263,16,265,331,267,16,269,16,16,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17,17
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83,83
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19,19
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,151,277,151,151,151,281,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151,151
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,87,277,87,87,87,281,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87,87
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,23,277,23,23,23,281,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23,23
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,224,281,347,283,224,285,224,224,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232,232
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25,25
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,91,277,91,91,91,281,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91,91
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,27,277,27,27,27,281,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27,27
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,160,281,347,283,160,285,160,160,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164,164
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,29,277,29,29,29,281,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29,29
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,96,281,347,283,96,285,96,96,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98,98
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,32,281,347,283,32,285,32,32,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33,33
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99,99
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35,35
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,167,293,167,167,167,297,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167,167
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37,37
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,103,293,103,103,103,297,103,103,103,103,103,103,103,103,103,103,103,103,103,103,103,103,103,103,103,103,103,103,103
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,39,293,39,39,39,297,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39,39
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,240,297,363,299,240,301,240,240,248,248,248,248,248,248,248,248,248,248,248,248,248,248,248,248,248
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41,41
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,107,293,107,107,107,297,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107,107
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,43,293,43,43,43,297,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43,43
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,176,297,363,299,176,301,176,176,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180,180
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,45,293,45,45,45,297,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45,45
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,112,297,363,299,112,301,112,112,114,114,114,114,114,114,114,114,114,114,114,114,114,114,114,114,114
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,48,297,363,299,48,301,48,48,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49,49
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,115,115,115,115,115,115,115,115,115,115,115,115,115,115,115
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,51,51,51,51,51,51,51,51,51,51,51,51,51,51,51
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,183,309,183,183,183,313,183,183,183,183,183,183,183
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,53,53,53,53,53,53,53,53,53,53,53,53,53,53,53
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,119,309,119,119,119,313,119,119,119,119,119,119,119
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,55,309,55,55,55,313,55,55,55,55,55,55,55
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,439,309,375,311,255,313,379,315,255,317,255,255,255
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,57,57,57,57,57,57,57,57,57,57,57,57,57,57,57
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,123,309,123,123,123,313,123,123,123,123,123,123,123
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,59,309,59,59,59,313,59,59,59,59,59,59,59
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,439,309,375,311,191,313,379,315,191,317,191,191,191
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,61,309,61,61,61,313,61,61,61,61,61,61,61
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,439,309,375,311,127,313,379,315,127,317,127,127,127
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,439,309,375,311,63,313,379,315,63,317,63,63,63
	,257,323,259,391,261,327,263,472,265,331,267,404,269,338,273,273,273,339,275,407,277,343,279,488,281,347,283,420,285,354,289,289,289,355,291,423,293,359,295,504,297,363,299,436,301,370,305,305,305,371,307,439,309,375,311,511,313,379,315,447,317,383,319,575
};

template<typename T> class Cuckoo_2x4_1t_lookup
{
private:
	T* data;
	uint8_t* labels;		// Use an uint8_t for a bucket (for fast code): 2 bits per element
	uint64_t num_buckets;
	uint32_t num_elems;

public:
	Cuckoo_2x4_1t_lookup(uint32_t expected_num_elems) noexcept : num_buckets((expected_num_elems + 3) / 4), num_elems(0)
	{
		data = (T*)malloc(num_buckets * 4 * sizeof(T));

		labels = (uint8_t*)malloc(num_buckets * sizeof(uint8_t));
		memset(labels, 0, num_buckets * sizeof(uint8_t));
	}
	~Cuckoo_2x4_1t_lookup() noexcept
	{
		free(data);
		free(labels);
	}
	uint32_t GetCapacity() const noexcept
	{
		return 4 * div.divisor;
	}
	uint32_t GetNumElems() const noexcept
	{
		return num_elems;
	}
	double GetFillRate() const noexcept
	{
		return 100. * num_elems / GetCapacity();
	}
	inline uint64_t HashElem(T elem) const noexcept
	{
		return elem;
	}

	bool Insert(T elem) noexcept
	{
		while (true)
		{
			uint64_t hash = HashElem(elem);

			// Calculate positions given hash
			uint32_t pos1 = (((hash & 0xffffffff) * num_buckets) >> 32);
			uint32_t pos2 = (((hash >> 32) * num_buckets) >> 32);

			// Find minimun label and other data
			uint32_t result = lt_2x4[(((uint32_t)labels[pos1]) << 6) | labels[pos2]];

			// Terminating condition
			if (result & 0x200)
				return false;

			// Update label
			uint32_t pos = result & 0x100 ? pos2 : pos1;
			uint32_t old_label = labels[pos];
			labels[pos] = (uint8_t)(result & 0x3f);
			pos = pos * 4 + ((result >> 6) & 3);
			// Put elem
			if (old_label >= 0x10)
			{
				T victim = data[pos];
				data[pos] = elem;
				elem = victim;
			}
			else// bin is empty
			{
				data[pos] = elem;
				num_elems++;
				return true;
			}
		}
	}
};

// Label coded as followed:
//
// Minimum    Displacement
// ------    --------------
// |    |    |            |
// b5  b4    b3  b2  b1  b0
static void Create_Lookup_Table_2x4()
{
	FILE* file = fopen("cuckoo_LSAmax_lookup_2x4.cpp", "w");

	fprintf(file, "static uint16_t lt_2x4[] = {\n");

	// Iterate for all posible combinations of labels
	for (uint32_t label1 = 0; label1 < 64; label1++)
	{
		for (uint32_t label2 = 0; label2 < 64; label2++)
		{
			if (label1 == 61 && label2 == 60)
				label1 = 61;

			uint32_t labels[8];
			uint32_t load_bucket1 = 0;
			uint32_t load_bucket2 = 0;

			for (int i = 0; i < 4; i++)
			{
				// Decode labels
				labels[i + 0] = (label1 >> 4) + ((label1 >> i) & 1);
				labels[i + 4] = (label2 >> 4) + ((label2 >> i) & 1);

				load_bucket1 += labels[i + 0];
				load_bucket2 += labels[i + 4];
			}
			// Find minimun label
			uint32_t min_label = UINT32_MAX;
			uint32_t min_label2 = UINT32_MAX;

			for (int i = 0; i < 8; i++)
			{
				if (labels[i] < min_label)
				{
					min_label2 = min_label;
					min_label = labels[i];
				}
				else if (labels[i] < min_label2)
					min_label2 = labels[i];
			}

			// Resolves ties
			uint32_t pos;
			uint32_t displacement = load_bucket1 <= load_bucket2 ? 0 : 4;
			for (int i = 0; i < 4; i++)
				if (labels[i + displacement] == min_label)
				{
					labels[i + displacement] = min_label2 + 1;
					pos = i + displacement;
					break;
				}

			// Encode
			// Find minimun label
			min_label = UINT32_MAX;
			for (int i = 0; i < 4; i++)
				if (labels[i + displacement] < min_label)
					min_label = labels[i + displacement];

			// TODO: High bit to check if algorithm terminate isn't needed
			uint32_t result = ((label1 == 63 && label2 == 63) ? 0x200 : 0) | ((label1 == 63 && label2 == 63) ? 0 : pos << 6) | (std::min(min_label, 3u) << 4);
			min_label = std::min(min_label, 3u);
			for (int i = 0; i < 4; i++)
				result |= (labels[i + displacement] == min_label ? 0 : 1) << i;

			fprintf(file, "%s%u", label1 == 0 && label2 == 0 ? "" : ",", result);
		}
		fprintf(file, "\n");
	}

	fprintf(file, "};");

	fclose(file);
}

#include <chrono>

#define CUCKOO_CAPACITY 100'000
#define MAX_REPETITIONS 1000
#define MAX_L 7
struct TestData {
	double min;
	double max;
	double sum;
	uint32_t sum_labels;
	uint32_t unnecesary_moves;
	uint32_t num_moves;
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
		data[i].unnecesary_moves = 0;
		data[i].num_moves = 0;
	}

	// Seed with a real random value, if available
	std::random_device good_random;
	// Change first template parameter from 4 to other k to see other schemes
	Cuckoo_1t<4, uint64_t, true> cuckoo_table(CUCKOO_CAPACITY, good_random());

	auto start = std::chrono::high_resolution_clock::now();

	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS; repeat++)
	{
		uint32_t l_max = 0;
		std::mt19937_64 r(good_random());
		cuckoo_table.Clear();

		// Fill the table
		while (true)
		{
			//if (!cuckoo_table.Insert_RW_k2(r()))
			//if (!cuckoo_table.Insert_RW_k3(r()))
			if (!cuckoo_table.Insert_LSA_max_k2(r(), l_max + 1))
			//if (!cuckoo_table.Insert_LSA_max_k3(r(), l_max + 1))
			{
				double rate = cuckoo_table.GetTableUse();
				data[l_max].sum += rate;
				data[l_max].sum_labels += cuckoo_table.sum_labels;
				data[l_max].unnecesary_moves += cuckoo_table.same_bucket_count;
				data[l_max].num_moves += cuckoo_table.num_moves;
				if (data[l_max].min > rate)
					data[l_max].min = rate;
				if (data[l_max].max < rate)
					data[l_max].max = rate;

				// break;// For RW uncomment this break
				l_max++;
				if (l_max >= MAX_L)
					break;
			}
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("Time: %ums\n\n", (uint32_t)elapsed.count());

	printf("------------------------------------------------------------\n");
	printf("L_max          Table_Use           Labels   Total  Unnecesary\n");
	printf("     (min_diff-  avg   +max_diff)           Moves    Moves\n");
	printf("------------------------------------------------------------\n");
	for (uint32_t l_max = 0; l_max < MAX_L; l_max++)
	{
		double avg_use = data[l_max].sum / MAX_REPETITIONS;
		printf("%02u     %05.2f%% - %05.2f%% + %05.2f%%     %05.2f   %.2f     %.2f\n"
			, l_max + 1
			, avg_use - data[l_max].min, avg_use, data[l_max].max - avg_use
			, data[l_max].sum_labels*1. / MAX_REPETITIONS / CUCKOO_CAPACITY,
			data[l_max].num_moves*1. / MAX_REPETITIONS / CUCKOO_CAPACITY,
			data[l_max].unnecesary_moves*1. / MAX_REPETITIONS / CUCKOO_CAPACITY);
	}
}

#define MAX_MOVES_STEPS 20
void test_moves_done()
{
	double rates[MAX_MOVES_STEPS];
	memset(rates, 0, sizeof(rates));
	const uint32_t STEP = CUCKOO_CAPACITY / 10;// Iterate by 10%

	// Seed with a real random value, if available
	std::random_device good_random;
	// Change first template parameter from 4 to other k to see other schemes
	Cuckoo_1t<4, uint64_t, true> cuckoo_table(CUCKOO_CAPACITY, good_random());

	auto start = std::chrono::high_resolution_clock::now();

	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS; repeat++)
	{
		std::mt19937_64 r(good_random());
		cuckoo_table.Clear();
		uint32_t num_moves = 0;

		// Fill the table
		while (true)
		{
			if (!cuckoo_table.Insert_LSA_max_k2(r(), 4))
			{
				double rate = cuckoo_table.GetTableUse();
				while (num_moves < MAX_MOVES_STEPS)
				{
					rates[num_moves] += rate;
					num_moves++;
				}
				break;
			}
			else if (cuckoo_table.num_moves > (num_moves + 1)*STEP)
			{
				rates[num_moves] += cuckoo_table.GetTableUse();
				num_moves++;
				if (num_moves >= MAX_MOVES_STEPS)
					break;
			}
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("Time: %ums\n\n", (uint32_t)elapsed.count());

	// Show the histogram
	for (int i = 0; i < MAX_MOVES_STEPS; i++)
		printf("%.1f  %.2f\n", (i+1)/10., rates[i] / MAX_REPETITIONS);
}

#define MAX_REPETITIONS_ERRORS 10'000
void test_error()
{
	double rate[MAX_REPETITIONS_ERRORS];
	double avg_rate = 0;

	// Seed with a real random value, if available
	std::random_device good_random;
	// Change first template parameter from 4 to other k to see other schemes
	Cuckoo_1t<4, uint64_t, true> cuckoo_table(CUCKOO_CAPACITY, good_random());

	auto start = std::chrono::high_resolution_clock::now();

	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS_ERRORS; repeat++)
	{
		std::mt19937_64 r(good_random());
		cuckoo_table.Clear();

		// Fill the table
		while (true)
		{
			if (!cuckoo_table.Insert_LSA_max_k2(r(), 4))
			{
				rate[repeat] = cuckoo_table.GetTableUse();
				avg_rate += rate[repeat];
				break;
			}
		}
	}

	auto end = std::chrono::high_resolution_clock::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
	printf("Time: %ums\n\n", (uint32_t)elapsed.count());

	avg_rate /= MAX_REPETITIONS_ERRORS;
	double points[] = {0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 200};
	uint32_t counters[] = { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

	// Calculate differences
	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS_ERRORS; repeat++)
	{
		double diff = avg_rate - rate[repeat];

		for (int i = 0; i < (sizeof(points)/sizeof(points[0])); i++)
			if (diff < points[i])
			{
				counters[i]++;
				break;
			}
	}

	// Show the histogram
	for (int i = 0; i < (sizeof(points) / sizeof(points[0])); i++)
		printf("%03.1f  %u\n", points[i] , counters[i]);
}

void test_insertion_time()
{
	uint64_t times[100];
	memset(times, 0, sizeof(times));

	// Seed with a real random value, if available
	std::random_device good_random;

	for (uint32_t repeat = 0; repeat < MAX_REPETITIONS; repeat++)
	{
		std::mt19937_64 r(good_random());
		Cuckoo_2x4_1t_lookup<uint64_t> cuckoo_table_lookup(100'000);
		//Cuckoo_1t<4, uint64_t> cuckoo_table(100'000, good_random());

		// Fill the table
		int per_cent = 0;
		bool continue_cycle = true;
		while (continue_cycle)
		{
			auto start = std::chrono::high_resolution_clock::now();

			for (int i = 0; i < 1000; i++)
				if (!cuckoo_table_lookup.Insert(r()))
				//if (!cuckoo_table.Insert_RW_k2(r()))
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


void main()
{
	test_lmax();
	//test_moves_done();
	//test_error();
	//Create_Lookup_Table_2x4();
	//test_insertion_time();

	// Wait for one keystroke
	printf("\nPress any key to exit...");
	char c;
	scanf("%c", &c);
}