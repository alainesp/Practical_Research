///////////////////////////////////////////////////////////////////////////////
// Production code of Cuckoo Breeding Ground (CBG) hash table.
// See "research_cuckoo_cbg.md" for scientific paper explaining the options.
///////////////////////////////////////////////////////////////////////////////
//
// Written by Alain Espinosa <alainesp at gmail.com> in 2018 and placed
// under the MIT license (see LICENSE file for a full definition).
//
///////////////////////////////////////////////////////////////////////////////
//
// CBG version 0.1 (consider alpha version).
// Common operations are implemented but many, particularly special ones, are
// not. Many optimizations opportunities remains unexploited. This is intended as
// an early prototype for early adopters.
//
// Tested on Windows 64 bits with VC++ 2017. Requires C++ 14.
// Easy to support 32 bits, Linux and other compilers, but currently not done.
// Comments and PR are welcomed.
//
///////////////////////////////////////////////////////////////////////////////
// Code Architecture
//-----------------------------------------------------------------------------
//
// Class 'cbg::cbg_internal::CBG_IMPL' contains the main algorithm with common
// operations. It doesn't depends in any particular data layout. Class
// 'cbg::cbg_internal::CBG_MAP_IMPL' only add mapping operations.
//
// Data management is provided by the DATA and METADATA template parameters.
// Currently we provide 3 different layouts:
//
// - "Struct of Arrays" (SoA): Metadata overhead of 2 bytes. Each metadata, 
//    keys and values are on a different array. The fastest for negative 
//    queries.
//
// - "Array of Structs" (AoS): Metadata overhead of one byte. Each metadata, 
//    keys and values are in the same array each one after the other. Heavy use
//    of unaligned memory access. The fastest for positive queries.
//
// - "Array of Blocks" (AoB): Metadata overhead of one byte. Each metadata, 
//    keys and values are in the same array on blocks. Don't use unaligned
//    memory access but performance suffers. Intended for positive queries.
//
// Namespace 'cbg::hashing' contains non-cryptographic hashing based on t1ha2
// [https://github.com/leo-yuriev/t1ha]. This is for 64 bits only. For 32 bits
// we need a different hash function for performance reason.
//
///////////////////////////////////////////////////////////////////////////////
// NUM_ELEMS_BUCKET parameter selection:
//
// - Fastest possible : NUM_ELEMS_BUCKET=2 with load_factor < 80%
// - No memory waste  : NUM_ELEMS_BUCKET=4 with load_factor > 95%
// - Balanced approach: NUM_ELEMS_BUCKET=3 with 80% < load_factor < 95%
///////////////////////////////////////////////////////////////////////////////

#pragma once

// Disable assert() and obtain best performance
#define NDEBUG

#include <cstdint>
#include <tuple>
#include <cassert>
#include <random>

#if defined(_MSC_VER) && defined (_WIN64)
#include <intrin.h>// should be part of all recent Visual Studio
#pragma intrinsic(_umul128)
#pragma intrinsic(__umulh)
#define mul_64x64_128(a, b, ph) _umul128(a, b, ph)
#define rot64(v, s) _rotr64(v, s)
#endif // defined(_MSC_VER) && defined (_WIN64)

namespace cbg
{
///////////////////////////////////////////////////////////////////////////////
// Non-cryptographic hashing for hashtables. A clearer C++ version of t1ha2
///////////////////////////////////////////////////////////////////////////////
namespace hashing
{
namespace t1ha2_internal
{
///////////////////////////////////////////////////////////////////////////////
// Data Access
///////////////////////////////////////////////////////////////////////////////
struct x86//Little Endian, Fast unaligned memory access
{
	static __forceinline uint64_t fetch64(const void* v)
	{
		return *((const uint64_t*)v);
	}
	static __forceinline uint64_t tail64(const void *v, size_t tail)
	{
		// On some systems (e.g. x86) we can perform a 'oneshot' read,
		// which is a little bit faster.
		const unsigned offset = (8 - tail) & 7;
		const unsigned shift = offset << 3;
		return fetch64(v) & (UINT64_MAX >> shift);
	}
};

///////////////////////////////////////////////////////////////////////////////
// t1ha2 base
///////////////////////////////////////////////////////////////////////////////
template<class DATA_ACCESS = x86> struct t1ha2_IMPL : public DATA_ACCESS
{
	// 'Magic' primes
	static constexpr uint64_t prime_0 = UINT64_C(0xEC99BF0D8372CAAB);
	static constexpr uint64_t prime_1 = UINT64_C(0x82434FE90EDCEF39);
	static constexpr uint64_t prime_2 = UINT64_C(0xD4F06DB99D67BE4B);
	static constexpr uint64_t prime_3 = UINT64_C(0xBD9CACC22C6E9571);
	static constexpr uint64_t prime_4 = UINT64_C(0x9C06FAF4D023E3AB);
	static constexpr uint64_t prime_5 = UINT64_C(0xC060724A8424F345);
	static constexpr uint64_t prime_6 = UINT64_C(0xCB5AF53AE3AAAC31);

	uint64_t seed;

	t1ha2_IMPL() noexcept : DATA_ACCESS()
	{
		// Seed with a real random value, if available
		std::random_device good_random;
		seed = static_cast<uint64_t>(good_random()) | static_cast<uint64_t>(good_random()) << 32;
	}
	t1ha2_IMPL(uint64_t seed) noexcept : DATA_ACCESS(), seed(seed)
	{}

	uint64_t operator()(const void* data, size_t length) const noexcept
	{
		// Init a,b
		uint64_t a = seed;
		uint64_t b = length;

		if (length > 32)
		{
			// Init c,d
			uint64_t c = rot64(length, 23) + ~seed;
			uint64_t d = ~length + rot64(seed, 19);

			// T1HA2_LOOP
			const void* detent = (const uint8_t*)data + length - 31;
			do {
				const uint64_t* v = (const uint64_t*)data;
				data = (const uint64_t*)data + 4;
				//prefetch(data);
				// T1HA2_UPDATE
				const uint64_t w0 = DATA_ACCESS::fetch64(v + 0);
				const uint64_t w1 = DATA_ACCESS::fetch64(v + 1);
				const uint64_t w2 = DATA_ACCESS::fetch64(v + 2);
				const uint64_t w3 = DATA_ACCESS::fetch64(v + 3);

				const uint64_t d02 = w0 + rot64(w2 + d, 56);
				const uint64_t c13 = w1 + rot64(w3 + c, 19);
				d ^= b + rot64(w1, 38);
				c ^= a + rot64(w0, 57);
				b ^= prime_6 * (c13 + w2);
				a ^= prime_5 * (d02 + w3);
			} while (data < detent);

			// Squash
			a ^= prime_6 * (c + rot64(d, 23));
			b ^= prime_5 * (rot64(c, 19) + d);
			length &= 31;
		}
		// T1HA2_TAIL_AB
		const uint64_t *v = (const uint64_t *)data;
		uint64_t l, h;
		switch (length)
		{
		default:// mixup64
			a ^= mul_64x64_128(b + DATA_ACCESS::fetch64(v++), prime_4, &h);
			b += h;
			//[[fallthrough]];
		case 24: case 23: case 22: case 21: case 20: case 19: case 18: case 17:
			// mixup64
			b ^= mul_64x64_128(a + DATA_ACCESS::fetch64(v++), prime_3, &h);
			a += h;
			//[[fallthrough]];
		case 16: case 15: case 14: case 13: case 12: case 11: case 10: case 9:
			// mixup64
			a ^= mul_64x64_128(b + DATA_ACCESS::fetch64(v++), prime_2, &h);
			b += h;
			//[[fallthrough]];
		case 8: case 7: case 6: case 5: case 4: case 3: case 2: case 1:
			// mixup64
			b ^= mul_64x64_128(a + DATA_ACCESS::tail64(v, length), prime_1, &h);
			a += h;
			//[[fallthrough]];
		case 0: break;
		}
		// final64(a, b);
		uint64_t x = (a + rot64(b, 41)) * prime_0;
		uint64_t y = (rot64(a, 23) + b) * prime_6;
		l = mul_64x64_128(x ^ y, prime_5, &h);
		return l ^ h;
	}
};
}// end namespace t1ha2_internal

///////////////////////////////////////////////////////////////////////////////
// Common t1ha2 for general use
///////////////////////////////////////////////////////////////////////////////
template<class T, class DATA_ACCESS = t1ha2_internal::x86> struct t1ha2 : public t1ha2_internal::t1ha2_IMPL<DATA_ACCESS>
{
	t1ha2() noexcept : t1ha2_internal::t1ha2_IMPL<DATA_ACCESS>()
	{}
	t1ha2(uint64_t seed) noexcept : t1ha2_internal::t1ha2_IMPL<DATA_ACCESS>(seed)
	{}

	template<size_t length = sizeof(T)> uint64_t operator()(const T& elem) const noexcept
	{
		const uint64_t* v = (uint64_t*)(&elem);
		// Init a,b
		uint64_t a = seed;
		uint64_t b = length;

		if (length > 32)
		{
			// Init c,d
			uint64_t c = rot64(length, 23) + ~seed;
			uint64_t d = ~length + rot64(seed, 19);

			// T1HA2_LOOP
			const void* detent = (const uint8_t*)v + length - 31;
			do {
				v += 4;
				//prefetch(v);
				// T1HA2_UPDATE
				const uint64_t w0 = DATA_ACCESS::fetch64(v + 0);
				const uint64_t w1 = DATA_ACCESS::fetch64(v + 1);
				const uint64_t w2 = DATA_ACCESS::fetch64(v + 2);
				const uint64_t w3 = DATA_ACCESS::fetch64(v + 3);

				const uint64_t d02 = w0 + rot64(w2 + d, 56);
				const uint64_t c13 = w1 + rot64(w3 + c, 19);
				d ^= b + rot64(w1, 38);
				c ^= a + rot64(w0, 57);
				b ^= prime_6 * (c13 + w2);
				a ^= prime_5 * (d02 + w3);
			} while (v < detent);

			// Squash
			a ^= prime_6 * (c + rot64(d, 23));
			b ^= prime_5 * (rot64(c, 19) + d);
		}
		// T1HA2_TAIL_AB
		uint64_t l, h;
		switch (length & 31)
		{
		default:// mixup64
			a ^= mul_64x64_128(b + DATA_ACCESS::fetch64(v++), prime_4, &h);
			b += h;
			//[[fallthrough]];
		case 24: case 23: case 22: case 21: case 20: case 19: case 18: case 17:
			// mixup64
			b ^= mul_64x64_128(a + DATA_ACCESS::fetch64(v++), prime_3, &h);
			a += h;
			//[[fallthrough]];
		case 16: case 15: case 14: case 13: case 12: case 11: case 10: case 9:
			// mixup64
			a ^= mul_64x64_128(b + DATA_ACCESS::fetch64(v++), prime_2, &h);
			b += h;
			//[[fallthrough]];
		case 8: case 7: case 6: case 5: case 4: case 3: case 2: case 1:
			// mixup64
			b ^= mul_64x64_128(a + DATA_ACCESS::tail64(v, length & 31), prime_1, &h);
			a += h;
			//[[fallthrough]];
		case 0: break;
		}
		// final64(a, b);
		uint64_t x = (a + rot64(b, 41)) * prime_0;
		uint64_t y = (rot64(a, 23) + b) * prime_6;
		l = mul_64x64_128(x ^ y, prime_5, &h);
		return l ^ h;
	}
};
// Partial specialization
template<class DATA_ACCESS> struct t1ha2<std::string, DATA_ACCESS> : public t1ha2_internal::t1ha2_IMPL<DATA_ACCESS>
{
	t1ha2() noexcept : t1ha2_internal::t1ha2_IMPL<DATA_ACCESS>()
	{}
	t1ha2(uint64_t seed) noexcept : t1ha2_internal::t1ha2_IMPL<DATA_ACCESS>(seed)
	{}

	uint64_t operator()(const std::string& data) const noexcept
	{
		return operator(data.c_str(), data.length);
	}
};
template<class DATA_ACCESS> struct t1ha2<char*, DATA_ACCESS> : public t1ha2_internal::t1ha2_IMPL<DATA_ACCESS>
{
	t1ha2() noexcept : t1ha2_internal::t1ha2_IMPL<DATA_ACCESS>()
	{}
	t1ha2(uint64_t seed) noexcept : t1ha2_internal::t1ha2_IMPL<DATA_ACCESS>(seed)
	{}

	uint64_t operator()(const char* data) const noexcept
	{
		return operator(data, strlen(data));
	}
};
}// end namespace hashing

// Search Hint
enum class Search_Hint
{
	Unknow,
	Expect_Positive,
	Expect_Negative
};
// Internal implementations
namespace cbg_internal
{
///////////////////////////////////////////////////////////////////////////////
// Data layout is "Struct of Arrays"
///////////////////////////////////////////////////////////////////////////////
// Metadata layout
struct MetadataLayout_SoA
{
	uint16_t* metadata;

	MetadataLayout_SoA() noexcept : metadata(nullptr)
	{}
	MetadataLayout_SoA(size_t num_bins) noexcept
	{
		metadata = (uint16_t*)malloc(num_bins * sizeof(uint16_t));
		memset(metadata, 0, num_bins * sizeof(uint16_t));
	}
	~MetadataLayout_SoA() noexcept
	{
		free(metadata);
		metadata = nullptr;
	}
	__forceinline void Clear(size_t initial_pos, size_t size_in_bins) noexcept
	{
		memset(metadata + initial_pos, 0, size_in_bins * sizeof(uint16_t));
	}
	__forceinline void ReallocMetadata(size_t new_num_bins) noexcept
	{
		metadata = (uint16_t*)realloc(metadata, new_num_bins * sizeof(uint16_t));
	}

	/////////////////////////////////////////////////////////////////////
	// Metadata coded utilities
	/////////////////////////////////////////////////////////////////////
	//
	// Unlucky  Bucket is   Element        
	// bucket   Reversed    Distance     Labels
	// -------  ----------  --------    -------- 
	//   |          |       |      |    |      |
	//   b7         b6      b5 b4 b3    b2 b1 b0
	// 0b00'000'000
	__forceinline uint16_t at(size_t pos) const noexcept
	{
		return metadata[pos];
	}
	__forceinline uint16_t Get_Label(size_t pos) const noexcept
	{
		return metadata[pos] & 0b00'000'111u;
	}
	__forceinline bool Is_Empty(size_t pos) const noexcept
	{
		return Get_Label(pos) == 0;
	}
	__forceinline void Set_Empty(size_t pos) noexcept
	{
		metadata[pos] &= 0b11'000'000;
	}
	__forceinline uint16_t Get_Hash(size_t pos) const noexcept
	{
		return metadata[pos] & 0xFF00u;
	}
	__forceinline void Update_Bin_At(size_t pos, size_t distance_to_base, bool is_reverse_item, uint_fast16_t label, size_t hash) noexcept
	{
		metadata[pos] = uint16_t((hash & 0xFF00) | (metadata[pos] & 0b11'000'000) | (is_reverse_item ? 0b00'100'000 : 0) | (distance_to_base << 3) | label);
	}
	__forceinline bool Is_Item_In_Reverse_Bucket(size_t pos) const noexcept
	{
		return metadata[pos] & 0b00'100'000;
	}
	__forceinline uint16_t Distance_to_Entry_Bin(size_t pos) const noexcept
	{
		return (metadata[pos] >> 3u) & 0b11u;
	}
	/*__forceinline bool Is_Unlucky_Bucket(size_t pos) const noexcept
	{
	return metadata[pos] & 0b10'000'000;
	}*/
	__forceinline void Set_Unlucky_Bucket(size_t pos) noexcept
	{
		metadata[pos] |= 0b10'000'000;
	}
	__forceinline bool Is_Bucket_Reversed(size_t pos) const noexcept
	{
		return metadata[pos] & 0b01'000'000;
	}
	__forceinline void Set_Bucket_Reversed(size_t pos) noexcept
	{
		metadata[pos] |= 0b01'000'000;
	}

	//// Cache line aware
	//__forceinline bool Benefit_With_Reversal(size_t pos, size_t size_bucket) const noexcept
	//{
	//	uintptr_t begin_cache_line = reinterpret_cast<uintptr_t>(metadata + pos) & (~static_cast<uintptr_t>(63));
	//	uintptr_t end_cache_line = reinterpret_cast<uintptr_t>(metadata + pos + size_bucket - 1) & (~static_cast<uintptr_t>(63));

	//	return begin_cache_line < end_cache_line;
	//}
};
// Data layouts
template<class KEY> struct KeyLayout_SoA : public MetadataLayout_SoA
{
	KEY* keys;

	// Constructors
	KeyLayout_SoA() noexcept : keys(nullptr), MetadataLayout_SoA()
	{}
	KeyLayout_SoA(size_t num_buckets) noexcept : MetadataLayout_SoA(num_buckets)
	{
		keys = (KEY*)malloc(num_buckets * sizeof(KEY));
	}
	~KeyLayout_SoA() noexcept
	{
		free(keys);
		keys = nullptr;
	}

	__forceinline void MoveElem(size_t dest, size_t orig) noexcept
	{
		keys[dest] = keys[orig];
	}
	__forceinline void SaveElem(size_t pos, const KEY& elem) noexcept
	{
		keys[pos] = elem;
	}

	__forceinline const KEY& GetKey(size_t pos) const noexcept
	{
		return keys[pos];
	}
	__forceinline const KEY& GetKeyFromValue(const KEY& elem) const noexcept
	{
		return elem;
	}
	__forceinline KEY GetElem(size_t pos) const noexcept
	{
		return keys[pos];
	}
	__forceinline KEY* GetValue(size_t pos) const noexcept
	{
		return keys + pos;
	}

	__forceinline void ReallocElems(size_t new_num_buckets) noexcept
	{
		keys = (KEY*)realloc(keys, new_num_buckets * sizeof(KEY));
	}
};
template<class KEY, class T> struct MapLayout_SoA : public MetadataLayout_SoA
{
	using INSERT_TYPE = std::pair<KEY, T>;

	KEY* keys;
	T* data;

	// Constructors
	MapLayout_SoA() noexcept : keys(nullptr), data(nullptr), MetadataLayout_SoA()
	{}
	MapLayout_SoA(size_t num_buckets) noexcept : MetadataLayout_SoA(num_buckets)
	{
		keys = (KEY*)malloc(num_buckets * sizeof(KEY));
		data = (T*)malloc(num_buckets * sizeof(T));
	}
	~MapLayout_SoA() noexcept
	{
		free(keys);
		free(data);
		keys = nullptr;
		data = nullptr;
	}

	__forceinline void MoveElem(size_t dest, size_t orig) noexcept
	{
		keys[dest] = keys[orig];
		data[dest] = data[orig];
	}
	__forceinline void SaveElem(size_t pos, const INSERT_TYPE& elem) noexcept
	{
		keys[pos] = elem.first;
		data[pos] = elem.second;
	}

	__forceinline const KEY& GetKey(size_t pos) const noexcept
	{
		return keys[pos];
	}
	__forceinline const KEY& GetKeyFromValue(const INSERT_TYPE& elem) const noexcept
	{
		return elem.first;
	}
	__forceinline INSERT_TYPE GetElem(size_t pos) const noexcept
	{
		return std::make_pair(keys[pos], data[pos]);
	}
	__forceinline T* GetValue(size_t pos) const noexcept
	{
		return data + pos;
	}

	__forceinline void ReallocElems(size_t new_num_buckets) noexcept
	{
		keys = (KEY*)realloc(keys, new_num_buckets * sizeof(KEY));
		data = (T*)realloc(data, new_num_buckets * sizeof(T));
	}
};

///////////////////////////////////////////////////////////////////////////////
// Data layout is "Array of Structs"
///////////////////////////////////////////////////////////////////////////////
// Auxiliary class
template<size_t ELEM_SIZE> struct ElemLayout
{
	uint8_t metadata;
	uint8_t elem[ELEM_SIZE];
};
// Metadata layout
template<size_t ELEM_SIZE> struct MetadataLayout_AoS
{
	ElemLayout<ELEM_SIZE>* all_data;

	MetadataLayout_AoS() noexcept : all_data(nullptr)
	{
	}
	MetadataLayout_AoS(size_t num_bins) noexcept
	{
		all_data = (ElemLayout<ELEM_SIZE>*)malloc(num_bins * sizeof(ElemLayout<ELEM_SIZE>));
		for (size_t i = 0; i < num_bins; i++)
			all_data[i].metadata = 0;
	}
	~MetadataLayout_AoS() noexcept
	{
		free(all_data);
		all_data = nullptr;
	}
	__forceinline void Clear(size_t initial_pos, size_t size_in_bins) noexcept
	{
		for (size_t i = initial_pos; i < (initial_pos + size_in_bins); i++)
			all_data[i].metadata = 0;
	}
	__forceinline void ReallocMetadata(size_t new_num_bins) noexcept
	{
		all_data = (ElemLayout<ELEM_SIZE>*)realloc(all_data, new_num_bins * sizeof(ElemLayout<ELEM_SIZE>));
	}

	/////////////////////////////////////////////////////////////////////
	// Metadata coded utilities
	/////////////////////////////////////////////////////////////////////
	//
	// Unlucky  Bucket is   Element        
	// bucket   Reversed    Distance     Labels
	// -------  ----------  --------    -------- 
	//   |          |       |      |    |      |
	//   b7         b6      b5 b4 b3    b2 b1 b0
	// 0b00'000'000
	__forceinline uint8_t at(size_t pos) const noexcept
	{
		return all_data[pos].metadata;
	}
	__forceinline uint8_t Get_Label(size_t pos) const noexcept
	{
		return all_data[pos].metadata & 0b00'000'111;
	}
	__forceinline bool Is_Empty(size_t pos) const noexcept
	{
		return Get_Label(pos) == 0;
	}
	__forceinline void Set_Empty(size_t pos) noexcept
	{
		all_data[pos].metadata &= 0b11'000'000;
	}
	__forceinline uint16_t Get_Hash(size_t pos) const noexcept
	{
		return 0;
	}
	__forceinline void Update_Bin_At(size_t pos, size_t distance_to_base, bool is_reverse_item, uint_fast16_t label, size_t hash) noexcept
	{
		all_data[pos].metadata = (all_data[pos].metadata & 0b11'000'000) | (is_reverse_item ? 0b00'100'000 : 0) | (uint8_t(distance_to_base) << 3) | label;
	}
	__forceinline bool Is_Item_In_Reverse_Bucket(size_t pos) const noexcept
	{
		return all_data[pos].metadata & 0b00'100'000;
	}
	__forceinline uint16_t Distance_to_Entry_Bin(size_t pos) const noexcept
	{
		return (all_data[pos].metadata >> 3) & 0b11;
	}
	/*__forceinline bool Is_Unlucky_Bucket(size_t pos) const noexcept
	{
	return all_data[pos].metadata & 0b10'000'000;
	}*/
	__forceinline void Set_Unlucky_Bucket(size_t pos) noexcept
	{
		all_data[pos].metadata |= 0b10'000'000;
	}
	__forceinline bool Is_Bucket_Reversed(size_t pos) const noexcept
	{
		return all_data[pos].metadata & 0b01'000'000;
	}
	__forceinline void Set_Bucket_Reversed(size_t pos) noexcept
	{
		all_data[pos].metadata |= 0b01'000'000;
	}

	//// Cache line aware
	//__forceinline bool Benefit_With_Reversal(size_t pos, size_t size_bucket) const noexcept
	//{
	//	uintptr_t begin_cache_line = reinterpret_cast<uintptr_t>(all_data + pos) & (~static_cast<uintptr_t>(63));
	//	uintptr_t end_cache_line = reinterpret_cast<uintptr_t>(all_data + pos + size_bucket - 1) & (~static_cast<uintptr_t>(63));
	//	uintptr_t reverse_cache_line = reinterpret_cast<uintptr_t>(all_data + pos - (size_bucket - 1)) & (~static_cast<uintptr_t>(63));

	//	// TODO: Do more checks, this is too simple
	//	return begin_cache_line < end_cache_line && reverse_cache_line == begin_cache_line;
	//}
};
// Data layouts
template<class KEY> struct KeyLayout_AoS : public MetadataLayout_AoS<sizeof(KEY)>
{
	// Constructors
	KeyLayout_AoS() noexcept : MetadataLayout_AoS<sizeof(KEY)>()
	{}
	KeyLayout_AoS(size_t num_bins) noexcept : MetadataLayout_AoS<sizeof(KEY)>(num_bins)
	{}

	__forceinline void MoveElem(size_t dest, size_t orig) noexcept
	{
		memcpy(all_data[dest].elem, all_data[orig].elem, sizeof(KEY));
	}
	__forceinline void SaveElem(size_t pos, const KEY& elem) noexcept
	{
		memcpy(all_data[pos].elem, &elem, sizeof(KEY));
	}

	__forceinline const KEY& GetKey(size_t pos) const noexcept
	{
		return *((KEY*)(all_data[pos].elem));
	}
	__forceinline const KEY& GetKeyFromValue(const KEY& elem) const noexcept
	{
		return elem;
	}
	__forceinline KEY GetElem(size_t pos) const noexcept
	{
		return *((KEY*)(all_data[pos].elem));
	}
	__forceinline KEY* GetValue(size_t pos) const noexcept
	{
		return (KEY*)(all_data[pos].elem);
	}

	__forceinline void ReallocElems(size_t new_num_buckets) noexcept
	{
		// Nothing
	}
};
template<class KEY, class T> struct MapLayout_AoS : public MetadataLayout_AoS<sizeof(KEY) + sizeof(T)>
{
	using INSERT_TYPE = std::pair<KEY, T>;

	// Constructors
	MapLayout_AoS() noexcept : MetadataLayout_AoS<sizeof(KEY) + sizeof(T)>()
	{}
	MapLayout_AoS(size_t num_buckets) noexcept : MetadataLayout_AoS<sizeof(KEY) + sizeof(T)>(num_buckets)
	{}

	__forceinline void MoveElem(size_t dest, size_t orig) noexcept
	{
		memcpy(all_data[dest].elem, all_data[orig].elem, sizeof(KEY) + sizeof(T));
	}
	__forceinline void SaveElem(size_t pos, const INSERT_TYPE& elem) noexcept
	{
		memcpy(all_data[pos].elem, &elem.first, sizeof(KEY));
		memcpy(all_data[pos].elem + sizeof(KEY), &elem.second, sizeof(T));
	}

	__forceinline const KEY& GetKey(size_t pos) const noexcept
	{
		return *((KEY*)(all_data[pos].elem));
	}
	__forceinline const KEY& GetKeyFromValue(const INSERT_TYPE& elem) const noexcept
	{
		return elem.first;
	}
	__forceinline INSERT_TYPE GetElem(size_t pos) const noexcept
	{
		return std::make_pair(*((KEY*)(all_data[pos].elem)), *((T*)(all_data[pos].elem + sizeof(KEY))));
	}
	__forceinline T* GetValue(size_t pos) const noexcept
	{
		return (T*)(all_data[pos].elem + sizeof(KEY));
	}

	__forceinline void ReallocElems(size_t new_num_buckets) noexcept
	{
		// Nothing
	}
};

///////////////////////////////////////////////////////////////////////////////
// Data layout is "Array of Blocks"
///////////////////////////////////////////////////////////////////////////////
// Auxiliary classes
template<class T> struct BlockKey
{
	static constexpr size_t BLOCK_SIZE = alignof(T);

	uint8_t metadata[BLOCK_SIZE];
	T data[BLOCK_SIZE];
};
template<class KEY, class T> struct MaxAlignOf
{
	static constexpr size_t BLOCK_SIZE = std::max(alignof(KEY), alignof(T));
};
template<class KEY, class T> struct BlockMap
{
	uint8_t metadata[MaxAlignOf<KEY, T>::BLOCK_SIZE];
	KEY keys[MaxAlignOf<KEY, T>::BLOCK_SIZE];
	T data[MaxAlignOf<KEY, T>::BLOCK_SIZE];
};
// Metadata layout
template<size_t BLOCK_SIZE, class BLOCK> struct MetadataLayout_AoB
{
	BLOCK* all_data;

	MetadataLayout_AoB() noexcept : all_data(nullptr)
	{}
	MetadataLayout_AoB(size_t num_bins) noexcept
	{
		num_bins = (num_bins + BLOCK_SIZE - 1) / BLOCK_SIZE;

		all_data = (BLOCK*)malloc(num_bins * sizeof(BLOCK));
		for (size_t i = 0; i < num_bins; i++)
			for (size_t j = 0; j < BLOCK_SIZE; j++)
				all_data[i].metadata[j] = 0;
	}
	~MetadataLayout_AoB() noexcept
	{
		free(all_data);
		all_data = nullptr;
	}
	__forceinline void Clear(size_t initial_pos, size_t size_in_bins) noexcept
	{
		for (size_t i = initial_pos; i < (initial_pos + size_in_bins); i++)
			all_data[i / BLOCK_SIZE].metadata[i % BLOCK_SIZE] = 0;
	}
	__forceinline void ReallocMetadata(size_t new_num_bins) noexcept
	{
		new_num_bins = (new_num_bins + BLOCK_SIZE - 1) / BLOCK_SIZE;
		all_data = (BLOCK*)realloc(all_data, new_num_bins * sizeof(BLOCK));
	}

	/////////////////////////////////////////////////////////////////////
	// Metadata coded utilities
	/////////////////////////////////////////////////////////////////////
	//
	// Unlucky  Bucket is   Element        
	// bucket   Reversed    Distance     Labels
	// -------  ----------  --------    -------- 
	//   |          |       |      |    |      |
	//   b7         b6      b5 b4 b3    b2 b1 b0
	// 0b00'000'000
	__forceinline uint8_t at(size_t pos) const noexcept
	{
		return all_data[pos / BLOCK_SIZE].metadata[pos%BLOCK_SIZE];
	}
	__forceinline uint8_t Get_Label(size_t pos) const noexcept
	{
		return at(pos) & 0b00'000'111u;
	}
	__forceinline bool Is_Empty(size_t pos) const noexcept
	{
		return Get_Label(pos) == 0;
	}
	__forceinline void Set_Empty(size_t pos) noexcept
	{
		all_data[pos / BLOCK_SIZE].metadata[pos%BLOCK_SIZE] &= 0b11'000'000;
	}
	__forceinline uint16_t Get_Hash(size_t /*pos*/) const noexcept
	{
		return 0;
	}
	__forceinline void Update_Bin_At(size_t pos, size_t distance_to_base, bool is_reverse_item, uint_fast16_t label, size_t /*hash*/) noexcept
	{
		all_data[pos / BLOCK_SIZE].metadata[pos%BLOCK_SIZE] = uint8_t((at(pos) & 0b11'000'000) | (is_reverse_item ? 0b00'100'000 : 0) | (distance_to_base << 3) | label);
	}
	__forceinline bool Is_Item_In_Reverse_Bucket(size_t pos) const noexcept
	{
		return at(pos) & 0b00'100'000;
	}
	__forceinline uint16_t Distance_to_Entry_Bin(size_t pos) const noexcept
	{
		return (at(pos) >> 3u) & 0b11u;
	}
	/*__forceinline bool Is_Unlucky_Bucket(size_t pos) const noexcept
	{
	return at(pos) & 0b10'000'000;
	}*/
	__forceinline void Set_Unlucky_Bucket(size_t pos) noexcept
	{
		all_data[pos / BLOCK_SIZE].metadata[pos%BLOCK_SIZE] |= 0b10'000'000;
	}
	__forceinline bool Is_Bucket_Reversed(size_t pos) const noexcept
	{
		return at(pos) & 0b01'000'000;
	}
	__forceinline void Set_Bucket_Reversed(size_t pos) noexcept
	{
		all_data[pos / BLOCK_SIZE].metadata[pos%BLOCK_SIZE] |= 0b01'000'000;
	}

	//// Cache line aware
	//__forceinline bool Benefit_With_Reversal(size_t pos, size_t size_bucket) const noexcept
	//{
	//	// TODO: Make calculations here
	//	return false;
	//}
};
// Data layouts
template<class KEY> struct KeyLayout_AoB : public MetadataLayout_AoB<alignof(KEY), BlockKey<KEY>>
{
	static constexpr size_t BLOCK_SIZE = alignof(KEY);

	// Constructors
	KeyLayout_AoB() noexcept : MetadataLayout_AoB()
	{}
	KeyLayout_AoB(size_t num_bins) noexcept : MetadataLayout_AoB(num_bins)
	{}

	__forceinline void MoveElem(size_t dest, size_t orig) noexcept
	{
		all_data[dest / BLOCK_SIZE].data[dest%BLOCK_SIZE] = all_data[orig / BLOCK_SIZE].data[orig%BLOCK_SIZE];
	}
	__forceinline void SaveElem(size_t pos, const KEY& elem) noexcept
	{
		all_data[pos / BLOCK_SIZE].data[pos%BLOCK_SIZE] = elem;
	}

	__forceinline const KEY& GetKey(size_t pos) const noexcept
	{
		return all_data[pos / BLOCK_SIZE].data[pos%BLOCK_SIZE];
	}
	__forceinline const KEY& GetKeyFromValue(const KEY& elem) const noexcept
	{
		return elem;
	}
	__forceinline KEY GetElem(size_t pos) const noexcept
	{
		return all_data[pos / BLOCK_SIZE].data[pos%BLOCK_SIZE];
	}
	__forceinline KEY* GetValue(size_t pos) const noexcept
	{
		return all_data[pos / BLOCK_SIZE].data + pos%BLOCK_SIZE;
	}

	__forceinline void ReallocElems(size_t new_num_buckets) noexcept
	{
		// Nothing
	}
};
template<class KEY, class T> struct MapLayout_AoB : public MetadataLayout_AoB<MaxAlignOf<KEY, T>::BLOCK_SIZE, BlockMap<KEY, T>>
{
	using INSERT_TYPE = std::pair<KEY, T>;
	static constexpr size_t BLOCK_SIZE = std::max(alignof(KEY), alignof(T));

	// Constructors
	MapLayout_AoB() noexcept : MetadataLayout_AoB()
	{}
	MapLayout_AoB(size_t num_buckets) noexcept : MetadataLayout_AoB(num_buckets)
	{}

	__forceinline void MoveElem(size_t dest, size_t orig) noexcept
	{
		all_data[dest / BLOCK_SIZE].keys[dest%BLOCK_SIZE] = all_data[orig / BLOCK_SIZE].keys[orig%BLOCK_SIZE];
		all_data[dest / BLOCK_SIZE].data[dest%BLOCK_SIZE] = all_data[orig / BLOCK_SIZE].data[orig%BLOCK_SIZE];
	}
	__forceinline void SaveElem(size_t pos, const INSERT_TYPE& elem) noexcept
	{
		all_data[pos / BLOCK_SIZE].keys[pos%BLOCK_SIZE] = elem.first;
		all_data[pos / BLOCK_SIZE].data[pos%BLOCK_SIZE] = elem.second;
	}

	__forceinline const KEY& GetKey(size_t pos) const noexcept
	{
		return all_data[pos / BLOCK_SIZE].keys[pos%BLOCK_SIZE];
	}
	__forceinline const KEY& GetKeyFromValue(const INSERT_TYPE& elem) const noexcept
	{
		return elem.first;
	}
	__forceinline INSERT_TYPE GetElem(size_t pos) const noexcept
	{
		return std::make_pair(all_data[pos / BLOCK_SIZE].keys[pos%BLOCK_SIZE], all_data[pos / BLOCK_SIZE].data[pos%BLOCK_SIZE]);
	}
	__forceinline T* GetValue(size_t pos) const noexcept
	{
		return all_data[pos / BLOCK_SIZE].data + pos%BLOCK_SIZE;
	}

	__forceinline void ReallocElems(size_t /*new_num_buckets*/) noexcept
	{
		// Nothing
	}
};

///////////////////////////////////////////////////////////////////////////////
// Basic implementation of CBG.
//
// TODO: Iterator, Handle move semantic of elems, destructor call when removed,
//       Memory_Allocator
///////////////////////////////////////////////////////////////////////////////
template<size_t NUM_ELEMS_BUCKET, class INSERT_TYPE, class KEY_TYPE, class VALUE_TYPE, class HASHER, class EQ, class DATA, class METADATA, bool USE_METADATA_SEARCH> class CBG_IMPL : private HASHER, private EQ, protected DATA
{
protected:
	// Pointers -> found in DATA and METADATA through inheritance
	// Counters
	size_t num_elems;
	size_t num_buckets;
	// Parameters
	float _max_load_factor = 0.9001f;// 90% -> When this load factor is reached the table grow
	float _grow_factor = 1.2f;// 20% -> How much to grow the table
	// Constants
	static constexpr uint_fast16_t L_MAX = 7;
	static constexpr size_t MIN_BUCKETS_COUNT = 2 * NUM_ELEMS_BUCKET - 2;

	static_assert(NUM_ELEMS_BUCKET >= 2 && NUM_ELEMS_BUCKET <= 4, "To use only 2 bits");
	static_assert(std::is_unsigned<size_t>::value, "size_t required to be unsigned");

	///////////////////////////////////////////////////////////////////////////
	// Utilities
	///////////////////////////////////////////////////////////////////////////
	__forceinline bool cmp_elems(size_t pos, const KEY_TYPE& r) const noexcept
	{
		return EQ::operator()(DATA::GetKey(pos), r);
	}
	// Note: We use a hash of only 64 bits for the two hash functions as:
	//
	// hash0 = hash;
	// hash1 = rot64(hash, 32);
	//
	// Hashtable capacity of more than 32 bits may give problems. Cuckoo filter
	// [[1]] demonstrates that cuckoo algorithm works similarly with 5-6 bits
	// for the secondary hash than with two fully independent hash functions.
	// Therefore we expect (TODO: check it) CBG to work well until the number
	// of buckets is 6 bytes or 255 Tera (of buckets, if each element is a
	// uint64_t this will use 2040 TB of memory). This is plenty for all normal
	// use of CBG. If your use-case needs more number of buckets, change the
	// following 2 functions, probably returning 128 bits of hash.
	//
	// [1] 2014 - "Cuckoo Filter Practically Better Than Bloom"
	// by Bin Fan, Dave Andersen, Michael Kaminsky and Michael D. Mitzenmacher
	//
	// Note on performance: The first idea that may occur to you is return a
	// std::pair<size_t, size_t> instead of uint64_t. This is clearer but the
	// performance is 60% worse on lookup when the structure is larger than the
	// cache. Given that a hashtable is used mainly for performance reasons we
	// use the faster version.
	__forceinline uint64_t hash_elem(const KEY_TYPE& elem) const noexcept
	{
		return HASHER::operator()(elem);
	}
	__forceinline uint64_t get_hash1(uint64_t hash) const noexcept
	{
		return rot64(hash, 32);
	}

	/////////////////////////////////////////////////////////////////////
	// Given a value "word", produces an integer in [0,p) without division.
	// The function is as fair as possible in the sense that if you iterate
	// through all possible values of "word", then you will generate all
	// possible outputs as uniformly as possible.
	/////////////////////////////////////////////////////////////////////
	static __forceinline uint32_t fastrange32(uint32_t word, uint32_t p)
	{
		return (uint32_t)(((uint64_t)word * (uint64_t)p) >> 32);
	}
	static __forceinline uint64_t fastrange64(uint64_t word, uint64_t p)
	{
#ifdef __SIZEOF_INT128__ // then we know we have a 128-bit int
		return (uint64_t)(((__uint128_t)word * (__uint128_t)p) >> 64);
#elif defined(_MSC_VER) && defined(_WIN64)
		// supported in Visual Studio 2005 and better
		return __umulh(word, p);
#else
		return word % p; // fallback
#endif // __SIZEOF_INT128__
	}
	static __forceinline size_t fastrange(size_t word, size_t p) {
#if (SIZE_MAX == UINT32_MAX)
		return (size_t)fastrange32(word, p);
#else // assume 64-bit
		return (size_t)fastrange64(word, p);
#endif // SIZE_MAX == UINT32_MAX
	}

	/////////////////////////////////////////////////////////////////////
	// Insertion algorithm utilities
	/////////////////////////////////////////////////////////////////////
	__forceinline void Update_Bin_At_Debug(size_t pos, size_t distance_to_base, bool is_reverse_item, uint_fast16_t label, size_t hash) noexcept
	{
		assert(distance_to_base < NUM_ELEMS_BUCKET);
		assert(label <= L_MAX);
		assert(pos < num_buckets);
		
		METADATA::Update_Bin_At(pos, distance_to_base, is_reverse_item, label, hash);
	}
	std::pair<uint16_t, size_t> Calculate_Minimum(size_t bucket_pos) const noexcept
	{
		uint16_t minimum = METADATA::Get_Label(bucket_pos);
		size_t pos = bucket_pos;

		for (size_t i = 1; minimum && i < NUM_ELEMS_BUCKET; i++)
		{
			uint16_t label_value = METADATA::Get_Label(bucket_pos + i);
			if (minimum > label_value)
			{
				minimum = label_value;
				pos = bucket_pos + i;
			}
		}

		return std::make_pair(minimum, pos);
	}
	__forceinline size_t Belong_to_Bucket(size_t elem_pos) const noexcept
	{
		if (METADATA::Is_Empty(elem_pos))
			return SIZE_MAX;

		return elem_pos + (METADATA::Is_Item_In_Reverse_Bucket(elem_pos) ? NUM_ELEMS_BUCKET - 1 : 0) - METADATA::Distance_to_Entry_Bin(elem_pos);
	}
	size_t Count_Empty(size_t pos) const noexcept
	{
		size_t count = METADATA::Is_Empty(pos) ? 1u : 0u;

		for (size_t i = 1; i < NUM_ELEMS_BUCKET; i++)
			if (METADATA::Is_Empty(pos + i))
				count++;

		return count;
	}
	size_t Count_Elems_In_Bucket_Non_Reversed(size_t bucket_pos) const noexcept
	{
		size_t count = 0;

		for (size_t i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			size_t pos = bucket_pos + i;

			if (!METADATA::Is_Item_In_Reverse_Bucket(pos) && i == METADATA::Distance_to_Entry_Bin(pos))// Faster Belong_to_Bucket(elem_pos)
				count++;
		}

		return count;
	}
	std::pair<size_t, size_t> Count_Elems_In_Bucket_Outside_Range(size_t bucket_pos, size_t range_init) const noexcept
	{
		size_t count_outsize_range = 0;
		size_t count = 0;

		for (size_t i = 0; i < NUM_ELEMS_BUCKET; i++)
		{
			size_t pos = bucket_pos + i;

			if (!METADATA::Is_Item_In_Reverse_Bucket(pos) && i == METADATA::Distance_to_Entry_Bin(pos))// Faster Belong_to_Bucket(elem_pos)
			{
				count++;
				if ((pos - range_init) >= NUM_ELEMS_BUCKET)// Faster check for: !(range_init <= pos < (range_init+NUM_ELEMS_BUCKET)))
					count_outsize_range++;
			}
		}

		return std::make_pair(count, count_outsize_range);
	}
	// TODO: Do this reversing maintaining elems near the entry bin
	void Reverse_Bucket(size_t bucket_pos) noexcept
	{
		METADATA::Set_Bucket_Reversed(bucket_pos);

		size_t j = NUM_ELEMS_BUCKET - 1;
		for (size_t i = NUM_ELEMS_BUCKET - 1; i < NUM_ELEMS_BUCKET; i--)
			if (Belong_to_Bucket(bucket_pos + i) == bucket_pos)// Elems belong to our bucket
			{
				for (; j < NUM_ELEMS_BUCKET && !METADATA::Is_Empty(bucket_pos - j); j--)
				{
				}// Find empty space
				if (j < NUM_ELEMS_BUCKET)
				{
					Update_Bin_At_Debug(bucket_pos - j, NUM_ELEMS_BUCKET - 1 - j, true, METADATA::Get_Label(bucket_pos + i), METADATA::Get_Hash(bucket_pos + i));
					METADATA::Set_Empty(bucket_pos + i);
					DATA::MoveElem(bucket_pos - j, bucket_pos + i);
				}
				else
				{
					assert(i == 0);
					Update_Bin_At_Debug(bucket_pos, NUM_ELEMS_BUCKET - 1, true, METADATA::Get_Label(bucket_pos), METADATA::Get_Hash(bucket_pos));
				}
			}
	}
	// Reverse or Hopscotch for an empty bin. No change if no empty bin es found
	size_t Find_Empty_Pos_Hopscotch(size_t bucket_pos, size_t bucket_init) noexcept
	{
		//////////////////////////////////////////////////////////////////
		// TODO: Try to Re-reverse buckets. Useful to put by default
		// reversed buckets near the end of a cache line. Also needed when
		// erasing elems.
		// TODO: Consider using more sliding windows positions than only
		// normal and reversal
		//////////////////////////////////////////////////////////////////
		// Then try to reverse the bucket
		//////////////////////////////////////////////////////////////////
		if (!METADATA::Is_Bucket_Reversed(bucket_pos) && bucket_pos >= NUM_ELEMS_BUCKET)
		{
			size_t count_empty = Count_Empty(bucket_pos + 1 - NUM_ELEMS_BUCKET);
			if (count_empty)
			{
				size_t count_elems = Count_Elems_In_Bucket_Non_Reversed(bucket_pos);

				// If we can reverse
				if (count_empty > count_elems || (count_empty == count_elems && Belong_to_Bucket(bucket_pos) == bucket_pos))
				{
					if(count_elems)// Some elems
						Reverse_Bucket(bucket_pos);
					else// No elem
						METADATA::Set_Bucket_Reversed(bucket_pos);

					// TODO: Remove this code
					uint16_t min1;
					size_t pos1;
					bucket_init = bucket_pos + (METADATA::Is_Bucket_Reversed(bucket_pos) ? (size_t(1) - NUM_ELEMS_BUCKET) : size_t(0));
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
			for (size_t i = 0; i < NUM_ELEMS_BUCKET; i++)
			{
				size_t pos_elem = bucket_init + i;
				if (!METADATA::Is_Item_In_Reverse_Bucket(pos_elem))
				{
					size_t bucket_elem = pos_elem - METADATA::Distance_to_Entry_Bin(pos_elem);// Faster Belong_to_Bucket(pos_elem)

					if (bucket_elem != bucket_pos)
					{
						size_t count_empty = Count_Empty(bucket_elem + 1 - NUM_ELEMS_BUCKET);// Clearly none of them are inside the bucket range
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

								// TODO: Remove this code
								uint16_t min1;
								size_t pos1;
								std::tie(min1, pos1) = Calculate_Minimum(bucket_init);
								assert(min1 == 0);
								return pos1;
							}
						}
					}
				}
			}

		//////////////////////////////////////////////////////////////////
		// TODO: Try to hopscotch in the other direction
		//////////////////////////////////////////////////////////////////
		// Then try to hopscotch for an empty space
		//////////////////////////////////////////////////////////////////
		size_t max_dist_to_move = NUM_ELEMS_BUCKET - 1;
		for (size_t i = 0; i <= max_dist_to_move && (bucket_init + i) < num_buckets; i++)
		{
			if (METADATA::Is_Empty(bucket_init + i))
			{
				// Find element to move
				size_t pos_blank = bucket_init + i;
				while ((pos_blank - bucket_init) >= NUM_ELEMS_BUCKET)
				{
					size_t pos_swap = pos_blank + 1 - NUM_ELEMS_BUCKET;

					for (; (pos_blank - pos_swap) > (NUM_ELEMS_BUCKET - 1 - METADATA::Distance_to_Entry_Bin(pos_swap)); pos_swap++)
					{
					}// TODO: Use a list with the options to not recalculate again

					// Swap elements
					DATA::MoveElem(pos_blank, pos_swap);
					Update_Bin_At_Debug(pos_blank, METADATA::Distance_to_Entry_Bin(pos_swap) + (pos_blank - pos_swap), METADATA::Is_Item_In_Reverse_Bucket(pos_swap), METADATA::Get_Label(pos_swap), METADATA::Get_Hash(pos_swap));

					pos_blank = pos_swap;
				}

				return pos_blank;
			}
			size_t current_max_move = i + NUM_ELEMS_BUCKET - 1 - METADATA::Distance_to_Entry_Bin(bucket_init + i);
			if (current_max_move > max_dist_to_move)
				max_dist_to_move = current_max_move;
		}

		return SIZE_MAX;// Not found
	}

	// Grow table
	void rehash(size_t new_num_buckets) noexcept
	{
		// Minimum grow space for the in-place algorithm to work
		assert(new_num_buckets >= num_buckets + MIN_BUCKETS_COUNT);

		// TODO: Do the rehash with less additional memory.
		//       Now it uses ~18% when the load_factor is 90%
		// TODO: Consider saving the hash also to speed-up insertion.
		std::vector<INSERT_TYPE> secondary_tmp;
		secondary_tmp.reserve(std::max(1ull, num_elems / 8));// reserve 12.5%
		bool need_rehash = true;

		while (need_rehash)
		{
			need_rehash = false;

			size_t old_num_buckets = num_buckets;
			num_buckets = new_num_buckets;
			new_num_buckets += std::max(static_cast<size_t>(1), new_num_buckets >> 5);// add 3.1% if fails

			// Realloc data
			DATA::ReallocElems(num_buckets);
			METADATA::ReallocMetadata(num_buckets);

			// Initialize metadata
			if (old_num_buckets)
				METADATA::Clear(old_num_buckets, num_buckets - old_num_buckets);
			else
				METADATA::Clear(0, num_buckets);
			num_elems = 0;
			for (size_t i = 0; i < (NUM_ELEMS_BUCKET - 1); i++)
				METADATA::Set_Bucket_Reversed(num_buckets - 1 - i);

			// Moves items from old end to new end
			for (size_t i = old_num_buckets - 1; i > 0; i--)
			{
				if (!METADATA::Is_Empty(i))
				{
					uint64_t hash0 = hash_elem(DATA::GetKey(i));
					size_t bucket1_pos_init = fastrange(hash0, num_buckets);
					bool is_bucket1_reversed = METADATA::Is_Bucket_Reversed(bucket1_pos_init);
					bucket1_pos_init += is_bucket1_reversed ? (1ull - NUM_ELEMS_BUCKET) : 0;
					bool item_is_moved = false;

					// Try to insert primary
					if (bucket1_pos_init > i)
					{
						uint16_t min1;
						size_t pos1;
						std::tie(min1, pos1) = Calculate_Minimum(bucket1_pos_init);
						if (min1 == 0)
						{
							Update_Bin_At_Debug(pos1, pos1 - bucket1_pos_init, is_bucket1_reversed, 1, hash0);
							// Put elem
							DATA::MoveElem(pos1, i);
							num_elems++;
							item_is_moved = true;
						}
					}

					// Not moved -> put in temporary list
					if (!item_is_moved)
						secondary_tmp.push_back(DATA::GetElem(i));
				}
				// Clear position
				METADATA::Clear(i, 1);
			}

			// First element
			if (!METADATA::Is_Empty(0))
				secondary_tmp.push_back(DATA::GetElem(0));
			METADATA::Clear(0, 1);

			// Insert other elements
			while (!secondary_tmp.empty() && !need_rehash)
			{
				if (try_insert(secondary_tmp.back()))
					secondary_tmp.pop_back();
				else
					need_rehash = true;
			}
		}
	}
	size_t get_grow_size() const noexcept
	{
		// Last buckets will be reverted, so they need to be outsize the old buckets
		size_t new_num_buckets = std::max(num_buckets + MIN_BUCKETS_COUNT, size_t(num_buckets*_grow_factor));
		if (new_num_buckets < num_buckets)
			new_num_buckets = SIZE_MAX;

		return new_num_buckets;
	}

	bool try_insert(INSERT_TYPE& elem) noexcept
	{
		while (true)
		{
			// Note on performance: Using only hash0 until needed (Secondary added)
			// will improve preformance by 10-15%. Not enough given that it
			// modify behavior with labels (assuming min2==0 and new label==1).
			uint64_t hash0 = hash_elem(DATA::GetKeyFromValue(elem));
			uint64_t hash1 = get_hash1(hash0);

			// Calculate positions given hash
			size_t bucket1_pos = fastrange(hash0, num_buckets);
			size_t bucket2_pos = fastrange(hash1, num_buckets);

			bool is_reversed_bucket1 = METADATA::Is_Bucket_Reversed(bucket1_pos);
			bool is_reversed_bucket2 = METADATA::Is_Bucket_Reversed(bucket2_pos);
			size_t bucket1_init = bucket1_pos + (is_reversed_bucket1 ? (1ull - NUM_ELEMS_BUCKET) : 0);
			size_t bucket2_init = bucket2_pos + (is_reversed_bucket2 ? (1ull - NUM_ELEMS_BUCKET) : 0);

			// Find minimun label
			uint_fast16_t min1 = METADATA::Get_Label(bucket1_init);
			uint_fast16_t min2 = METADATA::Get_Label(bucket2_init);
			size_t pos1 = bucket1_init;
			size_t pos2 = bucket2_init;
			for (size_t i = 1; i < NUM_ELEMS_BUCKET /*&& (min1 || min2)*/; i++)
			{
				size_t current_pos1 = bucket1_init + i;
				size_t current_pos2 = bucket2_init + i;
				uint_fast16_t label_value1 = METADATA::Get_Label(current_pos1);
				uint_fast16_t label_value2 = METADATA::Get_Label(current_pos2);
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
				Update_Bin_At_Debug(pos1, pos1 - bucket1_init, is_reversed_bucket1, std::min(min2 + 1, L_MAX), hash0);
				// Put elem
				DATA::SaveElem(pos1, elem);
				num_elems++;
				return true;
			}

			size_t empty_pos = Find_Empty_Pos_Hopscotch(bucket1_pos, bucket1_init);
			if (empty_pos != SIZE_MAX)
			{
				is_reversed_bucket1 = METADATA::Is_Bucket_Reversed(bucket1_pos);
				bucket1_init = bucket1_pos + (is_reversed_bucket1 ? (1 - NUM_ELEMS_BUCKET) : 0);
				Update_Bin_At_Debug(empty_pos, empty_pos - bucket1_init, is_reversed_bucket1, std::min(min2 + 1, L_MAX), hash0);

				// Put elem
				DATA::SaveElem(empty_pos, elem);
				num_elems++;
				return true;
			}

			///////////////////////////////////////////////////////////////////
			// Secondary added, Unlucky bucket added
			//////////////////////////////////////////////////////////////////
			if (min2 == 0)
			{
				METADATA::Set_Unlucky_Bucket(bucket1_pos);
				Update_Bin_At_Debug(pos2, pos2 - bucket2_init, is_reversed_bucket2, std::min(min1 + 1, L_MAX), hash1);
				// Put elem
				DATA::SaveElem(pos2, elem);
				num_elems++;
				return true;
			}

			//if (num_elems * 10 > 9 * num_buckets)// > 90%
			{
				empty_pos = Find_Empty_Pos_Hopscotch(bucket2_pos, bucket2_init);

				if (empty_pos != SIZE_MAX)
				{
					METADATA::Set_Unlucky_Bucket(bucket1_pos);
					is_reversed_bucket2 = METADATA::Is_Bucket_Reversed(bucket2_pos);
					bucket2_init = bucket2_pos + (is_reversed_bucket2 ? (1 - NUM_ELEMS_BUCKET) : 0);
					Update_Bin_At_Debug(empty_pos, empty_pos - bucket2_init, is_reversed_bucket2, std::min(min1 + 1, L_MAX), hash1);

					// Put elem
					DATA::SaveElem(empty_pos, elem);
					num_elems++;
					return true;
				}
			}

			// Terminating condition
			if (std::min(min1, min2) >= L_MAX)
				return false;

			if (min1 <= min2)// Selected pos in first bucket
			{
				Update_Bin_At_Debug(pos1, pos1 - bucket1_init, is_reversed_bucket1, std::min(min2 + 1, L_MAX), hash0);
				// Put elem
				INSERT_TYPE victim = DATA::GetElem(pos1);
				DATA::SaveElem(pos1, elem);
				elem = victim;
			}
			else
			{
				METADATA::Set_Unlucky_Bucket(bucket1_pos);
				Update_Bin_At_Debug(pos2, pos2 - bucket2_init, is_reversed_bucket2, std::min(min1 + 1, L_MAX), hash1);
				// Put elem
				INSERT_TYPE victim = DATA::GetElem(pos2);
				DATA::SaveElem(pos2, elem);
				elem = victim;
			}
		}
	}

	///////////////////////////////////////////////////////////////////////////////
	// Find an element
	///////////////////////////////////////////////////////////////////////////////
	size_t find_position_negative(const KEY_TYPE& elem) const noexcept
	{
		uint64_t hash = hash_elem(elem);

		// Check first bucket
		size_t pos = fastrange(hash, num_buckets);

		uint_fast16_t c0 = METADATA::at(pos);

		uint_fast16_t h = static_cast<uint_fast16_t>(hash);
		if (((c0 ^ h) & 0xFF00) == 0 && (c0 & 0b111) && cmp_elems(pos, elem))
			return pos;

		if (c0 & 0b01'000'000)// Is_Reversed_Window(pos)
		{
			uint_fast16_t cc = METADATA::at(pos-1);
			if (((cc ^ h) & 0xFF00) == 0 && (cc & 0b111) && cmp_elems(pos-1, elem))
				return pos-1;

			if (NUM_ELEMS_BUCKET > 2)
			{
				cc = METADATA::at(pos-2);
				if (((cc ^ h) & 0xFF00) == 0 && (cc & 0b111) && cmp_elems(pos-2, elem))
					return pos-2;
			}
			if (NUM_ELEMS_BUCKET > 3)
			{
				cc = METADATA::at(pos-3);
				if (((cc ^ h) & 0xFF00) == 0 && (cc & 0b111) && cmp_elems(pos-3, elem))
					return pos-3;
			}
		}
		else// Normal
		{
			uint_fast16_t cc = METADATA::at(pos+1);
			if (((cc ^ h) & 0xFF00) == 0 && (cc & 0b111) && cmp_elems(pos+1, elem))
				return pos+1;

			if (NUM_ELEMS_BUCKET > 2)
			{
				cc = METADATA::at(pos+2);
				if (((cc ^ h) & 0xFF00) == 0 && (cc & 0b111) && cmp_elems(pos+2, elem))
					return pos+2;
			}
			if (NUM_ELEMS_BUCKET > 3)
			{
				cc = METADATA::at(pos+3);
				if (((cc ^ h) & 0xFF00) == 0 && (cc & 0b111) && cmp_elems(pos+3, elem))
					return pos+3;
			}
		}

		// Check second bucket
		if (c0 & 0b10'000'000)// Is_Unlucky_Bucket(pos)
		{
			hash = get_hash1(hash);
			pos = fastrange(hash, num_buckets);

			uint_fast16_t cc = METADATA::at(pos);

			h = static_cast<uint_fast16_t>(hash);
			if (((cc ^ h) & 0xFF00) == 0 && (cc & 0b111) && cmp_elems(pos, elem))
				return pos;

			size_t reverse_sum = /*Is_Reversed_Window(pos)*/cc & 0b01'000'000 ? static_cast<size_t>(-1) : static_cast<size_t>(1);

			pos += reverse_sum;
			cc = METADATA::at(pos);
			if (((cc ^ h) & 0xFF00) == 0 && (cc & 0b111) && cmp_elems(pos, elem))
				return pos;
			if (NUM_ELEMS_BUCKET > 2)
			{
				pos += reverse_sum;
				cc = METADATA::at(pos);
				if (((cc ^ h) & 0xFF00) == 0 && (cc & 0b111) && cmp_elems(pos, elem))
					return pos;
			}
			if (NUM_ELEMS_BUCKET > 3)
			{
				pos += reverse_sum;
				cc = METADATA::at(pos);
				if (((cc ^ h) & 0xFF00) == 0 && (cc & 0b111) && cmp_elems(pos, elem))
					return pos;
			}
		}

		return SIZE_MAX;
	}
	size_t find_position_negative_no_metadata(const KEY_TYPE& elem) const noexcept
	{
		uint64_t hash = hash_elem(elem);

		// Check first bucket
		size_t pos = fastrange(hash, num_buckets);

		uint_fast16_t c0 = METADATA::at(pos);

		if (cmp_elems(pos, elem) && (c0 & 0b111))
			return pos;

		if (c0 & 0b01'000'000)// Is_Reversed_Window(pos)
		{
			uint_fast16_t cc = METADATA::at(pos-1);
			if (cmp_elems(pos-1, elem) && (cc & 0b111))
				return pos-1;

			if (NUM_ELEMS_BUCKET > 2)
			{
				cc = METADATA::at(pos-2);
				if (cmp_elems(pos-2, elem) && (cc & 0b111))
					return pos-2;
			}
			if (NUM_ELEMS_BUCKET > 3)
			{
				cc = METADATA::at(pos-3);
				if (cmp_elems(pos-3, elem) && (cc & 0b111))
					return pos-3;
			}
		}
		else// Normal
		{
			uint_fast16_t cc = METADATA::at(pos+1);
			if (cmp_elems(pos+1, elem) && (cc & 0b111))
				return pos+1;

			if (NUM_ELEMS_BUCKET > 2)
			{
				cc = METADATA::at(pos+2);
				if (cmp_elems(pos+2, elem) && (cc & 0b111))
					return pos+2;
			}
			if (NUM_ELEMS_BUCKET > 3)
			{
				cc = METADATA::at(pos+3);
				if (cmp_elems(pos+3, elem) && (cc & 0b111))
					return pos+3;
			}
		}

		// Check second bucket
		if (c0 & 0b10'000'000)//Is_Unlucky_Bucket(pos)
		{
			hash = get_hash1(hash);
			pos = fastrange(hash, num_buckets);

			uint_fast16_t cc = METADATA::at(pos);

			if (cmp_elems(pos, elem) && (cc & 0b111))
				return pos;

			size_t reverse_sum = /*Is_Reversed_Window(pos)*/cc & 0b01'000'000 ? static_cast<size_t>(-1) : static_cast<size_t>(1);

			pos += reverse_sum;
			cc = METADATA::at(pos);
			if (cmp_elems(pos, elem) && (cc & 0b111))
				return pos;
			if (NUM_ELEMS_BUCKET > 2)
			{
				pos += reverse_sum;
				cc = METADATA::at(pos);
				if (cmp_elems(pos, elem) && (cc & 0b111))
					return pos;
			}
			if (NUM_ELEMS_BUCKET > 3)
			{
				pos += reverse_sum;
				cc = METADATA::at(pos);
				if (cmp_elems(pos, elem) && (cc & 0b111))
					return pos;
			}
		}

		return SIZE_MAX;
	}
	size_t find_position_positive(const KEY_TYPE& elem) const noexcept
	{
		uint64_t hash = hash_elem(elem);

		// Check first bucket
		size_t pos = fastrange(hash, num_buckets);

		uint_fast16_t c0 = METADATA::at(pos);

		if (cmp_elems(pos, elem) && (c0 & 0b111))
			return pos;

		if (c0 & 0b01'000'000)// Is_Reversed_Window(pos)
		{
			uint_fast16_t cc = METADATA::at(pos - 1);
			if (cmp_elems(pos - 1, elem) && (cc & 0b111))
				return pos - 1;

			if (NUM_ELEMS_BUCKET > 2)
			{
				cc = METADATA::at(pos - 2);
				if (cmp_elems(pos - 2, elem) && (cc & 0b111))
					return pos - 2;
			}
			if (NUM_ELEMS_BUCKET > 3)
			{
				cc = METADATA::at(pos - 3);
				if (cmp_elems(pos - 3, elem) && (cc & 0b111))
					return pos - 3;
			}
		}
		else// Normal
		{
			uint_fast16_t cc = METADATA::at(pos + 1);
			if (cmp_elems(pos + 1, elem) && (cc & 0b111))
				return pos + 1;

			if (NUM_ELEMS_BUCKET > 2)
			{
				cc = METADATA::at(pos + 2);
				if (cmp_elems(pos + 2, elem) && (cc & 0b111))
					return pos + 2;
			}
			if (NUM_ELEMS_BUCKET > 3)
			{
				cc = METADATA::at(pos + 3);
				if (cmp_elems(pos + 3, elem) && (cc & 0b111))
					return pos + 3;
			}
		}

		// Check second bucket
		//if (c0 & 0b10'000'000)//Is_Unlucky_Bucket(pos)
		{
			hash = get_hash1(hash);
			pos = fastrange(hash, num_buckets);

			uint_fast16_t cc = METADATA::at(pos);

			if (cmp_elems(pos, elem) && (cc & 0b111))
				return pos;

			size_t reverse_sum = /*Is_Reversed_Window(pos)*/cc & 0b01'000'000 ? static_cast<size_t>(-1) : static_cast<size_t>(1);

			pos += reverse_sum;
			cc = METADATA::at(pos);
			if (cmp_elems(pos, elem) && (cc & 0b111))
				return pos;
			if (NUM_ELEMS_BUCKET > 2)
			{
				pos += reverse_sum;
				cc = METADATA::at(pos);
				if (cmp_elems(pos, elem) && (cc & 0b111))
					return pos;
			}
			if (NUM_ELEMS_BUCKET > 3)
			{
				pos += reverse_sum;
				cc = METADATA::at(pos);
				if (cmp_elems(pos, elem) && (cc & 0b111))
					return pos;
			}
		}

		return SIZE_MAX;
	}
	__forceinline size_t find_position(const KEY_TYPE& elem, Search_Hint hint) const noexcept
	{
		if (hint == Search_Hint::Expect_Positive)
			return find_position_positive(elem);// Positive queries prefered

		// Negative or unknow queries
		if (USE_METADATA_SEARCH)
			return find_position_negative(elem);
		else
			return find_position_negative_no_metadata(elem);
	}

public:
	///////////////////////////////////////////////////////////////////////////////
	// Constructors
	///////////////////////////////////////////////////////////////////////////////
	// Default
	CBG_IMPL() noexcept : num_elems(0), num_buckets(0), HASHER(), EQ(), DATA()
	{}
	explicit CBG_IMPL(size_t expected_num_elems) noexcept : HASHER(), EQ(), DATA(std::max(MIN_BUCKETS_COUNT, expected_num_elems)),
		num_elems(0), num_buckets(std::max(MIN_BUCKETS_COUNT, expected_num_elems))
	{
		// TODO: Call clear() here instead of this code?
		// Last buckets are always reversed
		for (size_t i = 0; i < (NUM_ELEMS_BUCKET - 1); i++)
			METADATA::Set_Bucket_Reversed(num_buckets - 1 - i);

		// TODO: Check why this does not improve performance
		// This need to be repetead in clear() and rehash() methods
		/*for (size_t i = NUM_ELEMS_BUCKET - 1; i < num_buckets; i++)
		if(METADATA::Benefit_With_Reversal(i, NUM_ELEMS_BUCKET))
		METADATA::Set_Bucket_Reversed(i);*/
	}
	// TODO: Do all this constructors
	// Initializer List constructor
	//CBG_IMPL(std::initializer_list<INSERT_TYPE> init) noexcept : CBG_IMPL(init.size() / _max_load_factor)
	//{
	//	// TODO...
	//}
	// Range constructor
	//CBG_IMPL(InputIter first, InputIter last) noexcept : CBG_IMPL(init.size() / _max_load_factor)
	//{
	//	// TODO...
	//}
	// Copy constructor
	CBG_IMPL(const CBG_IMPL& that) = delete;
	// Move constructor
	CBG_IMPL(const CBG_IMPL&& that) = delete;
	// Copy assignment operator
	CBG_IMPL& operator=(const CBG_IMPL& that) = delete;
	// Move assignment operator
	CBG_IMPL& operator=(CBG_IMPL&& that) = delete;

	~CBG_IMPL() noexcept
	{
		num_elems = 0;
		num_buckets = 0;
	}
	///////////////////////////////////////////////////////////////////////////////

	size_t capacity() const noexcept
	{
		return num_buckets;
	}
	// Same as capacity()
	size_t bucket_count() const noexcept
	{
		return num_buckets;
	}
	size_t size() const noexcept
	{
		return num_elems;
	}
	bool empty() const noexcept
	{
		return num_elems == 0;
	}
	void clear() noexcept
	{
		num_elems = 0;
		METADATA::Clear(0, num_buckets);

		// Last buckets are always reversed
		for (size_t i = 0; i < (NUM_ELEMS_BUCKET - 1); i++)
			METADATA::Set_Bucket_Reversed(num_buckets - 1 - i);

		/*for (size_t i = NUM_ELEMS_BUCKET - 1; i < num_buckets; i++)
			if (METADATA::Benefit_With_Reversal(i, NUM_ELEMS_BUCKET))
				METADATA::Set_Bucket_Reversed(i);*/
	}
	float load_factor() const noexcept
	{
		return size() * 100f / capacity();
	}
	void max_load_factor(float value) noexcept
	{
		_max_load_factor = value;
	}
	float max_load_factor() const noexcept
	{
		return _max_load_factor;
	}
	void grow_factor(float value) noexcept
	{
		_grow_factor = value;
	}
	float grow_factor() const noexcept
	{
		return _grow_factor;
	}

	void reserve(size_t new_capacity) noexcept
	{
		if(new_capacity >= num_buckets + MIN_BUCKETS_COUNT)
			rehash(new_capacity);
	}

	void insert(const INSERT_TYPE& to_insert_elem) noexcept
	{
		if (num_elems >= num_buckets * _max_load_factor)
			rehash(get_grow_size());

		INSERT_TYPE elem = to_insert_elem;
		// TODO: break infinity cycle when -> num_buckets=SIZE_MAX
		while (!try_insert(elem))
			rehash(get_grow_size());
	}

	// Check if an element exist
	uint32_t count(const KEY_TYPE& elem, Search_Hint hint = Search_Hint::Unknow) const noexcept
	{
		return find_position(elem, hint) != SIZE_MAX ? 1u : 0u;
	}
	bool contains(const KEY_TYPE& elem, Search_Hint hint = Search_Hint::Unknow) const noexcept
	{
		return find_position(elem, hint) != SIZE_MAX;
	}

	// TODO: As currently implemented the performance may
	//       degrade if many erase operations are done.
	uint32_t erase(const KEY_TYPE& elem) noexcept
	{
		size_t elem_pos = find_position(elem);
		if (elem_pos != SIZE_MAX)
		{
			METADATA::Set_Empty(elem_pos);
			num_elems--;
			return 1;
		}

		return 0;
	}
};

// Map. Only added simple mapping operations.
template<size_t NUM_ELEMS_BUCKET, class KEY, class T, class HASHER, class EQ, class DATA, class METADATA, bool IS_NEGATIVE> class CBG_MAP_IMPL : 
	public CBG_IMPL<NUM_ELEMS_BUCKET, std::pair<KEY, T>, KEY, T, HASHER, EQ, DATA, METADATA, IS_NEGATIVE>
{
	using Base = typename CBG_MAP_IMPL::CBG_IMPL;

public:
	using Base::Base;

	// Map operations
	T& operator[](const KEY& key) noexcept
	{
		size_t key_pos = find_position(key, Search_Hint::Expect_Positive);
		if (key_pos == SIZE_MAX)
		{
			insert(std::make_pair(key, T()));
			key_pos = find_position(key);
		}

		return *DATA::GetValue(key_pos);
	}
	T& operator[](KEY&& key) noexcept
	{
		size_t key_pos = find_position(key, Search_Hint::Expect_Positive);
		if (key_pos == SIZE_MAX)
		{
			insert(std::make_pair(std::move(key), T()));
			key_pos = find_position(key);
		}

		return *DATA::GetValue(key_pos);
	}
	T& at(const KEY& key)
	{
		size_t key_pos = find_position(key, Search_Hint::Expect_Positive);
		if (key_pos == SIZE_MAX)
			throw std::out_of_range("Argument passed to at() was not in the map.");

		return *DATA::GetValue(key_pos);
	}
	const T& at(const KEY& key) const
	{
		size_t key_pos = find_position(key, Search_Hint::Expect_Positive);
		if (key_pos == SIZE_MAX)
			throw std::out_of_range("Argument passed to at() was not in the map.");

		return *DATA::GetValue(key_pos);
	}
};
}// end namespace cbg_internal

///////////////////////////////////////////////////////////////////////////////
// CBG Sets
///////////////////////////////////////////////////////////////////////////////
// (Struct of Arrays)
template<size_t NUM_ELEMS_BUCKET, class T, class HASHER = hashing::t1ha2<T>, class EQ = std::equal_to<T>> class Set_SoA :
	public cbg_internal::CBG_IMPL<NUM_ELEMS_BUCKET, T, T, T, HASHER, EQ, cbg_internal::KeyLayout_SoA<T>, cbg_internal::MetadataLayout_SoA, true>
{
	using Base = typename Set_SoA::CBG_IMPL;

public:
	using Base::Base;
};
// (Array of structs)
template<size_t NUM_ELEMS_BUCKET, class T, class HASHER = hashing::t1ha2<T>, class EQ = std::equal_to<T>> class Set_AoS :
	public cbg_internal::CBG_IMPL<NUM_ELEMS_BUCKET, T, T, T, HASHER, EQ, cbg_internal::KeyLayout_AoS<T>, cbg_internal::MetadataLayout_AoS<sizeof(T)>, false>
{
	using Base = typename Set_AoS::CBG_IMPL;

public:
	using Base::Base;
};
// (Array of blocks)
template<size_t NUM_ELEMS_BUCKET, class T, class HASHER = hashing::t1ha2<T>, class EQ = std::equal_to<T>> class Set_AoB :
	public cbg_internal::CBG_IMPL<NUM_ELEMS_BUCKET, T, T, T, HASHER, EQ, cbg_internal::KeyLayout_AoB<T>, cbg_internal::MetadataLayout_AoB<alignof(T), cbg_internal::BlockKey<T>>, false>
{
	using Base = typename Set_AoB::CBG_IMPL;

public:
	using Base::Base;
};
///////////////////////////////////////////////////////////////////////////////
// CBG Maps
///////////////////////////////////////////////////////////////////////////////
// (Struct of Arrays)
template<size_t NUM_ELEMS_BUCKET, class KEY, class T, class HASHER = hashing::t1ha2<KEY>, class EQ = std::equal_to<KEY>> class Map_SoA :
	public cbg_internal::CBG_MAP_IMPL<NUM_ELEMS_BUCKET, KEY, T, HASHER, EQ, cbg_internal::MapLayout_SoA<KEY, T>, cbg_internal::MetadataLayout_SoA, true>
{
	using Base = typename Map_SoA::CBG_MAP_IMPL;

public:
	using Base::Base;
};
// (Array of structs)
template<size_t NUM_ELEMS_BUCKET, class KEY, class T, class HASHER = hashing::t1ha2<KEY>, class EQ = std::equal_to<KEY>> class Map_AoS :
	public cbg_internal::CBG_MAP_IMPL<NUM_ELEMS_BUCKET, KEY, T, HASHER, EQ, cbg_internal::MapLayout_AoS<KEY, T>, cbg_internal::MetadataLayout_AoS<sizeof(KEY) + sizeof(T)>, false>
{
	using Base = typename Map_AoS::CBG_MAP_IMPL;

public:
	using Base::Base;
};
// (Array of blocks)
template<size_t NUM_ELEMS_BUCKET, class KEY, class T, class HASHER = hashing::t1ha2<KEY>, class EQ = std::equal_to<KEY>> class Map_AoB :
	public cbg_internal::CBG_MAP_IMPL<NUM_ELEMS_BUCKET, KEY, T, HASHER, EQ, cbg_internal::MapLayout_AoB<KEY, T>, cbg_internal::MetadataLayout_AoB<cbg_internal::MaxAlignOf<KEY, T>::BLOCK_SIZE, cbg_internal::BlockMap<KEY, T>>, false>
{
	using Base = typename Map_AoB::CBG_MAP_IMPL;

public:
	using Base::Base;
};
}// end namespace cbg
