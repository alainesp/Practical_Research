## Preamble

Practical research on computer science is incredible important to the advance of our field but with scientists generating enormous amounts of research/papers is very easy for a *state-of-the-art* algorithm to remain under the radar for years, even decades.

The reason of existence of this repository is to list the *state-of-the-art* algorithms/libraries/resources for common tasks in the programming field.

## Hash Tables

There are many hash tables types out there, each one with different characteristics. Hash tables provided by libraries in programming languages are very conservatives and have mediocre performance. it is almost always better to use an external library.

**Recommended:**

- **Cuckoo Hash Table**: High table load (until **99%**) with deterministic worst case access time. See [Cuckoo Breeding Ground](research_cuckoo_cbg.md) for our variant of it.

### Non-cryptographic hashing

Hash functions typically use between 32-bit and 64 or 128-bit hashes. Typical median key size in perl5 is **20**, the most common **4**. Similar for all other dynamic languages.

**Recommended:**

- **t1ha2**: "Fast Positive Hash" for 64-bit CPUs.
- [SMHasher](https://github.com/rurban/smhasher): The most well known test suite designed to test the distribution, collision, and performance properties of non-cryptographic hash functions.
- [A Seven-Dimensional Analysis of Hashing Methods and its Implications on Query Processing](https://infosys.cs.uni-saarland.de/publications/p249-richter.pdf) for a concise overview of the best hash table strategies, confirming that the simplest Mult hashing (bernstein, FNV*, x17, sdbm) always beat "better" hash functions (Tabulation, Murmur, Farm, ...) when used in a hash table.

The hash table attacks described in [SipHash](https://131002.net/siphash/) against City, Murmur or Perl JenkinsOAAT can be safely ignored in most cases. If you insist on the use of a *secure* hash function consider [HighwayHash](https://arxiv.org/abs/1612.06257).
