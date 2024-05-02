#include <cstddef>

#include <algorithm>
#include <atomic>
#include <optional>
#include <utility>

#include <parlay/primitives.h>
#include <parlay/sequence.h>

// **************************************************************
// A simple concurrent hash-based map
// Supports concurrent linearizable, insert, find and remove.
// size(), keys() don't linearize with updates.
// Requires the capacity to be  specified on construction.
// No more than the capacity distinct keys can ever be added.
// Once a key is added, removing it will empty the value and mark
// the key as deleted, but only a value with the same key can use the
// same slot (i.e. it still counts towards the capacity).
// It uses locks, but holds them very briefly.
// **************************************************************

template <typename K,
	  typename V,
	  typename Hash = parlay::hash<K>,
	  typename Equal = std::equal_to<>>
struct convex_hash_map {
 private:
  using KV = std::pair<K,V>;
  using index = unsigned long;
  long m;
  index start_index(K k) { return hash(k) % m;}
  index next_index(index h) { return (h + 1 == m) ? 0 : h + 1; }
  Hash hash;
  Equal equal;

  struct entry {
    std::atomic<bool> taken;
    std::atomic<bool> check;
    std::atomic<bool> removed;
    K key;
    V value;
    entry() : taken(false), check(false), removed(false) {}
  };
  parlay::sequence<entry> R;

 public:
  convex_hash_map(long size, Hash&& hash = {}, Equal&& equal = {}) 
    : m(100 + static_cast<index>(1.5 * size)),
      hash(hash), equal(equal),
      R(parlay::sequence<entry>(m)) {}

  bool insert_and_set(const K& k, const V& v) {
    index start = start_index(k);
    index i = start_index(k);
    while (true) {
      bool expected = false;
      if (R[i].taken.compare_exchange_strong(expected, true)) {
        break;
      }
      i = next_index(i);
      if (i == start) {
        std::cout << "Hash table overfull" << std::endl;
        return false;
      }
    }
    R[i].key = k;
    R[i].value = v;

    i = start;
    while (R[i].taken) {
      if (R[i].key == k) {
        bool expected_check = false;
        if (false == R[i].check.compare_exchange_strong(expected_check, true)) {
          return false;
        }
      }
      i = next_index(i);
    }
    return true;
  }

  // used for convex hull
  bool remove(const K& k) {
    index i = start_index(k);
    while (true) {
      if (R[i].taken == false) return false;
      if (R[i].removed == false && R[i].key == k) {
        R[i].removed = true;
        return true;
      }
      i = next_index(i);
    }
  }

  std::optional<V> get_value(const K& k, const V& v) {
    index i = start_index(k);
    while (R[i].taken) {
      if (R[i].key == k) {
        if (R[i].value != v) {
          return R[i].value;
        }
      }
      i = next_index(i);
    }

    return {};
  }

  parlay::sequence<K> keys() {
    return parlay::map_maybe(R, [] (const entry &x) {
	return (x.taken && !x.removed) ? std::optional{x.key} : std::optional<K>{};});
  }
};

