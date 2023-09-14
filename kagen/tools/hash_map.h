#pragma once

#ifdef KAGEN_SPARSEHASH_FOUND
    #include <google/dense_hash_map>
#else // KAGEN_SPARSEHASH_FOUND
    #include <unordered_map>
#endif // KAGEN_SPARSEHASH_FOUND

namespace kagen {
#ifdef KAGEN_SPARSEHASH_FOUND
template <typename Key, typename Value>
using HashMap = google::dense_hash_map<Key, Value>;
#else
template <typename Key, typename Value>
struct HashMap {
    using Underlaying = std::unordered_map<Key, Value>;

    using key_type             = typename Underlaying::key_type;
    using mapped_type          = typename Underlaying::mapped_type;
    using value_type           = typename Underlaying::value_type;
    using size_type            = typename Underlaying::size_type;
    using difference_type      = typename Underlaying::difference_type;
    using hasher               = typename Underlaying::hasher;
    using key_equal            = typename Underlaying::key_equal;
    using allocator_type       = typename Underlaying::allocator_type;
    using reference            = typename Underlaying::reference;
    using const_reference      = typename Underlaying::const_reference;
    using iterator             = typename Underlaying::iterator;
    using const_iterator       = typename Underlaying::const_iterator;
    using local_iterator       = typename Underlaying::local_iterator;
    using const_local_iterator = typename Underlaying::const_local_iterator;

    //
    // Additional functions
    //

    void set_empty_key(Key) {}

    //
    // Iterators
    //

    inline iterator begin() noexcept {
        return map.begin();
    }

    inline iterator end() noexcept {
        return map.end();
    }

    inline const_iterator begin() const noexcept {
        return map.begin();
    }

    inline const_iterator end() const noexcept {
        return map.end();
    }

    inline const_iterator cbegin() const noexcept {
        return map.cbegin();
    }

    inline const_iterator cend() const noexcept {
        return map.cend();
    }

    //
    // Capacity
    //

    inline bool empty() const noexcept {
        return map.empty();
    }

    inline size_type size() const noexcept {
        return map.size();
    }

    inline size_type max_size() const noexcept {
        return map.max_size();
    }

    //
    // Modifiers
    //

    inline void clear() noexcept {
        return map.clear();
    }

    //
    // Lookup
    //

    inline Value& at(const Key& key) {
        return map.at(key);
    }

    inline const Value& at(const Key& key) const {
        return map.at(key);
    }

    inline Value& operator[](const Key& key) {
        return map[key];
    }

    inline Value& operator[](Key&& key) {
        return map[std::forward<Key>(key)];
    }

    inline size_type count(const Key& key) const {
        return map.count(key);
    }

    template <typename K>
    inline size_type count(const K& x) const {
        return map.count(x);
    }

    iterator find(const Key& key) {
        return map.find(key);
    }

    const_iterator find(const Key& key) const {
        return map.find(key);
    }

    template <typename K>
    iterator find(const K& x) {
        return map.find(x);
    }

    template <typename K>
    const_iterator find(const K& x) const {
        return map.find(x);
    }

    Underlaying map;
};
#endif
} // namespace kagen
