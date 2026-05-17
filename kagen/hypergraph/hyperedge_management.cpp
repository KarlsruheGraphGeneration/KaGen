
#include "kagen/kagen.h"

namespace kagen {

std::vector<PinRange> RemoveEmpty(std::vector<PinRange> hyperedge_ranges) {
    hyperedge_ranges.erase(
        std::remove_if(
            hyperedge_ranges.begin(), hyperedge_ranges.end(), [](const PinRange& r) { return r.begin >= r.end; }),
        hyperedge_ranges.end());
    return hyperedge_ranges;
}

std::vector<PinRange> SortRanges(std::vector<PinRange> hyperedge_ranges) {
    std::sort(hyperedge_ranges.begin(), hyperedge_ranges.end(), [](const PinRange& a, const PinRange& b) {
        if (a.begin != b.begin) {
            return a.begin < b.begin;
        }
        return a.end < b.end;
    });

    return hyperedge_ranges;
}

std::vector<PinRange> MergeAdjacentRanges(std::vector<PinRange> hyperedge_ranges) {
    std::vector<PinRange> merged_ranges;
    for (const PinRange& range: hyperedge_ranges) {
        if (merged_ranges.empty() || merged_ranges.back().end < range.begin) {
            merged_ranges.push_back(range);
        } else {
            merged_ranges.back().end = std::max(merged_ranges.back().end, range.end);
        }
    }
    hyperedge_ranges = std::move(merged_ranges);

    return hyperedge_ranges;
}

std::vector<SInt> SortDeduplicatePins(std::vector<SInt> hyperedge_pins) {
    std::sort(hyperedge_pins.begin(), hyperedge_pins.end());
    hyperedge_pins.erase(std::unique(hyperedge_pins.begin(), hyperedge_pins.end()), hyperedge_pins.end());

    return hyperedge_pins;
}

std::vector<SInt> RemovePinsInRanges(std::vector<SInt> hyperedge_pins, const std::vector<PinRange>& hyperedge_ranges) {
    auto covered_by_range = [&hyperedge_ranges](const SInt pin) {
        auto it = std::upper_bound(
            hyperedge_ranges.begin(), hyperedge_ranges.end(), pin,
            [](const SInt value, const PinRange& range) { return value < range.begin; });

        if (it == hyperedge_ranges.begin()) {
            return false;
        }

        --it;
        return it->begin <= pin && pin < it->end;
    };

    hyperedge_pins.erase(
        std::remove_if(hyperedge_pins.begin(), hyperedge_pins.end(), covered_by_range), hyperedge_pins.end());

    return hyperedge_pins;
}

std::pair<std::vector<SInt>, std::vector<PinRange>> ConvertConsecutivePinsToRanges(
    const std::vector<SInt>& pins, std::vector<PinRange> ranges, const SInt min_run_length = 4) {
    std::vector<SInt> remaining_pins;

    for (SInt i = 0; i < static_cast<SInt>(pins.size());) {
        const SInt begin = pins[i];
        SInt       end   = begin + 1;

        SInt j = i + 1;
        while (j < static_cast<SInt>(pins.size()) && pins[j] == end) {
            ++end;
            ++j;
        }

        const SInt run_length = end - begin;

        if (run_length >= min_run_length) {
            ranges.push_back({begin, end});
        } else {
            for (SInt v = begin; v < end; ++v) {
                remaining_pins.push_back(v);
            }
        }

        i = j;
    }

    return {remaining_pins, ranges};
}

std::pair<std::vector<SInt>, std::vector<PinRange>>
NormalizeCurrentHyperedge(std::vector<SInt> pins, std::vector<PinRange> ranges, const SInt min_run_length = 4) {
    ranges = RemoveEmpty(std::move(ranges));

    ranges = SortRanges(std::move(ranges));

    ranges = MergeAdjacentRanges(ranges);

    pins = SortDeduplicatePins(std::move(pins));

    pins = RemovePinsInRanges(std::move(pins), ranges);

    auto converted = ConvertConsecutivePinsToRanges(pins, std::move(ranges), min_run_length);
    pins           = std::move(converted.first);
    ranges         = std::move(converted.second);

    ranges = RemoveEmpty(std::move(ranges));
    ranges = SortRanges(std::move(ranges));
    ranges = MergeAdjacentRanges(ranges);

    pins = RemovePinsInRanges(std::move(pins), ranges);

    return {pins, ranges};
}

} // namespace kagen