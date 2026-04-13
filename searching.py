from pathlib import Path
import json
import time

import matplotlib.pyplot as plt

from generators import unordered_sequence, ordered_sequence, dna_sequence


ALLOWED_FIELDS = ("unordered_numbers", "ordered_numbers", "dna_sequence")


def read_data(file_name, field):
    """
    Reads a JSON file and returns data for a given field.

    Args:
        file_name (str): Name of the JSON file.
        field (str): Key to retrieve from the JSON data.
            Must be one of: 'unordered_numbers', 'ordered_numbers' or 'dna_sequence'.

    Returns:
        list | str | None:
            - list: If data retrieved by the selected field contains numeric data.
            - str: If field is 'dna_sequence'.
            - None: If the field is not supported.
    """
    if field not in ALLOWED_FIELDS:
        return None

    file_path = Path.cwd() / file_name
    with file_path.open(encoding="utf-8") as f:
        data = json.load(f)

    return data.get(field)


def linear_search(sequence, target_number):
    """
    Sequentially searches for target_number in sequence.

    Args:
        sequence (list): Sequence to search in.
        target_number: Value to find.

    Returns:
        dict: {'positions': list[int], 'count': int}
    """
    positions = []
    for index, value in enumerate(sequence):
        if value == target_number:
            positions.append(index)
    return {"positions": positions, "count": len(positions)}


def binary_search(sorted_list, target_number):
    """
    Binary search in a sorted list.

    Args:
        sorted_list (list): Sorted sequence.
        target_number: Value to find.

    Returns:
        int | None: Index of target_number if found, otherwise None.
    """
    low = 0
    high = len(sorted_list) - 1
    while low <= high:
        mid = (low + high) // 2
        middle_value = sorted_list[mid]
        if middle_value == target_number:
            return mid
        if middle_value < target_number:
            low = mid + 1
        else:
            high = mid - 1
    return None


def pattern_search(sequence, pattern):
    """
    Naive pattern search. Returns set of starting indices where pattern occurs.

    Uses early exit on first mismatch inside each window.

    Args:
        sequence (str): Searched sequence.
        pattern (str): Pattern to find.

    Returns:
        set[int]: Set of starting indices where the pattern occurs.
    """
    positions = set()
    n = len(sequence)
    m = len(pattern)
    if m == 0 or m > n:
        return positions

    for i in range(n - m + 1):
        match = True
        for j in range(m):
            if sequence[i + j] != pattern[j]:
                match = False
                break
        if match:
            positions.add(i)
    return positions


def measure_time(func, *args):
    start = time.perf_counter()
    func(*args)
    end = time.perf_counter()
    return end - start


def compare_search_times(sizes=(100, 500, 1000, 5000, 10000, 50000)):
    """
    Compares linear, binary and set membership lookup on sequences of various sizes.
    Returns a dict of measurements and shows a matplotlib plot.
    """
    linear_times = []
    binary_times = []
    set_times = []

    for size in sizes:
        unordered = unordered_sequence(max_len=size)
        ordered = ordered_sequence(max_len=size)
        as_set = set(unordered)
        target = unordered[-1] if unordered else 0

        linear_times.append(measure_time(linear_search, unordered, target))
        binary_times.append(measure_time(binary_search, ordered, target))
        set_times.append(measure_time(lambda s, t: t in s, as_set, target))

    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(sizes, linear_times, marker="o", label="linear_search")
    ax.plot(sizes, binary_times, marker="s", label="binary_search")
    ax.plot(sizes, set_times, marker="^", label="set membership")
    ax.set_xlabel("Input size n")
    ax.set_ylabel("Time [s]")
    ax.set_title("Comparison of search algorithm running times")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.grid(True, which="both", ls="--", alpha=0.4)
    ax.legend()
    fig.tight_layout()
    fig.savefig("search_times.png", dpi=120)
    plt.close(fig)

    return {
        "sizes": list(sizes),
        "linear": linear_times,
        "binary": binary_times,
        "set": set_times,
    }


def main():
    unordered = read_data("sequential.json", "unordered_numbers")
    ordered = read_data("sequential.json", "ordered_numbers")
    dna = read_data("sequential.json", "dna_sequence")

    print("Unordered numbers:", unordered)
    print("Ordered numbers:", ordered)
    print("DNA sequence:", dna)

    target = 5
    linear_result = linear_search(unordered, target)
    print(f"\nlinear_search(unordered_numbers, {target}) = {linear_result}")

    for target in (14, 63, 100):
        index = binary_search(ordered, target)
        print(f"binary_search(ordered_numbers, {target}) = {index}")

    pattern = "ATA"
    positions = pattern_search(dna, pattern)
    print(f"\npattern_search(dna_sequence, {pattern!r}) = {sorted(positions)}")
    print(f"number of occurrences: {len(positions)}")

    print("\nInvalid field returns:", read_data("sequential.json", "nonsense"))

    print("\nRunning time measurements (this may take a moment)...")
    results = compare_search_times()
    print("sizes :", results["sizes"])
    print("linear:", [f"{t:.6f}" for t in results["linear"]])
    print("binary:", [f"{t:.6f}" for t in results["binary"]])
    print("set   :", [f"{t:.6f}" for t in results["set"]])
    print("Plot saved to search_times.png")


if __name__ == "__main__":
    main()
