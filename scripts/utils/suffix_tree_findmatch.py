class SuffixTreeNode:
    """Node class for the Suffix Tree."""
    def __init__(self):
        self.children = {}
        self.start = -1
        self.end = -1

class SuffixTree:
    """Simple suffix tree implementation."""
    def __init__(self, text):
        self.text = text + "$"  # Add unique end marker
        self.root = SuffixTreeNode()
        self.build_tree()

    def build_tree(self):
        """Builds the suffix tree."""
        for i in range(len(self.text)):
            current_node = self.root
            for char in self.text[i:]:
                if char not in current_node.children:
                    new_node = SuffixTreeNode()
                    current_node.children[char] = new_node
                current_node = current_node.children[char]
                if current_node.start == -1:
                    current_node.start = i

    def find_matches(self, pattern, min_len):
        """Finds matches for a pattern of minimum length."""
        matches = []
        current_node = self.root
        length = 0
        for i, char in enumerate(pattern):
            if char in current_node.children:
                current_node = current_node.children[char]
                length += 1
                if length >= min_len:
                    matches.append((i - length + 1, i))
            else:
                break
        return matches

def extend_alignment(ref, query, ref_start, que_start, max_mismatches):
    """Extends an alignment."""
    mismatches = 0
    ref_idx, que_idx = ref_start, que_start
    while ref_idx < len(ref) and que_idx < len(query):
        if ref[ref_idx] != query[que_idx]:
            mismatches += 1
            if mismatches > max_mismatches:
                break
        ref_idx += 1
        que_idx += 1
    return ref_start, ref_idx, que_start, que_idx

def find_and_extend(ref, query, min_match_len, max_mismatches):
    """Finds and extends matches using seed-and-extend."""
    suffix_tree = SuffixTree(ref)
    matches = []
    for i in range(len(query)):
        mems = suffix_tree.find_matches(query[i:], min_match_len)
        for start, end in mems:
            alignment = extend_alignment(ref, query, start, i, max_mismatches)
            matches.append(alignment)
    return matches

# Example Usage
ref = "ATCGTACGATCGTTACGGTACG"
query = "GATCGTACGTTACGGATCGTAC"
min_match_len = 20
max_mismatches = 50

matches = find_and_extend(ref, query, min_match_len, max_mismatches)
print("Matches:")
for match in matches:
    print(f"REF: {match[0]}-{match[1]}, QUE: {match[2]}-{match[3]}")
