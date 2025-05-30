核心思路简述
构建k-mer哈希表：将ref中的所有k-mer映射为哈希表，用于快速查询query中匹配的k-mer。

搜索匹配点：对query及其反向互补序列滑动窗口获取k-mer，查询哈希表得到匹配点。

分离正向和反向匹配点，分别处理。

构建匹配链（锚点链），寻找最长递增路径表示匹配的连续片段。

严格处理链之间重叠，筛选最佳匹配组合。
# 伪代码
function mymatch(data, k):
    if data is empty:
        return []

    forward_matches = []
    reverse_matches = []

    for match in data:
        (query_pos, ref_pos, strand, kmer_size) = match
        if strand == 1:
            forward_matches.append((query_pos, ref_pos, kmer_size))
        else:
            reverse_matches.append((query_pos, ref_pos, kmer_size))

    forward_chains = find_chains_fixed(forward_matches, k)
    reverse_chains = find_chains_fixed(reverse_matches, k)

    all_chains = forward_chains + reverse_chains

    final_chains = resolve_overlaps_strict(all_chains, k)

    return final_chains

function find_chains_fixed(matches, min_length):
    if matches empty:
        return []

    if len(matches) > 5000:
        matches = sample_matches(matches, 5000)

    sort matches by query_pos and ref_pos

    initialize DP array

    for i in range(len(matches)):
        for j in range(i):
            if matches[j] can connect to matches[i]:
                update DP[i] with max chain length

    chains = traceback longest chains with length >= min_length

    return chains

function resolve_overlaps_strict(chains, k):
    sort chains by start positions
    select chains greedily avoiding overlap more than k
    return selected chains


时空复杂度分析
设参考序列长度为 m，查询序列长度为 n。

1. seq2hashtable_multi_test构建哈希表
生成ref所有k-mer：约 O(m)

查询query所有k-mer：约 O(n)

每个k-mer查询哈希表时间平均 O(1)

总时间：O(m + n)

空间：存储哈希表大小 O(m)

2. 匹配点分类与分组
遍历所有匹配点，分类正向反向：O(L)，其中L为匹配点数量（≤n，实际远小于n*m）

3. find_chains_fixed最长链构建
最坏时间复杂度是DP两层循环：O(L²)，L为匹配点数

但算法中有采样限制（最多处理5000个匹配点），使得实际复杂度常数上限：O(5000²) = O(1) (常数时间)

空间复杂度：存储DP数组和matches，O(L)

4. resolve_overlaps_strict链重叠处理
主要是排序和线性扫描：O(C log C)，C为链的数量，通常远小于L

空间复杂度：O(C)
