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


运行结果
实验1
[(0, 3597, 0, 3597), (3609, 3754, 3609, 3754), (3766, 5208, 3766, 5208), (5221, 6810, 5221, 6810), (6824, 6908, 6824, 6908), (6920, 6980, 6920, 6980), (7000, 7083, 7001, 7084), (7094, 7197, 7095, 7198), (7290, 7326, 7288, 7324), (7414, 7475, 7412, 7472), (7569, 7603, 7568, 7602), (7733, 7781, 7726, 7774), (7807, 7851, 7800, 7844), (7896, 7945, 7889, 7938), (7996, 8078, 7989, 8071), (12333, 12363, 14996, 15027), (14403, 14435, 21224, 21256), (16411, 16474, 2150, 2213), (16534, 16576, 2273, 2315), (16589, 16808, 2328, 2547), (16861, 17113, 2600, 2852), (17148, 17205, 2887, 2944), (17227, 17267, 2962, 3002), (17345, 17402, 3078, 3135), (17414, 17623, 3147, 3356), (17626, 17775, 25660, 25809), (17840, 17982, 25874, 26016), (17998, 18076, 26031, 26109), (18110, 18142, 26143, 26175), (18183, 18228, 26216, 26261), (18282, 18335, 26315, 26368), (18372, 18424, 26405, 26457), (21683, 21728, 21677, 21722), (21759, 21799, 21753, 21793), (21888, 21935, 21882, 21929), (21976, 22031, 21970, 22025), (22054, 22105, 22048, 22099), (22336, 22414, 22334, 22413), (22634, 22732, 22636, 22734), (22744, 22823, 22746, 22825), (22850, 23001, 22851, 23002), (23017, 25739, 23018, 25740), (25750, 26941, 25751, 26942), (26954, 29510, 26955, 29511), (29523, 29827, 29524, 29828)]
得分16350
实验2
[(0, 301, 0, 301), (507, 721, 607, 821), (768, 801, 718, 751), (805, 849, 605, 649), (1000, 1200, 700, 900), (1308, 1394, 908, 994), (1402, 1445, 402, 445), (1458, 1501, 458, 501), (1598, 1651, 1298, 1351), (1694, 1736, 1194, 1236), (1758, 1819, 1258, 1319), (1892, 1998, 1392, 1498), (2299, 2500, 1499, 1700)]
得分1370
