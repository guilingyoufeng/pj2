import numpy as np
from collections import defaultdict

def mymatch(data, k=10):

    if len(data) == 0:
        return []
    
    # 1. 数据预处理 - 按链分组
    forward_matches = []
    reverse_matches = []
    
    for match in data:
        query_pos, ref_pos, strand, kmer_size = match
        query_pos, ref_pos, kmer_size = int(query_pos), int(ref_pos), int(kmer_size)
        
        if strand == 1:
            forward_matches.append((query_pos, ref_pos, kmer_size))
        else:
            reverse_matches.append((query_pos, ref_pos, kmer_size))
    
    # 2. 分别处理正向和反向匹配
    forward_chains = find_chains_fixed(forward_matches, k)
    reverse_chains = find_chains_fixed(reverse_matches, k)
    
    # 3. 合并所有链
    all_chains = forward_chains + reverse_chains
    
    # 4. 严格解决重叠并选择最佳组合
    final_chains = resolve_overlaps_strict(all_chains, k)
    
    return final_chains

def find_chains_fixed(matches, min_length):
    """
    修复的链查找算法，控制复杂度
    """
    if not matches:
        return []
    
    # 限制处理的匹配数量，避免过度计算
    if len(matches) > 5000:
        # 按密度采样，保持代表性
        matches = sample_matches_by_density(matches, 5000)
    
    # 按query位置排序
    matches.sort(key=lambda x: x[0])
    
    # 使用简化的动态规划
    n = len(matches)
    dp = [0] * n
    parent = [-1] * n
    
    for i in range(n):
        query_pos_i, ref_pos_i, kmer_size_i = matches[i]
        dp[i] = kmer_size_i
        
        # 只检查前面的一定数量的匹配，控制复杂度
        start_j = max(0, i - 100)  # 限制回溯范围
        
        for j in range(start_j, i):
            query_pos_j, ref_pos_j, kmer_size_j = matches[j]
            
            # 简化的链接条件
            if can_chain_simple(matches[j], matches[i]):
                score = dp[j] + kmer_size_i
                
                if score > dp[i]:
                    dp[i] = score
                    parent[i] = j
    
    # 简化的回溯
    chains = extract_chains_simple(matches, dp, parent, min_length)
    
    return chains

def sample_matches_by_density(matches, max_count):
    """
    按密度采样匹配，保持代表性
    """
    if len(matches) <= max_count:
        return matches
    
    # 按query位置排序
    matches.sort(key=lambda x: x[0])
    
    # 计算采样步长
    step = len(matches) / max_count
    sampled = []
    
    i = 0
    while len(sampled) < max_count and i < len(matches):
        sampled.append(matches[int(i)])
        i += step
    
    return sampled

def can_chain_simple(match1, match2):
    """
    简化的链接判断
    """
    q1, r1, k1 = match1
    q2, r2, k2 = match2
    
    # 基本顺序检查
    if q2 <= q1 or r2 <= r1:
        return False
    
    # 距离不能太远
    if q2 - q1 > 1000 or r2 - r1 > 1000:
        return False
    
    # 计算间隔差异
    query_gap = q2 - (q1 + k1)
    ref_gap = r2 - (r1 + k1)
    gap_diff = abs(query_gap - ref_gap)
    
    # 允许的最大差异
    max_diff = min(50, max(q2 - q1, r2 - r1) * 0.15)
    
    return gap_diff <= max_diff

def extract_chains_simple(matches, dp, parent, min_length):
    """
    简化的链提取算法
    """
    n = len(matches)
    used = [False] * n
    chains = []
    
    # 按得分排序，优先处理高得分的链
    candidates = [(dp[i], i) for i in range(n) if dp[i] >= min_length]
    candidates.sort(reverse=True)
    
    for score, end_idx in candidates:
        if used[end_idx]:
            continue
        
        # 回溯构建链
        chain_matches = []
        current = end_idx
        
        while current != -1:
            if not used[current]:
                chain_matches.append(matches[current])
                used[current] = True
            current = parent[current]
        
        if chain_matches:
            chain_matches.reverse()
            # 转换为区间
            intervals = convert_to_intervals(chain_matches, min_length)
            chains.extend(intervals)
    
    return chains

def convert_to_intervals(chain_matches, min_length):
    """
    将匹配链转换为区间
    """
    if not chain_matches:
        return []
    
    intervals = []
    chain_matches.sort(key=lambda x: x[0])
    
    i = 0
    while i < len(chain_matches):
        # 开始新区间
        query_start = chain_matches[i][0]
        query_end = chain_matches[i][0] + chain_matches[i][2]
        ref_start = chain_matches[i][1]
        ref_end = chain_matches[i][1] + chain_matches[i][2]
        
        # 尝试扩展区间
        j = i + 1
        while j < len(chain_matches):
            next_q, next_r, next_k = chain_matches[j]
            
            # 检查是否可以合并
            q_gap = next_q - query_end
            r_gap = next_r - ref_end
            
            if q_gap <= 10 and r_gap <= 10 and abs(q_gap - r_gap) <= 5:
                # 扩展区间
                query_end = max(query_end, next_q + next_k)
                ref_end = max(ref_end, next_r + next_k)
                j += 1
            else:
                break
        
        # 检查区间长度
        if query_end - query_start >= max(min_length, 20):
            intervals.append((query_start, query_end, ref_start, ref_end))
        
        i = j
    
    return intervals

def resolve_overlaps_strict(chains, min_length):
    """
    严格解决重叠问题，确保query上没有重叠
    """
    if not chains:
        return []
    
    # 过滤长度
    valid_chains = []
    for chain in chains:
        if chain[1] - chain[0] >= max(min_length, 30):  # 提高最小长度要求
            valid_chains.append(chain)
    
    if not valid_chains:
        return []
    
    # 按query起始位置排序
    valid_chains.sort(key=lambda x: x[0])
    
    # 使用贪心算法选择不重叠的区间
    selected = []
    last_query_end = -1
    
    for chain in valid_chains:
        query_start, query_end, ref_start, ref_end = chain
        
        # 严格检查：新区间必须在前一个区间结束之后
        if query_start >= last_query_end:
            selected.append(chain)
            last_query_end = query_end
    
    # 如果贪心结果太少，尝试动态规划优化
    if len(selected) < len(valid_chains) * 0.3:
        selected = dynamic_select_non_overlapping(valid_chains)
    
    # 最终转换为标准int类型
    result = []
    for chain in selected:
        result.append((int(chain[0]), int(chain[1]), int(chain[2]), int(chain[3])))
    
    return result

def dynamic_select_non_overlapping(chains):
    """
    动态规划选择最优的非重叠区间组合
    """
    if not chains:
        return []
    
    n = len(chains)
    if n == 1:
        return chains
    
    # 按query_end排序
    chains_sorted = sorted(chains, key=lambda x: x[1])
    
    # 动态规划
    dp = [0] * n
    parent = [-1] * n
    
    # 计算每个区间的权重（长度）
    weights = [c[1] - c[0] for c in chains_sorted]
    dp[0] = weights[0]
    
    for i in range(1, n):
        # 不选择当前区间
        dp[i] = dp[i-1]
        parent[i] = parent[i-1]
        
        # 选择当前区间
        current_weight = weights[i]
        best_prev = -1
        
        # 找到最近的不重叠区间
        for j in range(i-1, -1, -1):
            if chains_sorted[j][1] <= chains_sorted[i][0]:
                best_prev = j
                break
        
        prev_value = dp[best_prev] if best_prev >= 0 else 0
        if prev_value + current_weight > dp[i]:
            dp[i] = prev_value + current_weight
            parent[i] = i
    
    # 回溯构建解
    result = []
    i = n - 1
    while i >= 0:
        if parent[i] == i:  # 选择了当前区间
            result.append(chains_sorted[i])
            # 跳到前一个选择的区间
            while i >= 0 and chains_sorted[i][1] > chains_sorted[parent[i]][0]:
                i -= 1
        else:
            i -= 1
    
    result.reverse()
    return result


#评分部分

import numpy as np
from numba import njit
import edlib


def get_rc(s):
    map_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    l = []
    for c in s:
        l.append(map_dict[c])
    l = l[::-1]
    return ''.join(l)
def rc(s):
    map_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    l = []
    for c in s:
        l.append(map_dict[c])
    l = l[::-1]
    return ''.join(l)

def seq2hashtable_multi_test(refseq, testseq, kmersize=15, shift = 1):
    rc_testseq = get_rc(testseq)
    testseq_len = len(testseq)
    local_lookuptable = dict()
    skiphash = hash('N'*kmersize)
    for iloc in range(0, len(refseq) - kmersize + 1, 1):
        hashedkmer = hash(refseq[iloc:iloc+kmersize])
        if(skiphash == hashedkmer):
            continue
        if(hashedkmer in local_lookuptable):

            local_lookuptable[hashedkmer].append(iloc)
        else:
            local_lookuptable[hashedkmer] = [iloc]
    iloc = -1
    readend = testseq_len-kmersize+1
    one_mapinfo = []
    preiloc = 0
    while(True):
   
        iloc += shift
        if(iloc >= readend):
            break

        #if(hash(testseq[iloc: iloc + kmersize]) == hash(rc_testseq[-(iloc + kmersize): -iloc])):
            #continue
 
        hashedkmer = hash(testseq[iloc: iloc + kmersize])
        if(hashedkmer in local_lookuptable):

            for refloc in local_lookuptable[hashedkmer]:

                one_mapinfo.append((iloc, refloc, 1, kmersize))



        hashedkmer = hash(rc_testseq[-(iloc + kmersize): -iloc])
        if(hashedkmer in local_lookuptable):
            for refloc in local_lookuptable[hashedkmer]:
                one_mapinfo.append((iloc, refloc, -1, kmersize))
        preiloc = iloc

    

    return np.array(one_mapinfo)

def get_points(tuples_str):
    data = []
    num = 0
    for c in tuples_str:
        if(ord('0') <= c <= ord('9')):
            num = num * 10 + c - ord('0')
        elif(ord(',') == c):
            data.append(num)
            num = 0
    if(num != 0):
        data.append(num)
    return data

def calculate_distance(ref, query, ref_st, ref_en, query_st, query_en):
    A = ref[ref_st: ref_en]
    a = query[query_st: query_en]
    _a = rc(query[query_st: query_en])
    return min(edlib.align(A, a)['editDistance'], edlib.align(A, _a)['editDistance'])

def get_first(x):
    return x[0]


def calculate_value(tuples_str, ref, query):  

    slicepoints = np.array(get_points(tuples_str.encode()))
    if(len(slicepoints) > 0 and len(slicepoints) % 4 == 0):
        editdistance = 0
        aligned = 0
        preend = 0
        points = np.array(slicepoints).reshape((-1, 4)).tolist()
        points.sort(key=get_first)
        for onetuple in points:
            query_st, query_en, ref_st, ref_en = onetuple
            if(preend > query_st):
                return 0
            if(query_en - query_st < 30):
                continue
            preend = query_en
            if((calculate_distance(ref, query, ref_st, ref_en, query_st, query_en)/len(query[query_st:query_en])) > 0.1):
                continue
            editdistance += calculate_distance(ref, query, ref_st, ref_en, query_st, query_en)
            aligned += len(query[query_st:query_en])
        return max(aligned - editdistance, 0)
    else:
        return 0
data = seq2hashtable_multi_test(ref, query, kmersize=9, shift = 1)
data.shape
tuples_str = str(mymatch(data, k=10))
print(tuples_str)
#Score
calculate_value(tuples_str, ref, query)

