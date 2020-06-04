def oddeven_merge(lo: int, hi: int, r: int):
    step = r * 2
    if step <= hi - lo:
        yield from oddeven_merge(lo, hi, step)
        yield from oddeven_merge(lo + r, hi, step)
        yield from [(i, i + r) for i in range(lo + r, hi - r, step)]
    else:
        yield (lo, lo + r)

def oddeven_merge_sort_range(lo: int, hi: int):
    """ sort the part of x with indices between lo and hi.

    Note: endpoints (lo and hi) are included.
    """
    if (hi - lo) >= 1:
        # if there is more than one element, split the input
        # down the middle and first sort the first and second
        # half, followed by merging them.
        mid = lo + ((hi - lo) // 2)
        yield from oddeven_merge_sort_range(lo, mid)
        yield from oddeven_merge_sort_range(mid + 1, hi)
        yield from oddeven_merge(lo, hi, 1)

def oddeven_merge_sort(length: int):
    """ "length" is the length of the list to be sorted.
    Returns a list of pairs of indices starting with 0 """
    yield from oddeven_merge_sort_range(0, length - 1)

def compare_and_swap(x, a, b) -> None:
    if x[a] > x[b]:
        x[a], x[b] = x[b], x[a]

n=32

pairs_to_compare = list(oddeven_merge_sort(n))

print("\n####################################################")
print("###################### n =",n,"######################")
print("####################################################\n")

max_depth = 0
for i in range(0,n):
    depth = 0
    for j in pairs_to_compare:
        if (j[0] == i) or (j[1] == i):
            depth = depth + 1
    if depth > max_depth:
        max_depth = depth

print("depth of the network:",max_depth,"\n")

print(len(pairs_to_compare),"pairs to compare\n")

lvl = 1;
copy = pairs_to_compare[:]
while pairs_to_compare != []:
    level = []
    check_depth = [0]*n
    for i in pairs_to_compare:
        if check_depth[i[0]] == 0 and check_depth[i[1]] == 0:
            check_depth[i[0]] = 1
            check_depth[i[1]] = 1
            level.append(i)
            copy.remove(i)
            
    pairs_to_compare = copy[:]
    print("Level",lvl,":",len(level),"pairs","\n",level,"\n")
    lvl = lvl+1
