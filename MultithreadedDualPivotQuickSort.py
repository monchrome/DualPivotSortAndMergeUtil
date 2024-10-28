import queue
import random
import time
import concurrent.futures
import os


class DualPivot:
    def __init__(self, left, right):
        self.left = left
        self.right = right


def dual_pivot_quick_sort(nums: list):
    '''
    Top most function called by each worker thread with an input list to be sorted in place
    :param nums: list of ints to be sorted
    :return:
    '''
    if nums is None or len(nums) == 0:
        raise Exception("Input list is empty.")

    quick_sort_using_dual_pivots(nums, 0, len(nums) - 1)


def quick_sort_using_dual_pivots(nums: list, low_index: int, high_index: int):
    '''
     Partitions and the recursively sorts the given list using dual pivot quick sort logic
    :param nums: list of ints to be sorted
    :param low_index: left boundary, inclusive
    :param high_index: right boundary, inclusive
    :return:
    '''
    if (low_index >= high_index):
        return
    # Step 1 : Partition the input list
    pivot_record = partition(nums, low_index, high_index)
    # Step 2 : Sort such that all numbers < than left pivot get shifted to the left side of left pivot element
    quick_sort_using_dual_pivots(nums, low_index, pivot_record.left - 1)
    # Step 3 : Sort suc that numbers that fall in between left pivot and right pivot, get shifted in between those two pivots
    quick_sort_using_dual_pivots(nums, pivot_record.left + 1, pivot_record.right - 1)
    # Step 4 : Sort such that all numbers > than right pivot get shifted to the right side of right pivot element
    quick_sort_using_dual_pivots(nums, pivot_record.right + 1, high_index)


def partition(nums: list, low_index: int, high_index: int) -> DualPivot:
    '''
    Partitions the given list and does in place sorting based on left and right pivot
    :param nums: list with int items to be sorted
    :param low_index: index position to be treated as left pivot
    :param high_index: index position to be treated as right pivot
    :return:
    '''
    # Assume element at low_index to be left pivot and element at high _index to be the right pivot
    # Left pivot should be lower than right pivot, if not then swap
    if nums[low_index] > nums[high_index]:
        swap(nums, low_index, high_index)

    left_ptr, curr_ptr, right_ptr, = low_index + 1, low_index + 1, high_index - 1
    while curr_ptr <= right_ptr:
        # If current element is lower than left pivot then swap
        if nums[curr_ptr] < nums[low_index]:
            swap(nums, curr_ptr, left_ptr)
            curr_ptr += 1
            left_ptr += 1
        # If current element if higher than right pivot then swap
        elif nums[curr_ptr] > nums[high_index]:
            swap(nums, curr_ptr, right_ptr)
            right_ptr -= 1
        # If elements are equal then switch to next element
        else:
            curr_ptr += 1

    # Swap left pivot and leftPtr and likewise right pivot and rightPtr
    left_ptr -= 1
    swap(nums, low_index, left_ptr)
    right_ptr += 1
    swap(nums, high_index, right_ptr)
    return DualPivot(left_ptr, right_ptr)


def swap(nums: list, frm: int, to: int):
    tmp = nums[frm]
    nums[frm] = nums[to]
    nums[to] = tmp


# Total Time Taken for an array of size 10000000 is 17451.538801193237 ms
# Total Time Taken for an array of size 100000000 is 312287.6877784729 ms
# does not finish with  100000000 is 312287.6877784729 ms
# def test_dual_pivot_sort():
#     maxlimit = 1000000000
#     nums = []
#     for i in range(0, maxlimit):
#         nums.append(random.randint(1, 1000000))
#     start_time = time.time()
#     dual_pivot_quick_sort(nums)
#     elapsed_time = (time.time() - start_time) * 1000

# len_nums = len(nums)
# for i in range(1, len_nums):
#     #print(nums[i])
#     if nums[i] == nums[i-1]:
#         continue
#     if nums[i] < nums[i - 1]:
#         raise Exception(f"Not sorted completely. Found anomaly : {nums[i] ,} {nums[i-1]}")
# print(f"Total Time Taken for an array of size {maxlimit} is {elapsed_time} ms")

# test_dual_pivot_sort()

def worker(max_limit_per_list: int, output_file_name: str):
    '''
    Worker function does the generation of input data to be sorted. Each call to worker thread
    generates a list of random ints containing max_limit_per_list items.
    output_file_name is name of the file where sorted list of ints is stored after worker thread completes
    dual pivot quick sorting.
    In production, I had time stamped data via incoming files which had to be sorted.
    Here just doing simulation by generating list of ints.
    :param max_limit_per_list: Max number of int type items to be generated in a list.
    :param output_file_name: output file name with sorted items. File is stored in same directory as
    the executable
    :return:
    '''
    input_list = []
    for i in range(1, max_limit_per_list):
        input_list.append(random.randint(1, 10000000))
    dual_pivot_quick_sort(input_list)

    lines = [str(item) + '\n' for item in input_list]
    with open(f"./{output_file_name}", 'a') as res_file:
        res_file.writelines(lines)

# configurable properties
max_threads = 3
max_limit_per_list = 10000000
batch_size = 50000

'''
Total time taken for sorting 3 lists, containing 100000000 items each, is 878124 ms.
Total time taken for merging is 119543 ms.
'''

pool = concurrent.futures.ThreadPoolExecutor(max_workers=max_threads)
start_time = time.time()
output_files = []
for i in range(1, 4):
    output_file_name = "sorted_list_" + str(i) + ".txt"
    # delete files from last run, if they exist
    if os.path.exists(output_file_name):
        os.remove(output_file_name)
    # launch worker thread
    output_files.append(output_file_name)
    pool.submit(worker(max_limit_per_list, output_file_name))

pool.shutdown(wait=True)
elapsed_time = round((time.time() - start_time) * 1000)
print(
    f"Total time taken for sorting {max_threads} lists, containing {max_limit_per_list} items each, is {elapsed_time} ms.")

# Merge all the sorted files into one single file
# Number of sorted lists in result_list should not exceed number of worker threads launched earlier
result_lists = []
total_items_from_input_files = 0
for file_name in output_files:
    with open(file_name) as input_file:
        lines = [int(line.rstrip()) for line in input_file]
        total_items_from_input_files += len(lines)
        result_lists.append(lines)
assert len(result_lists) == max_threads, f"Expected {max_threads} files. Found just {len(result_lists)}"

final_merged_file = "./final_merged_list.txt"
if os.path.exists(final_merged_file):
    os.remove(final_merged_file)
total_items_merged = 0
print('Starting merge of sorted items.....')
i, j, k = 0, 0, 0
write_batch = []
list1_len = len(result_lists[0])
list2_len = len(result_lists[1])
list3_len = len(result_lists[2])
min_val = float('inf')
start_time = time.time()
while i < list1_len or j < list2_len or k < list3_len:
    if i == float('inf') and j == float('inf') and k == float('inf'):
        break

    if i < len(result_lists[0]):
        i_val = result_lists[0][i]
    else:
        i_val = float('inf')
    if j < len(result_lists[1]):
        j_val = result_lists[1][j]
    else:
        j_val = float('inf')
    if k < len(result_lists[2]):
        k_val = result_lists[2][k]
    else:
        k_val = float('inf')

    min_val = min(i_val, min(j_val, k_val))
    write_batch.append(min_val)
    # Since each input file can contain upto 1 Billion entries, creating a list to hold 3 Billion items at a time,
    # can cause memory issues if the container is memory strapped. Hence doing batched writes of merged items.
    if len(write_batch) >= batch_size:
        lines = [str(item) + '\n' for item in write_batch]
        total_items_merged += len(write_batch)
        with open(final_merged_file, 'a') as file:
            file.writelines(lines)
            file.flush()
        write_batch.clear()
    if i < list1_len and result_lists[0][i] == min_val:
        i += 1
    if j < list2_len and result_lists[1][j] == min_val:
        j += 1
    if k < list3_len and result_lists[2][k] == min_val:
        k += 1

# last batch may contain items fewer than batch_size, so just flushing them to disk finally
if len(write_batch) > 0:
    lines = [str(item) + '\n' for item in write_batch]
    total_items_merged += len(write_batch)
    with open(final_merged_file, 'a') as file:
        file.writelines(lines)
        file.flush()
elapsed_time = round((time.time() - start_time) * 1000)
print(f"Total time taken for merging is {elapsed_time} ms.")
