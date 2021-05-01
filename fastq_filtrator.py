import sys

# it is a block with default constant and flags
min_length = 0
keep_filtered = 0
gc_bounds_low = None
gc_bounds_high = 101
output_base_name = 'file_name'
four_lines = []
i = 0

# Usage
Inf = "This tool can filter your fastq file.\n\
Standard call of tools looks like this:\n\
python filter_fastq.py --some_argument some_values my_reads.fastq \n\n\
Optional arguments:\n\
    --min_length - minimal length of reads, which will be not filtered\n\
    --gc_bounds - range of necessary percent of GC content.\n\
        If one value is specified, this is the lower threshold, if 2 is specified, \n\
        then the lower and upper thresholds, respectively.\n\
    --output_base_name - common prefix name for all output files\n\
    --keep_filtered - flag for saving filtered reads."

# block for argument parsing
list_arg = [i for i in sys.argv]
output_base_name = list_arg[-1]
output_base_name = output_base_name.split('/')[-1].replace('.fastq', '')

help_list = ['-h', '--help']
possible_arguments = ['--min_length', '--gc_bounds', '--output_base_name', '--keep_filtered']

if not set(help_list).isdisjoint(list_arg) or len(list_arg) <= 2:
    print(Inf)
    sys.exit(0)
for i in list_arg:
    if "--" in i and i not in possible_arguments:
        print("Wrong arguments!")
        print(Inf)

for i in range(1, len(list_arg)):
    if i == len(list_arg) - 1:
        if list_arg[i].endswith('.fastq'):
            file = list_arg[i]
        else:
            print('The last argument must be file in fastq format!')
            sys.exit(1)
    elif list_arg[i] == '--min_length':
        try:
            int(list_arg[i + 1])
        except ValueError:
            print('Min_length value is not defined!')
            sys.exit(1)
        else:
            min_length = int(list_arg[i + 1])
    elif list_arg[i] == '--keep_filtered':
        keep_filtered = 1
    elif list_arg[i] == '--gc_bounds':
        try:
            int(list_arg[i + 1])
        except ValueError:
            print('GC_bounds value is not defined!')
            sys.exit(1)
        else:
            gc_bounds_low = int(list_arg[i + 1])
        try:
            int(list_arg[i + 2])
        except ValueError:
            pass
        else:
            if int(list_arg[i + 2]) > gc_bounds_low:
                gc_bounds_high = int(list_arg[i + 2])
            else:
                print('Wrong range of GC content!')
                sys.exit(1)
    elif list_arg[i] == '--output_base_name':
        if list_arg[i + 1].endswith('.fastq') or list_arg[i + 1].startswith('--'):
            print('Output base name is not specified!')
            sys.exit(1)
        else:
            output_base_name = list_arg[i + 1]
    elif list_arg[i] == '-h' or list_arg[i] == '--help':
        print(Inf)

# This part filter file according users setting
out_failed = str(output_base_name + '__failed' + '.fastq')
out_passed = str(output_base_name + '__passed' + '.fastq')

with open(out_passed, 'w') as f:
    pass
if keep_filtered == 1:
    with open(out_failed, 'w') as f:
        pass


def four_lines_read():
    with open(file, 'r') as f:
        res = []
        for line in f:
            res.append(line)
            if len(res) == 4:
                yield res
                res = []


def filter_(four_lines_):
    if keep_filtered == 1:
        with open(out_failed, 'a+') as f:
            for index_ in four_lines_:
                f.write(index_)


four_lines_iterator = four_lines_read()
while True:
    try:
        four_lines = next(four_lines_iterator)
    except StopIteration:
        break
    if len(four_lines[1][:-1]) < min_length:
        if keep_filtered == 1:
            filter_(four_lines)
        continue
    if gc_bounds_low is not None:
        nucl_content = {'G': 0, 'C': 0, 'A': 0, 'T': 0}
        for i in four_lines[1][:-1]:
            nucl_content[i] = nucl_content.get(i, 0) + 1
        GC_content = (nucl_content['C'] + nucl_content['G']) / len(four_lines[1][:-1])
        if GC_content < gc_bounds_low / 100 or GC_content > gc_bounds_high / 100:
            filter_(four_lines)
            continue
    with open(out_passed, 'a+') as f:
        for index in four_lines:
            f.write(index)

print('Process completed successfully')