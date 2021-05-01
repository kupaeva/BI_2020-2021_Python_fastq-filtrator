import sys

list_arg = [i for i in sys.argv]

inf = "This tool can filter your fastq file.\n\
Standard call of tools looks like this:\n\
python filter_fastq.py --some_argument some_values my_reads.fastq \n\n\
Optional arguments:\n\
    --min_length - minimal length of reads, which will be not filtered\n\
    --gc_bounds - range of necessary percent of GC content.\n\
        If one value is specified, this is the lower threshold, if 2 is specified, \n\
        then the lower and upper thresholds, respectively.\n\
    --output_base_name - common prefix name for all output files\n\
    --keep_filtered - flag for saving filtered reads."

if ('-h' in list_arg) or ('--help' in list_arg):
    print(inf)
    sys.exit(0)


def check_correctness(arg_list):
    possible_arguments = ['--min_length',
                          '--gc_bounds',
                          '--output_base_name',
                          '--keep_filtered']
    if len(arg_list) <= 4:
        raise ValueError('Not enough arguments! Run with -h or --help argument.')
    for i in arg_list:
        if "--" in i and i not in possible_arguments:
            raise ValueError('Wrong arguments! Run with -h or --help argument.')
    return True


def parse_file_name(arg_list):
    if arg_list[-1].endswith('.fastq'):
        return arg_list[-1]
    else:
        raise ValueError('The last argument must be file in fastq format!')


# must be int
def parse_min_length(arg_list):
    if '--min_length' in arg_list:
        i = arg_list.index('--min_length')
        if arg_list[i + 1].isnumeric():
            return int(arg_list[i + 1])
        else:
            raise ValueError('Min_length value is set, but not defined!')
    else:
        return 0


def parse_keep_filtered(arg_list):
    if '--keep_filtered' in arg_list:
        return 1
    else:
        return 0


# must be int
def parse_gc_bounds(arg_list):
    if '--gc_bounds' in arg_list:
        i = arg_list.index('--gc_bounds')
        if arg_list[i + 1].isnumeric():
            gc_low = int(arg_list[i + 1])
        else:
            raise ValueError('GC_bounds value is set, but not defined!')
        try:
            int(arg_list[i + 2])
        except ValueError:
            gc_high = 101
            pass
        else:
            if int(arg_list[i + 2]) > gc_low:
                gc_high = int(arg_list[i + 2])
            else:
                raise ValueError('Wrong range of GC content!')
        return gc_low, gc_high
    else:
        return None, 101


def parse_output_base_name(arg_list):
    if '--output_base_name' in arg_list:
        i = arg_list.index('--output_base_name')
        if arg_list[i + 1].endswith('.fastq') or arg_list[i + 1].startswith('--'):
            raise ValueError('Output base name is set, but not specified!')
        else:
            return arg_list[i + 1]
    else:
        name = arg_list[-1]
        return name.split('/')[-1].replace('.fastq', '')


def genere_output_names(output_base_name):
    out_failed = str(output_base_name + '__failed' + '.fastq')
    out_passed = str(output_base_name + '__passed' + '.fastq')
    return (out_failed, out_passed)


def genere_files(out_passed, out_failed, keep_filtered):
    with open(out_passed, 'w') as f:
        pass
    if keep_filtered == 1:
        with open(out_failed, 'w') as f:
            pass


# read lines
def line_reader(four_lines_iterator):
    try:
        four_lines = next(four_lines_iterator)
        return four_lines
    except StopIteration:
        return None


def four_lines_read(file):
    with open(file, 'r') as f:
        res = []
        for line in f:
            res.append(line)
            if len(res) == 4:
                yield res
                res = []


# filtration

def filter_(four_lines_, keep_filtered, out_failed):
    if keep_filtered == 1:
        with open(out_failed, 'a+') as f:
            for index_ in four_lines_:
                f.write(index_)


def count_gc(sequence):
    g = sequence.count('G')
    c = sequence.count('C')
    if len(sequence) > 0:
        return ((g + c) / len(sequence) * 100)
    else:
        return None


def write_file(out_passed, four_lines):
    with open(out_passed, 'a+') as f:
        for index in four_lines:
            f.write(index)


def filter_gc_content(four_lines, gc_bounds_low, gc_bounds_high):
    gc_content = count_gc(four_lines[1][:-1])
    if gc_content < gc_bounds_low or gc_content > gc_bounds_high:
        return 0
    return 1


if __name__ == '__main__':
    check = check_correctness(list_arg)
    file = parse_file_name(list_arg)
    min_length = parse_min_length(list_arg)
    keep_filtered = parse_keep_filtered(list_arg)
    gc_bounds_low, gc_bounds_high = parse_gc_bounds(list_arg)
    output_base_name = parse_output_base_name(list_arg)
    out_failed, out_passed = genere_output_names(output_base_name)
    genere_files(out_passed, out_failed, keep_filtered)
    four_lines_iterator = four_lines_read(file)
    while True:
        four_lines = line_reader(four_lines_iterator)
        print(four_lines)
        if four_lines is None:
            break
        if len(four_lines[1][:-1]) < min_length:
            if keep_filtered == 1:
                filter_(four_lines)
            continue
        if gc_bounds_low is not None:
            if filter_gc_content(four_lines, gc_bounds_low, gc_bounds_high) == 1:
                write_file(four_lines)
            continue
        else:
            filter_(four_lines, keep_filtered, out_failed)
            write_file(out_passed, four_lines)
            continue
    print('Process completed successfully')