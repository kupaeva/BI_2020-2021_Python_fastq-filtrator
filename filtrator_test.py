import unittest
import os.path
from fastq_filtrator import genere_output_names, genere_files, line_reader, four_lines_read, filter_, count_gc, write_file, \
    filter_gc_content

from filtrator import check_correctness, parse_file_name, parse_min_length, \
    parse_keep_filtered, parse_gc_bounds, parse_output_base_name


class TestParsing(unittest.TestCase):
    def setUp(self):
        self.arguments_1 = ['python', 'Fastq_filtrator.py',  # correct record
                            '--min_length', '10',
                            '--gc_bounds', '40', '60',
                            '--keep_filtered',
                            '--output_base_name', 'out',
                            'file.fastq']
        self.arguments_2 = ['python', 'Fastq_filtrator.py',  # less arguments
                            '--help']
        self.arguments_3 = ['python', 'Fastq_filtrator.py',  # mistake in arguments
                            '--min_len', '10',
                            '--keep_filtered',
                            'file.fast']
        self.arguments_4 = ['value', 'file.fastq']  # for check if no arg
        self.arguments_5 = ['--min_length', '10.5']  # incorrect/not defined value min_length
        self.arguments_6 = ['--gc_bounds', '40.5']  # incorrect/not defined value gc_bounds
        self.arguments_7 = ['--gc_bounds', '60', '40']  # incorrect order gc_bounds
        self.arguments_8 = ['--gc_bounds', '40', 'value']  # one value gc_bounds
        self.arguments_9 = ['--output_base_name',  # # output_base_name is not specified
                            'file.fastq']
        self.keep_filtered = 1
        self.not_keep_filtered = 0
        self.gc_bounds_low = None
        self.gc_bounds_low_ziro = 0
        self.gc_bounds_low_50 = 50
        self.gc_bounds_low_100 = 100
        self.gc_bounds_high = 101
        self.gc_bounds_high_ziro = 0
        self.gc_bounds_high_50 = 50
        self.gc_bounds_high_100 = 100
        self.min_length = 0
        self.min_length_20 = 20
        self.output_base_name = 'file_name'
        self.file = 'test.fastq'
        self.four_lines_1 = ['@1\n', 'GGTTGCAGGG\n', '+\n', '@?:=:;DB\n']
        self.four_lines_2 = ['@3\n',
                             'GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG\n', \
                             '+\n',
                             'CCCFFFFFHHHHHJJJJJJJJJJFFHIJJJJJJJJJJJJJJJJJJJJJJJIJHHHHHHFDEDF;AEEEEEEDDDDDBBACDDDCDDDDCCDDDDDDCCDC\n']
        self.four_lines_3 = ['@4\n',
                             'AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA\n', \
                             '+\n',
                             'CCCFFFFFHHHHHJHIIJIIIIJJJJJJGIJJJJJIJJIIGHJJJJJIIJJDHFFFFFEDACDDDCDDDDCCDDECACCDCCCDACDDDDCCDDDDDBD@A\n']

    def test_check_correctness(self):
        self.assertEqual(check_correctness(self.arguments_1), True)
        with self.assertRaises(ValueError):
            check_correctness(self.arguments_2)
        with self.assertRaises(ValueError):
            check_correctness(self.arguments_3)

    def test_parse_file_name(self):
        self.assertEqual(parse_file_name(self.arguments_1), 'file.fastq')
        with self.assertRaises(ValueError):
            parse_file_name(self.arguments_3)

    def test_parse_min_length(self):
        self.assertEqual((parse_min_length(self.arguments_1),
                          parse_min_length(self.arguments_4)), (10, 0))
        with self.assertRaises(ValueError):
            parse_min_length(self.arguments_5)

    def test_parse_keep_filtered(self):
        self.assertEqual((parse_keep_filtered(self.arguments_1),
                          parse_keep_filtered(self.arguments_4)), (1, 0))

    def test_parse_gc_bounds(self):
        self.assertEqual((parse_gc_bounds(self.arguments_1),
                          parse_gc_bounds(self.arguments_4),
                          parse_gc_bounds(self.arguments_8)),
                         ((40, 60), (None, 101), (40, 101)))
        with self.assertRaises(ValueError):
            parse_gc_bounds(self.arguments_6)
        with self.assertRaises(ValueError):
            parse_gc_bounds(self.arguments_7)

    def test_parse_output_base_name(self):
        self.assertEqual((parse_output_base_name(self.arguments_1),
                          parse_output_base_name(self.arguments_4)), ('out', 'file'))
        with self.assertRaises(ValueError):
            parse_output_base_name(self.arguments_9)

    def test_genere_output_names(self):
        self.out_failed, self.out_passed = genere_output_names(self.output_base_name)
        self.assertEqual(self.out_failed, 'file_name__failed.fastq')
        self.assertEqual(self.out_passed, 'file_name__passed.fastq')

    def test_genere_files(self):
        genere_files('file_name__passed.fastq', 'file_name__failed.fastq', self.keep_filtered)
        check_file = os.path.exists('file_name__failed.fastq')  # True
        self.assertEqual(check_file, True)

        size_passed = os.path.getsize('file_name__passed.fastq') > 0
        self.assertEqual(size_passed, False)

        size_failed = os.path.getsize('file_name__failed.fastq') > 0
        self.assertEqual(size_failed, False)

        os.remove('file_name__passed.fastq')
        os.remove('file_name__failed.fastq')

        genere_files('file_name__passed.fastq', 'file_name__failed.fastq', self.not_keep_filtered)
        check_file = os.path.exists('file_name__failed.fastq')  # True
        self.assertEqual(check_file, False)

        os.remove('file_name__passed.fastq')

    def test_write_file(self):
        genere_files('file_name__passed.fastq', 'file_name__failed.fastq', self.keep_filtered)
        filter_(self.four_lines_1, self.keep_filtered, 'file_name__failed.fastq')
        size_failed = os.path.getsize('file_name__failed.fastq') > 0
        self.assertEqual(size_failed, True)

        write_file('file_name__passed.fastq', self.four_lines_1)
        size_passed = os.path.getsize('file_name__failed.fastq') > 0
        self.assertEqual(size_passed, True)

        os.remove('file_name__passed.fastq')
        os.remove('file_name__failed.fastq')

    def test_four_lines_read(self):
        four_lines_iterator = four_lines_read(self.file)
        four_lines = line_reader(four_lines_iterator)
        self.assertEqual(four_lines, self.four_lines_1)

    def test_count_gc(self):
        gc_2 = count_gc(self.four_lines_2[1][:-1])
        self.assertEqual(gc_2, 100)
        gc_3 = count_gc(self.four_lines_3[1][:-1])
        self.assertEqual(gc_3, 0)
        gc_1 = count_gc(self.four_lines_1[1][:-1])
        self.assertEqual(gc_1, 70)

    def test_filter_gc_content(self):
        gc_2_yes = filter_gc_content(self.four_lines_2, self.gc_bounds_low_ziro, self.gc_bounds_high)
        self.assertEqual(gc_2_yes, 1)
        gc_2_not = filter_gc_content(self.four_lines_2, self.gc_bounds_low_ziro, self.gc_bounds_high_50)
        self.assertEqual(gc_2_not, 0)
        gc_3_not = filter_gc_content(self.four_lines_3, self.gc_bounds_low_50, self.gc_bounds_high)
        self.assertEqual(gc_3_not, 0)
        gc_3_yes = filter_gc_content(self.four_lines_3, self.gc_bounds_low_ziro, self.gc_bounds_high_ziro)
        self.assertEqual(gc_3_yes, 1)


if __name__ == "__main__":
    unittest.main()
