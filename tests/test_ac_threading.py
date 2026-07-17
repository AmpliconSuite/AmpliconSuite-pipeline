import unittest

from paalib.ac_threading import compute_ac_thread_allocation


class ACThreadAllocationTests(unittest.TestCase):
    def test_single_amplicon_gets_full_budget(self):
        self.assertEqual(compute_ac_thread_allocation(8, 1), (1, 8))

    def test_multiple_amplicons_favor_parallel_jobs(self):
        self.assertEqual(compute_ac_thread_allocation(8, 20), (4, 2))

    def test_jobs_do_not_exceed_amplicon_count(self):
        self.assertEqual(compute_ac_thread_allocation(16, 2), (2, 3))

    def test_one_thread_falls_back_to_one_by_one(self):
        self.assertEqual(compute_ac_thread_allocation(1, 10), (1, 1))

    def test_invalid_thread_budget_is_normalized(self):
        self.assertEqual(compute_ac_thread_allocation("invalid", 3), (1, 1))


if __name__ == "__main__":
    unittest.main()
