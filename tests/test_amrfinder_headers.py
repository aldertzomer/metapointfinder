import unittest
from itertools import permutations

from metapointfinder import pick_fields_from_hit


class TestAMRFinderHeaders(unittest.TestCase):
    def test_old_and_new_amrfinder_headers_parse_identically(self):
        suffix = "\t99\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200\tREFERENCESEQ\tTARGETSEQ\n"
        # Synthetic records: no particular gene or database ordering is required.
        records = (
            ("WP_000000001.1", "geneA", "example_protein_A"),
            ("NP_000000002.2", "geneB", "example_protein_B"),
            ("ABC00003.1", "geneC", "example_protein_C"),
        )

        for ordered_records in permutations(records):
            for prefix in ("0|", ""):
                for position, (accession, gene, description) in enumerate(ordered_records):
                    with self.subTest(
                        order=ordered_records, prefix=prefix, position=position
                    ):
                        read_id = f"read{position + 1}"
                        header = f"{accession}|1|1|{gene}|{gene}|mutation|2|||{description}"
                        line = f"{read_id}\t{prefix}{header}" + suffix
                        expected = (
                            read_id,
                            accession,
                            description,
                            "REFERENCESEQ",
                            "TARGETSEQ",
                        )
                        self.assertEqual(pick_fields_from_hit(line), expected)


if __name__ == "__main__":
    unittest.main()
