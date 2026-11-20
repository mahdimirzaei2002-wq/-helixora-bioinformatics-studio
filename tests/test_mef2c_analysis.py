import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation

from modules.mef2c_analysis import catalogue, primer_feasibility

class Mef2cAnalysisTests(unittest.TestCase):
    def test_catalogue_distinguishes_coding_and_non_coding_records(self):
        coding = SeqRecord(Seq("ATGAAATAA"), id="NM_TEST.1", description="MEF2C transcript variant 1")
        coding.features = [SeqFeature(FeatureLocation(0, 9), type="CDS", qualifiers={"protein_id": ["NP_TEST.1"], "translation": ["MK*"]})]
        noncoding = SeqRecord(Seq("A" * 50), id="NR_TEST.1", description="MEF2C non-coding transcript")
        df = catalogue([coding, noncoding])
        self.assertEqual(df.loc[0, "Protein accession"], "NP_TEST.1")
        self.assertEqual(df.loc[1, "Coding status"], "No CDS feature")

    def test_unique_window_screen(self):
        a = SeqRecord(Seq("ACGT" * 30), id="NM_A.1")
        b = SeqRecord(Seq("TTTT" * 30), id="NM_B.1")
        df = primer_feasibility([a, b])
        self.assertTrue((df["Unique 22-nt windows"] > 0).all())

if __name__ == "__main__":
    unittest.main()
