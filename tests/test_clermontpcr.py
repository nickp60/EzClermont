#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import re
import argparse
import sys
import os
import unittest

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import cpcr.run as clp

@unittest.skipIf((sys.version_info[0] != 3) or (sys.version_info[1] < 5),
                 "Subprocess.call among other things wont run if tried " +
                 " with less than python 3.5")
class clermontTestCase(unittest.TestCase):
    """
    """
    def test_interpret(self):
        """ Ensure that the decision tree for parsing the 4 primer sets works right
        See the clermont 2013 paper
        """

        ref = [
            'A', 'A/C', 'B1', 'A/C', 'D/E', 'D/E', 'E/cryptic', 'D/E', 'D/E',
            'F/G', 'B2', 'B2/G', 'B2', 'cryptic', 'E/cryptic', 'U/cryptic'
        ]
        test = []
        # if True:
        # A's
        test.append(
            clp.interpret_hits(arpA=True, chu=False, yjaA=False, TspE4=False))
        test.append(
            clp.interpret_hits(arpA=True, chu=False, yjaA=True, TspE4=False))
        # B1
        test.append(
            clp.interpret_hits(arpA=True, chu=False, yjaA=False, TspE4=True))
        # C
        test.append(
            clp.interpret_hits(arpA=True, chu=False, yjaA=True, TspE4=False))
        # E
        test.append(
            clp.interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=False))
        test.append(
            clp.interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=True))
        test.append(
            clp.interpret_hits(arpA=True, chu=True, yjaA=True, TspE4=False))
        # D
        test.append(
            clp.interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=False))
        test.append(
            clp.interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=True))
        # F/G
        test.append(
            clp.interpret_hits(arpA=False, chu=True, yjaA=False, TspE4=False))
        # B2
        test.append(
            clp.interpret_hits(arpA=False, chu=True, yjaA=True, TspE4=False))
        # B2/G
        test.append(
            clp.interpret_hits(arpA=False, chu=True, yjaA=False, TspE4=True))
        # B2
        test.append(
            clp.interpret_hits(arpA=False, chu=True, yjaA=True, TspE4=True))
        # cryptic
        test.append(
            clp.interpret_hits(arpA=False, chu=False, yjaA=True, TspE4=False))
        test.append(
            clp.interpret_hits(arpA=True, chu=True, yjaA=True, TspE4=False))
        # unknown
        test.append(
            clp.interpret_hits(arpA=False, chu=False, yjaA=False, TspE4=False))
        print(ref)
        print(test)
        self.assertEqual(ref, test)

    def test_get_matches(self):
        """ test creation of PCRhit object from two string matches
        """
        test_seq = [SeqRecord(
            Seq("GCACAGTCGATCAAAATTTTTGCAGTCGACTGGACTGACTGTCGGATCTCAGTCAT"))]
        test_fwd = "GCACAG"
        fwd = re.compile(test_fwd, re.IGNORECASE)
        test_rev = "TGACTG"
        self.assertEqual(fwd.search(str(test_seq[0].seq)).span(), (0, 6))
        forward_control_matches = clp.get_matches(
            allele="test",
            seq_list=test_seq,
            fwd_primer=test_fwd,
            rev_primer=test_rev,
            allow_partial=True,
            expected_size=60,
            strand='+')
        self.assertEqual(forward_control_matches[0].R_end, 55)

    def test_refine_hits(self):
        """ TODO:  write this test
        """
        pass

    def test_ambig_rc(self):
        """ test that we accurately are getting the reverse compliment of ambigous bases
        """
        res1 = clp.ambig_rc("[AT]ATCACTACTACTA[CA]TC", verbose=True)
        self.assertEqual(res1, "GA[TG]TAGTAGTAGTGAT[AT]")
        res2 = clp.ambig_rc("[AT]ATCACTAC[TACT]A[CA]TC", verbose=True)
        self.assertEqual(res2, "GA[TG]T[AGTA]GTAGTGAT[AT]")
        res3 = clp.ambig_rc("[AT]ATCACTAC[TACT]A[CA][TC]", verbose=True)
        self.assertEqual(res3, "[GA][TG]T[AGTA]GTAGTGAT[AT]")


    def test_ambig_to_regex(self):
        self.assertEqual(
            clp.ambig_to_regex("TRCCRCCAGTRCCAAAGACA"),
                               "T[AG]CC[AG]CCAGT[AG]CCAAAGACA",

        )
        self.assertEqual(
            clp.ambig_to_regex(
                "AGTTTTAYGCSSARTGCGAG"),
                "AGTTTTA[CT]GC[GC][GC]A[AG]TGCGAG"
        )

    def test_integration(self):
        """ execute ezclermont on connonical clermonted genome sequences
        These were selected from a lit review (primarily from clermont2013
        """
        A_ref = ["A", "NC_010468.1"]
        B1_ref = ["B1", "NC_011741.1"]
        E1_ref = ["E", "BA000007.2"]
        E2_ref = ["E", "AE005174.2"]
        D_ref = ["D", "NC_011751.1"]
        C_ref = ["C", "CP004009.1"]
        F1_ref = ["F", "NC_011750.1"]
        F2_ref = ["F", "LYBO00000000.1"]
        B2_ref = ["B2", "CU651637.1"]
        Ferg_ref = ["U/cryptic", "NC_011740.1"]
        G_ref = ["G", "GCF_900499975.1_921A_genomic.fna"]

        for ref in [A_ref, B1_ref, C_ref, E1_ref, E2_ref, D_ref, F1_ref,
                    F2_ref, B2_ref, Ferg_ref]:
            args = argparse.Namespace(
                contigs=os.path.join(os.path.dirname(__file__),
                                     "refs", ref[1] + ".fasta"),
                experiment_name=None,
                no_partial=False,
                min_length=500)
            result, profile = clp.main(args)
            self.assertEqual(
                ref[0], result)

    def test_integration2(self):
        """ execute ezclermont Big 6 and a O104 for good measure
        - https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5908948/
        - https://jcm.asm.org/content/50/9/2894
        - https://www.liebertpub.com/doi/10.1089/fpd.2018.2567
        """
        O26_ref  = ["B1", "CP028126.1"]
        O45_ref  = ["B1", "CP031916.1"]
        O103_ref = ["B1", "CP028112.1"]
        O111_ref = ["B1", "CP028117.1"]
        O121_ref = ["B1", "CP028117.1"]
        O145_ref = ["A",  "CP027105.1"]
        # no reference for this one
        O1452_ref = ["D", "CP028379.1"]
        O104_ref = ["B1", "CP003289.1"]
        for ref in [O26_ref, O45_ref, O103_ref, O111_ref, O121_ref, O145_ref, O1452_ref, O104_ref]:
            args = argparse.Namespace(
                contigs=os.path.join(os.path.dirname(__file__),
                                     "refs","big6", ref[1] + ".fasta"),
                experiment_name=None,
                no_partial=False,
                min_length=500)
            result, profile = clp.main(args)
            self.assertEqual(
                ref[0], result)

    def test_integration_noncoli(self):
        Crypt_ref = ["Cryptic", "AEMF01000001"]
        args = argparse.Namespace(
            contigs=os.path.join(os.path.dirname(__file__),
                                 "refs", Crypt_ref[1] + ".fasta"),
            experiment_name=None,
            no_partial=False,
            min_length=500)
        result, profile  = clp.main(args)
        self.assertEqual(result, "U/cryptic")


if __name__ == "__main__":
    unittest.main()
