#!/usr/bin/env python3
#-*- coding: utf-8 -*-
import re
import argparse
import sys
import os
import unittest
import itertools

from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# need this line for unittesting
sys.path.append(os.path.join('..', 'clermontpcr'))

class PcrHit(object):
    newid = itertools.count()

    def __init__(self, index=None,
                 allele_name=None,
                 template_orientation=None,
                 template_id=None,
                 true_hit=False,
                 F_start=None,
                 R_start=None,
                 F_end=None,
                 R_end=None,
                 F_hit=None,
                 R_hit=None,
                 allow_partial=False):
        # int: unique identifier for match
        self.index = next(PcrHit.newid)
        self.allele_name = allele_name
        self.template_orientation = template_orientation
        self.template_id = template_id
        self.true_hit = true_hit
        self.F_start = F_start
        self.F_end = F_end
        self.R_start = R_start
        self.R_end = R_end
        self.F_hit = F_hit
        self.R_hit = R_hit
        self.allow_partial = allow_partial
        self.parse_hits()

    def parse_hits(self):
        if self.F_start is not None and self.F_end is not None:
            self.F_hit = True
        if self.R_start is not None and self.R_end is not None:
            self.R_hit = True
        if self.F_hit and self.R_hit:
            self.true_hit = True
        elif (self.F_hit or self.R_hit) and self.allow_partial:
            self.true_hit = True
        # elif (self.F_hit or self.R_hit):
        else:
            self.true_hit = False


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="run a 'PCR' to get clermont types",
        add_help=False)  # to allow for custom help
    parser.add_argument("contigs", action="store",
                        help="FASTA formatted genome or set of contigs")

    # # taking a hint from http://stackoverflow.com/questions/24180527
    # requiredNamed = parser.add_argument_group('required named arguments')
    # requiredNamed.add_argument("-F", "--fastq1", dest='fastq1', action="store",
    #                            help="forward fastq reads, can be compressed",
    #                            type=str, default="", required=True)
    # # had to make this faux "optional" parse so that the named required ones
    # # above get listed first
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-p", "--partial", dest='partial',
                          action="store_true",
                          help="If scanning contigs, breaks between " +
                          "contigs could potentially contain your " +
                          "sequence of interest.  if --partial, partial " +
                          "matches that could be intereupted by contig " +
                          "breaks are reported",
                          default=False)
    optional.add_argument("-c", "--ignore_control", dest='ignore_control',
                          action="store_true",
                          help="if control PCR fails, continue",
                          default=False)
    # had to make this explicitly to call it a faux optional arg
    optional.add_argument("-h", "--help",
                          action="help", default=argparse.SUPPRESS,
                          help="Displays this help message")
    args = parser.parse_args()
    return args


def ambig_rc(rev_primer, verbose=False):
    """ This allows you to get the reverse complement when there are primers
    with ambiguities, as delineated by brackets for grepping.
    retruns the reverse compimented string
    """
    if verbose:
        print("starting seq")
        print(rev_primer)
    if(len(rev_primer.split("[")) > 1):
        # get the corrdinates of brackets in the reversed sequenced.  Note
        # that when the seqeuence is reversed "[" become "]", hence
        # "new forward brackets" refer to back brackets, etc
        new_for_bracks = [i for i, char in enumerate(rev_primer[::-1]) if
                          char == "]"]
        if verbose:
            print("new forward brackets")
            print(new_for_bracks)
        new_back_bracks = [i for i, char in enumerate(rev_primer[::-1]) if
                           char == "["]
        if verbose:
            print("new back brackets")
            print(new_back_bracks)
        assert len(new_for_bracks) == len(new_back_bracks), "unequal brackets!"
        clean_primer = rev_primer.replace("[", "").replace("]", "")
        rc_primer = str(SeqRecord(Seq(clean_primer).reverse_complement()).seq)
        if verbose:
            print("reverse comped sequence pre-bracket insertion")
            print(rc_primer)
        new_primer = ""
        # make up for the length of string after brackets are removed
        bufferX = "XX" * 2 * len(new_for_bracks)
        rc_primer = rc_primer + bufferX
        new_primer = rc_primer
        for idx in range(0, len(rc_primer)):
            if rc_primer[idx] is None:  # ignore end of string character
                break
            for i in range(0, len(new_for_bracks)):
                if idx == new_for_bracks[i]:
                    new_primer = new_primer[0:idx] + "[" + new_primer[idx: ]
                elif idx + 1 == new_back_bracks[i]:
                    new_primer = new_primer[
                        0:idx + 1] + "]" + new_primer[idx + 1: ]
                else:
                    pass
            if verbose:
                print("%i, %s" % (idx, new_primer))
        return(new_primer[:-len(bufferX)])
    else:
        return(str(SeqRecord(Seq(rev_primer).reverse_complement()).seq))


def get_matches(allele, seq_list, fwd_primer, rev_primer, expected_size,
                allow_partial=False, strand="+"):
    """given a seqence list  and regex compilations of your primers
    return the matches
    """
    # assert logger is not None, "must use logger!"
    assert strand in ["-", "+"], "strand must be either + or -"
    assert isinstance(seq_list[0], SeqRecord), "must submit list of SeqRecords"
    if strand == "+":
        fwd = re.compile(fwd_primer, re.IGNORECASE)
        rev = re.compile(ambig_rc(rev_primer), re.IGNORECASE)
    else:
        fwd = re.compile(rev_primer, re.IGNORECASE)
        rev = re.compile(ambig_rc(fwd_primer), re.IGNORECASE)
    matches = []
    for i in seq_list:
        coords_F = None
        coords_R = None
        try:
            coords_F = fwd.search(str(i.seq)).span()
            sys.stderr.write("F match!\n")
        except:
            pass
        try:
            coords_R = rev.search(str(i.seq)).span()
            sys.stderr.write("R match!\n")
        except:
            pass

        if coords_F is not None and coords_R is not None:
            sys.stderr.write("Match found on %s (%s)\n" % (i.id, strand))
            if abs(coords_R[1] - coords_F[0] - expected_size) < 20:
                matches.append(PcrHit(
                    allele_name=allele,
                    template_orientation=strand,
                    template_id=i.id,
                    F_start=coords_F[0],
                    R_start=coords_R[0],
                    F_end=coords_F[1],
                    R_end=coords_R[1],
                    allow_partial=allow_partial
                ))
            else:
                sys.stderr.write(
                    "Match found, but it is %ibp and it should only be %i\n" %
                    (coords_R[1] - coords_F[0], expected_size))
        elif coords_F is not None:
            if len(i.seq[coords_F[0]:]) > expected_size:
                sys.stderr.write("Possible match on %s (%s)\n" % (i.id, strand))
                matches.append(PcrHit(
                    allele_name=allele,
                    template_orientation=strand,
                    template_id=i.id,
                    F_start=coords_F[0],
                    R_start=None,
                    F_end=coords_F[1],
                    R_end=None,
                    allow_partial=allow_partial
                ))
            else:
                pass
        elif coords_R is not None:
            if not coords_R[0] < expected_size:
                sys.stderr.write("Possible match on %s (%s)\n" % (i.id, strand))
                matches.append(PcrHit(
                    allele_name=allele,
                    template_orientation=strand,
                    template_id=i.id,
                    R_start=coords_R[0],
                    F_start=None,
                    R_end=coords_R[1],
                    F_end=None,
                    allow_partial=allow_partial
                ))
        else:
            # sys.stderr.write("No hits on %s" % i.id)
            pass
    return(matches)


def interpret_hits(arpA, chu, TspE4, yjaA):
    if arpA:
        if chu:
            if TspE4 and yjaA:
                result = "U"
            elif not TspE4 and not yjaA:
                result = "E/D"
            elif TspE4:
                result = "E/D"
            else:
                assert yjaA, "error interpretting results!"
                result = "E/cryptic"
        else:
            if TspE4 and yjaA:
                result = "U"
            elif not TspE4 and not yjaA:
                result = "A"
            elif TspE4:
                result = "B1"
            else:
                assert yjaA, "error interpretting results"
                result = "A/C"
    else:
        if chu:
            if yjaA or TspE4:
                result = "B2"
            else:
                result = "F"
        else:
            if yjaA:
                result = "cryptic"
            else:
                result = "U"
    return(result)


def refine_hits(hit, c_primers, e_primers, allow_partial, seqs):
    if hit == "E/D":
        e_primers["arpA_e"], report_string = run_primer_pair(
            seqs=seqs, allele="arpA_e",
            vals=e_primers["arpA_e"],
            allow_partial=allow_partial)
        if e_primers["arpA_e"][3]:
            return "E"
        else:
            return "D"
    elif hit == "E/cryptic":
        e_primers["arpA_e"], report_string = run_primer_pair(
            seqs=seqs, allele="arpA_e",
            vals=e_primers["arpA_e"],
            allow_partial=allow_partial)
        if e_primers["arpA_e"][3]:
            return "E"
        else:
            return "cryptic"
    elif hit == "A/C":
        c_primers["trpA_c"], report_string = run_primer_pair(
            seqs=seqs, allele="trpA_c",
            vals=c_primers["trpA_c"],
            allow_partial=allow_partial)
        if c_primers["trpA_c"][3]:
            return "C"
        else:
            return "A"
    else:
        return hit


def run_primer_pair(seqs, allele, vals, allow_partial):
    """ expects a vals with the following structure:
    [forward, reverse, length]
    returns a list with True or False appended to value list,
    and a report string:
     ([forward, reverse, length, True], "allele: +" )

    """
    # forward match
    fwd_matches = get_matches(allele=allele,
                              seq_list=seqs,
                              fwd_primer=vals[0],
                              rev_primer=vals[1],
                              allow_partial=allow_partial,
                              expected_size=vals[2],
                              strand='+')
    # reverse compliment
    rev_matches = get_matches(allele=allele,
                              seq_list=seqs,
                              fwd_primer=vals[0],
                              rev_primer=vals[1],
                              allow_partial=allow_partial,
                              expected_size=vals[2],
                              strand='-')
    if (
            any([x.true_hit for x in fwd_matches]) or
            any([x.true_hit for x in rev_matches])):
        # sys.stderr.write("%s: +" % key)
        profile = "{0}: +".format(allele)
        vals.append(True)
    else:
        profile = "{0}: -".format(allele)
        vals.append(False)
    return(vals, profile)


def main(args):
    # chuA
    chuA_1b = "ATGGTACCGGACGAACCAAC"
    chuA_2  = "TGCCGCCAGTACCAAAGACA"
    # yjaH
    yjaA_1b = "CAAACGTGAAGTGTCAGGAG"
    yjaA_2b = "AATGCGTTCCTCAACCTGTG"
    # TspE4.C2
    TspE4_C2 = "CACTATTCGTAAGGTCATCC"
    TspE4C2_2b = "AGTTTATCGCTGCGGGTCGC"
    # arpA
    # AceK_f = "AACGCTATTCGCCAGCTTGC"
    AceK_f = "AA[CT]GC[TC]ATTCGCCAGCTTGC"
    ArpA1_r = "TCTCCCCATACCGTACGCTA"
    # arpA, for group e
    ArpAgpE_f = "GATTCCATCTTGTCAAAATATGCC"
    ArpAgpE_r = "GAAAAGAAAAAGAATTCCCAAGAG"
    # trpA, for group c
    trpAgpC_1 = "AGTTTTATGCCCAGTGCGAG"
    trpAgpC_2 = "TCTGCGCCGGTCACGCCC"
    # trpA, for control
    trpBA_f = "CGGCGATAAAGACATCTTCAC"
    trpBA_r = "GCAACGCGGCCTGGCGGAAG"

    quad_primers = {"chu": [chuA_1b, chuA_2, 288],
                    "yjaA": [yjaA_1b, yjaA_2b, 211],
                    "TspE4": [TspE4_C2, TspE4C2_2b, 152],
                    "arpA": [AceK_f, ArpA1_r, 400]}
    c_primers = {"trpA_c": [trpAgpC_1, trpAgpC_2, 219]}
    e_primers = {"arpA_e": [ArpAgpE_f, ArpAgpE_r, 301]}

    controls = {"trpBA_control": [trpBA_f, trpBA_r, 489]}
    sys.stderr.write("Reading in sequence(s)\n")
    with open(args.contigs, 'r') as fasta:
        seqs = list(SeqIO.parse(fasta, 'fasta'))

    # Start with the Control primers before trying anythong else
    sys.stderr.write("Running Control PCR\n")
    controls["trpBA_control"], report_string = run_primer_pair(
        seqs=seqs, allele="trpBA_control",
        vals=controls["trpBA_control"],
        allow_partial=args.partial)
    sys.stderr.write(report_string + "\n")
    if not controls["trpBA_control"][3]:
        if args.ignore_control:
            sys.stderr.write(
                "No matches found for control PCR, but continuing analysis\n")
        else:
            sys.stderr.write("No matches found for control PCR.  Exiting\n")
            sys.stdout.write(
                "{0}\t{1}\n".format(
                    os.path.splitext(os.path.basename(args.contigs))[0],
                    "control_fail"
                ))
            sys.exit(1)
    else:
        pass
    # run Clermont Typing
    sys.stderr.write("Running Quadriplex PCR\n")
    profile = ""
    for key, val in sorted(quad_primers.items()):
        sys.stderr.write("Scanning %s\n" % key)
        quad_primers[key], report_string = run_primer_pair(
            seqs=seqs, allele=key,
            vals=val,
            allow_partial=args.partial)
        profile = profile + report_string + "\n"
    sys.stderr.write("\n-------- Results --------\n")
    sys.stderr.write(profile)
    sys.stderr.write("\n--------   ---   --------\n")
    Clermont_type = interpret_hits(arpA=quad_primers['arpA'][3],
                                   chu=quad_primers['chu'][3],
                                   TspE4=quad_primers['TspE4'][3],
                                   yjaA=quad_primers['yjaA'][3])
    Clermont_type = refine_hits(hit=Clermont_type,
                                c_primers=c_primers,
                                e_primers=e_primers,
                                allow_partial=args.partial,
                                seqs=seqs)
    sys.stderr.write("Clermont type: %s\n" % Clermont_type)
    sys.stderr.write("-------------------------\n")
    # This should be the only thing being written to stdout
    sys.stdout.write(
        "{0}\t{1}\n".format(
            os.path.splitext(os.path.basename(args.contigs))[0],
            Clermont_type
        ))
    return(Clermont_type)


class clermontTestCase(unittest.TestCase):
    """
    """
    def test_interpret(self):

        ref = ['A', 'A/C', 'B1', 'A/C', 'E/D', 'E/D', 'E/cryptic', 'E/D',
               'E/D', 'F', 'B2', 'B2', 'B2', 'cryptic', 'E/cryptic', 'U']
        test = []
        # if True:
        # A's
        test.append(
            interpret_hits(arpA=True, chu=False, yjaA=False, TspE4=False))
        test.append(
            interpret_hits(arpA=True, chu=False, yjaA=True, TspE4=False))
        # B1
        test.append(
            interpret_hits(arpA=True, chu=False, yjaA=False, TspE4=True))
        # C
        test.append(
            interpret_hits(arpA=True, chu=False, yjaA=True, TspE4=False))
        # E
        test.append(
            interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=False))
        test.append(
            interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=True))
        test.append(
            interpret_hits(arpA=True, chu=True, yjaA=True, TspE4=False))
        # D
        test.append(
            interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=False))
        test.append(
            interpret_hits(arpA=True, chu=True, yjaA=False, TspE4=True))
        # F
        test.append(
            interpret_hits(arpA=False, chu=True, yjaA=False, TspE4=False))
        # B2
        test.append(
            interpret_hits(arpA=False, chu=True, yjaA=True, TspE4=False))
        test.append(
            interpret_hits(arpA=False, chu=True, yjaA=False, TspE4=True))
        test.append(
            interpret_hits(arpA=False, chu=True, yjaA=True, TspE4=True))
        # cryptic
        test.append(
            interpret_hits(arpA=False, chu=False, yjaA=True, TspE4=False))
        test.append(
            interpret_hits(arpA=True, chu=True, yjaA=True, TspE4=False))
        # unknown
        test.append(
            interpret_hits(arpA=False, chu=False, yjaA=False, TspE4=False))
        print(ref)
        print(test)
        self.assertEqual(ref, test)

    def test_get_matches(self):
        test_seq = [SeqRecord(
            Seq("GCACAGTCGATCAAAATTTTTGCAGTCGACTGGACTGACTGTCGGATCTCAGTCAT"))]
        test_fwd = "GCACAG"
        fwd = re.compile(test_fwd, re.IGNORECASE)
        test_rev = "TGACTG"
        self.assertEqual(fwd.search(str(test_seq[0].seq)).span(), (0, 6))
        forward_control_matches = get_matches(
            allele="test",
            seq_list=test_seq,
            fwd_primer=test_fwd,
            rev_primer=test_rev,
            allow_partial=True,
            expected_size=60,
            strand='+')
        self.assertEqual(forward_control_matches[0].R_end, 55)

    def test_refine_hits(self):
        pass

    def test_ambig_rc(self):
        res1 = ambig_rc("[AT]ATCACTACTACTA[CA]TC", verbose=True)
        self.assertEqual(res1, "GA[TG]TAGTAGTAGTGAT[AT]")
        res2 = ambig_rc("[AT]ATCACTAC[TACT]A[CA]TC", verbose=True)
        self.assertEqual(res2, "GA[TG]T[AGTA]GTAGTGAT[AT]")
        res3 = ambig_rc("[AT]ATCACTAC[TACT]A[CA][TC]", verbose=True)
        self.assertEqual(res3, "[GA][TG]T[AGTA]GTAGTGAT[AT]")

    def test_integration(self):
        A_ref = ["A", "NC_010468.1"]
        B1_ref = ["B1", "NC_011741.1"]
        E1_ref = ["E", "BA000007.2"]
        E2_ref = ["E", "AE005174.2"]
        D_ref = ["D", "NC_011751.1"]
        F1_ref = ["F", "NC_011750.1"]
        F2_ref = ["F", "LYBO00000000.1"]
        B2_ref = ["B2", "CU651637.1"]
        Ferg_ref = ["U", "NC_011740.1"]

        for ref in [A_ref, B1_ref, E1_ref, E2_ref, D_ref, F1_ref,
                    F2_ref, B2_ref, Ferg_ref]:
            args = argparse.Namespace(
                contigs=os.path.join(os.path.dirname(__file__),
                                     "refs", ref[1] + ".fasta"),
                partial=False,
                ignore_control=False)
            result = main(args)
            self.assertEqual(
                ref[0], result)

    def test_integration_noncoli(self):
        Crypt_ref = ["Cryptic", "AEMF01000001"]
        args = argparse.Namespace(
            contigs=os.path.join(os.path.dirname(__file__),
                                 "refs", Crypt_ref[1] + ".fasta"),
            partial=False,
            ignore_control=False)
        with self.assertRaises(SystemExit):
            result = main(args)


if __name__ == "__main__":
    args = get_args()
    main(args)
