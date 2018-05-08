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
        else:
            self.true_hit = False


def get_args():  # pragma: no cover
    parser = argparse.ArgumentParser(
        description="run a 'PCR' to get clermont types",
        add_help=False)  # to allow for custom help
    parser.add_argument("contigs", action="store",
                        help="FASTA formatted genome or set of contigs")

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("-m", "--min_length", dest='min_length',
                          help="minimum contig length to consider." +
                          "default: %(default)s",
                          default=500)
    optional.add_argument("-n", "--no_partial", dest='no_partial',
                          action="store_true",
                          help="If scanning contigs, breaks between " +
                          "contigs could potentially contain your " +
                          "sequence of interest.  if --np_partial, partial " +
                          "matches that could be intereupted by contig " +
                          "breaks are reported; default " +
                          "behaviour is to consider partial hits if the " +
                          "sequence is an assembly of more than 4 sequences" +
                          "(ie, no partial matches for complete genomes, " +
                          "allowing for 1 chromasome and several plasmids)",
                          default=False)
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
            sys.stderr.write("  F match! on %s \n" % i.id)
        except:
            pass
        try:
            coords_R = rev.search(str(i.seq)).span()
            sys.stderr.write("  R match! on %s \n" % i.id)
        except:
            pass

        if coords_F is not None and coords_R is not None:
            sys.stderr.write("  Match found on %s (%s)\n" % (i.id, strand))
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
                    "  Match found, but it is %ibp and it should only be %i\n" %
                    (coords_R[1] - coords_F[0], expected_size))
        elif coords_F is not None:
            if (
                    (strand == "+" and len(i.seq) - coords_F[0] < expected_size) or
                    (strand == "-" and coords_F[0] < expected_size)
            ):
                sys.stderr.write("  Possible match on %s (%s) %d\n" % (i.id, strand, coords_F[0]))
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
                sys.stderr.write("  Only forward primer hit found\n")
                pass
        elif coords_R is not None:
            if (
                    (strand == "+" and len(i.seq) - coords_R[0] < expected_size) or
                    (strand == "-" and coords_R[0] < expected_size)
            ):
                sys.stderr.write("  Possible match on %s (%s)\n" % (i.id, strand))
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
                sys.stderr.write("  Possible match on %s (%s) %d\n" % (i.id, strand, coords_R[0]))
            else:
                sys.stderr.write("  Only reverse primer hit found\n")
        else:
            # sys.stderr.write("No hits on %s" % i.id)
            pass
    return(matches)


def interpret_hits(arpA, chu, yjaA, TspE4):
    # {result: [[arpA,chuA, yjaA, TspE4.C2],
    #           [arpA,chuA, yjaA, TspE4.C2]], etc},
    # also, screw pep8 spacing, I need this to look pretty
    clermont_key = \
        {
            "A":         [[1, 0, 0, 0]],
            "B1":        [[1, 0, 0, 1]],
            "F":         [[0, 1, 0, 0]],
            "B2":        [[0, 1, 1, 0],
                          [0, 1, 1, 1],
                          [0, 1, 0, 1]],
            "A/C":       [[1, 0, 1, 0]],
            "D/E":       [[1, 1, 0, 0],
                          [1, 1, 0, 1]],
            "E/cryptic": [[1, 1, 1, 0]],
            "cryptic":   [[0, 0, 1, 0]],
            # add this in when we learn the frag sequence, see note below
            "U/cryptic": [[0, 0, 0, 0]],
            # "U":         [[0, 0, 0, 0],
            "U":         [[0, 0, 0, 1],
                          [0, 0, 1, 1],
                          [1, 0, 1, 1],
                          [1, 1, 1, 1]]
        }
    isolate_profile = [arpA, chu, yjaA, TspE4]
    for key, value in clermont_key.items():
        for profile in value:
            if profile == isolate_profile:
                return key
    raise ValueError("results profile could not be interpretted! Exiting")


def refine_hits(hit, c_primers, e_primers, cryptic_chu_primers, EC_control_fail,
                allow_partial, seqs):

    if hit == "D/E":
        sys.stderr.write("Clermont type is D/E; running ArpAgpE primers\n")
        e_primers["arpA_e"], report_string = run_primer_pair(
            seqs=seqs, allele="arpA_e",
            vals=e_primers["arpA_e"],
            allow_partial=allow_partial)
        if e_primers["arpA_e"][3]:
            return "E"
        else:
            return "D"
    elif hit == "E/cryptic":
        sys.stderr.write("Clermont type is E/cryptic; running ArpAgpE primers\n")
        e_primers["arpA_e"], report_string = run_primer_pair(
            seqs=seqs, allele="arpA_e",
            vals=e_primers["arpA_e"],
            allow_partial=allow_partial)
        if e_primers["arpA_e"][3]:
            return "E"
        else:
            return "cryptic"
    elif hit == "A/C":
        sys.stderr.write("Clermont type is A/C; running the trpAgpC primers\n")
        c_primers["trpA_c"], report_string = run_primer_pair(
            seqs=seqs, allele="trpA_c",
            vals=c_primers["trpA_c"],
            allow_partial=allow_partial)
        if c_primers["trpA_c"][3]:
            return "C"
        else:
            return "A"
    # this is to resolve instances where a large, 476bp chuA fragment can
    # differentiate between cryptic clades and the unknown.
    # "Instead yield a PCR product that is 476 bp in size (Fig. 1).
    # This product is a fragment of the chuA gene that is a consequence of
    # an 11 bp match between the primer AceK.f and the chuA gene and the
    # presence of the chuA.2 primer. This has been confirmed by sequencing of
    # the PCR product (data not shown)."

    #  I dont have a test case for this, so we will just report abiguously
    # elif hit == "U/cryptic":
    #     cryptic_chu_primers["476_chu"], report_string = run_primer_pair(
    #         seqs=seqs, allele="476_chu",
    #         vals=cryptic_chu_primers["476_chu"],
    #         allow_partial=allow_partial)
    #     if cryptic_chu_primers["476_chu"][3]:
    #         return "cryptic"
    #     else:
    #         return "U"
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


def main(args=None):
    if args is None:
        args = get_args()
    ####### Quadriplex PCR ########
    # chuA
    chuA_1b = "ATGGTACCGGACGAACCAAC"
    chuA_2  = "TGCC[GA]CCAGTACCAAAGACA"
    # yjaH
    yjaA_1b = "CAAACGTGAAGTGTCAGGAG"
    yjaA_2b = "AAT[GA]CGTTCCTCAACCTGTG"
    # TspE4.C2
    # TspE4_C2 = "CACTATTCGTAAGGTCATC[CG]"
    TspE4C2_1b =   "CACTATTCGTAAGG[TC]CATCC"
    TspE4C2_2b = "AGTTTATCGCTGCGGGTCGC"
    # arpA
    AceK_f =  "AA[CT]GC[TC]ATTCGCCAGCTTGC"
    ArpA1_r = "TCTCC[CA]CATA[CT][CA]G[TC]ACGCTA"
    # ArpA1_r = "TCTCCCCATACCGTACGCTA"
    #################################
    # arpA, for group e
    ArpAgpE_f = "GAT[GT]CCAT[CT]TTGTC[AG]AAATATGCC"
    ArpAgpE_r = "GAAAA[GT]AAAAAGA[AC]TT[CT][CAT]CAAGAG"
    # trpA, for group c
    trpAgpC_1 = "AGTTTTATGCC[CG]A[GA]TGCGAG"
    trpAgpC_2 = "TC[TA]GC[GT]C[CT]GGTCACGCCC"
    # trpAgpC_2 = "TC[TA]GC[GT]C[CT]GGTCA[CT][GA]CC[CT]"
    # trpA, for control
    trpBA_f = "CGGCGATAAAGACAT[CT]TTCAC"
    trpBA_r = "GCAACGCGGC[CT]TGGCGGAAG"

    quad_primers = {"chu": [chuA_1b, chuA_2, 288],
                    "yjaA": [yjaA_1b, yjaA_2b, 211],
                    "TspE4": [TspE4C2_1b, TspE4C2_2b, 152],
                    "arpA": [AceK_f, ArpA1_r, 400]}
    c_primers = {"trpA_c": [trpAgpC_1, trpAgpC_2, 219]}
    e_primers = {"arpA_e": [ArpAgpE_f, ArpAgpE_r, 301]}
    cryptic_chu_primers = {"476_chu": [AceK_f, chuA_2, 476]}

    controls = {"trpBA_control": [trpBA_f, trpBA_r, 489]}
    sys.stderr.write("Reading in sequence(s) from %s\n" %
                     os.path.basename(args.contigs))
    seqs = []
    rejected_contigs = 0
    with open(args.contigs, 'r') as fasta:
        for seq in SeqIO.parse(fasta,"fasta"):
            if len(seq.seq) > args.min_length:
                seqs.append(seq)
            else:
                rejected_contigs = rejected_contigs + 1
    if rejected_contigs > 0:
        sys.stderr.write(str(
            "Ignoring %d contigs less than the set " +
            "minimum contig length (%d)\n") %
                         (rejected_contigs, args.min_length))
    allow_partial = False
    if len(seqs) > 4 and not args.no_partial:
        allow_partial = True
    else:
        sys.stderr.write("rejecting potetial partial matches\n")

    # Start with the Control primers before trying anything else
    # these are used primarily for differentiating the C/E groups
    sys.stderr.write("Running Control PCR\n")
    controls["trpBA_control"], report_string = run_primer_pair(
        seqs=seqs, allele="trpBA_control",
        vals=controls["trpBA_control"],
        allow_partial=allow_partial)
    sys.stderr.write(report_string + "\n")
    EC_control_fail = False
    if not controls["trpBA_control"][3]:
        sys.stderr.write(
                "No matches found for control PCR, but continuing analysis\n")
        EC_control_fail = True
    #     else:
    #         sys.stderr.write("No matches found for control PCR.  Exiting\n")
    #         sys.stdout.write(
    #             "{0}\t{1}\n".format(
    #                 os.path.splitext(os.path.basename(args.contigs))[0],
    #                 "control_fail"
    #             ))
    #         raise ValueError("Control (trpBA) failed")
    # else:
    #     pass
    # run Clermont Typing
    sys.stderr.write("Running Quadriplex PCR\n")
    profile = ""
    for key, val in sorted(quad_primers.items()):
        sys.stderr.write("Scanning %s\n" % key)
        quad_primers[key], report_string = run_primer_pair(
            seqs=seqs, allele=key,
            vals=val,
            allow_partial=allow_partial)
        profile = profile + report_string + "\n"
    sys.stderr.write("\n-------- Results --------\n")
    sys.stderr.write(profile)
    sys.stderr.write("\n--------   ---   --------\n")
    Clermont_type = interpret_hits(arpA=quad_primers['arpA'][3],
                                   chu=quad_primers['chu'][3],
                                   TspE4=quad_primers['TspE4'][3],
                                   yjaA=quad_primers['yjaA'][3])
    refine = True
    if Clermont_type in ["A/C", "D/E", "E/cryptic"]:
        # All these reactions use the internal control.  If that fails,
        # we call it "EC_control_fail"
        if EC_control_fail:
            sys.stderr.write("No matches found for C/E control PCR. Exiting\n")
            Clermont_type = "EC_control_fail"
            refine = False
    if refine:
        Clermont_type = refine_hits(
            hit=Clermont_type,
            c_primers=c_primers,
            e_primers=e_primers,
            cryptic_chu_primers=cryptic_chu_primers,
            allow_partial=allow_partial,
            EC_control_fail=EC_control_fail,
            seqs=seqs)
    sys.stderr.write("Clermont type: %s\n" % Clermont_type)
    sys.stderr.write("-------------------------\n")
    # This should be the only thing being written to stdout
    sys.stdout.write(
        "{0}\t{1}\n".format(
            os.path.splitext(os.path.basename(args.contigs))[0],
            Clermont_type
        ))
    return(Clermont_type, profile)


if __name__ == "__main__":
    args = get_args()
    main(args)
