#!/usr/bin/env python2
#-*- coding: utf-8 -*-
import os
import requests
from requests.auth import HTTPBasicAuth
import urllib.request, urllib.error, urllib.parse
import json
import base64
import sys
from urllib.error import HTTPError
import logging
import argparse
# You must have a valid API Token
API_TOKEN = os.getenv('ENTEROBASE_API_TOKEN', None)
API_USER = os.getenv('ENTEROBASE_API_USER', None)
if API_TOKEN is None:
    raise OSError("must have $ENTEROBASE_API_TOKEN variable set " +
                  "to your enterobase token!")
if API_USER is None:
    raise OSError("must have $ENTEROBASE_API_USER variable set " +
                  "to your enterobase Username!")
SERVER_ADDRESS = 'http://enterobase.warwick.ac.uk'

def get_args():  # pragma: no cover
    """#TODO:     for cli mods:
    http://stackoverflow.com/questions/18025646/
         python-argparse-conditional-requirements
    make this able to handle different library types such as two unpaired runs
    """
    parser = argparse.ArgumentParser(
        description="Give it a Enterbase barcode, it gives you the assembly ")

    parser.add_argument("-b", "--barcode_list",
                        dest='barcode',
                        action="store", type=str,
                        help="file containing single colmn of enterobase "+
                        "barcodes",
                        required=True)
    parser.add_argument("-o", "--outdir", dest='outdir', action="store",
                       help="output directory; " +
                       "default: %(default)s", default=os.getcwd(),
                        type=str, required=False)
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument("--max_tries", dest='max_tries', action="store",
                          help="how many times to retry a bad  https query",
                          type=int,  default=3)
    optional.add_argument("-m", "--manifest", dest='manifest', action="store",
                          help="file to append results to",
                          type=str, default="./fetchEnteroAssemblies.log")
    optional.add_argument("-d", "--database", dest='database', action="store",
                          help="name of the enterobase database",
                          type=str, default="ecoli")
    optional.add_argument("-v", "--verbose", dest='verbose', action="store_true",
                          help="log verbosely")
    args = parser.parse_args()
    return args


def process_barcode(args, barcode, outfile):
    address = str(SERVER_ADDRESS + '/api/v2.0/%s/straindata?barcode=%s' + \
              '&assembly_status=Assembled&only_fields=strain_name,d' + \
              'ownload_fasta_link') % (args.database, barcode)
    if args.verbose:
        print("API call: \n%s"  % address)
    with open(args.manifest, "a") as manif:
        try:
            code = 0
            this_try=1
            while code != 200 and this_try < args.max_tries:
                response = requests.get(address, auth=HTTPBasicAuth(API_USER,API_TOKEN))
                code = response.status_code
                this_try = this_try + 1
            if code != 200:
                return
            print(response.status_code)
            if response.status_code != 200:
                manif.write("%s\t%s\t%s\n" % (barcode, "UNKNOWN", "FAIL"))
                logging.error('%d %s' %(response.status_code, response.text))
                return 1
            data = response.json()
            if args.verbose:
                print(data)
            for record in data['straindata']:
                record_values = data['straindata'][record]
                if args.verbose:
                    print("Getting" + record_values['download_fasta_link'])
                r = requests.get(record_values['download_fasta_link'], auth=HTTPBasicAuth(API_USER, API_TOKEN), stream=True, allow_redirects=True)
                print(len(r.text))
                with open(outfile, 'w') as out_ass:
                    # for chunk in r.iter_content(chunk_size=128):
                    #     out_ass.write(chunk)
                    out_ass.write(r.text)
                manif.write("%s\t%s\t%s\n" % (barcode, record_values['strain_name'], "PASS"))
        except HTTPError as Response_error:
            manif.write("%s\t%s\t%s\n" % (barcode, "UNKNOWN", "FAIL"))
            logging.error('%d %s. <%s>\n Reason: %s' %(Response_error.code,
                                                       Response_error.msg,
                                                       Response_error.geturl(),
                                                       Response_error.read()))


if __name__ == "__main__":
    args = get_args()
    if not os.path.exists(args.outdir):
        os.mkdir(args.outdir)
    lenbarf = 0
    with open(args.barcode) as barf:
        for line in barf:
            lenbarf = lenbarf + 1
    with open(args.barcode) as barf:
        for idx, barcode in enumerate(barf):
            barcode = barcode.strip()
            outfile = os.path.join(args.outdir, '%s.fasta' % barcode)
            if os.path.exists(outfile) and os.path.getsize(outfile) > 0:
                print("skipping %s, item %i of %i" %(barcode, idx + 1, lenbarf))
            else:
                print("processing %s, item %i of %i" %(barcode, idx + 1, lenbarf))
                process_barcode(args=args, barcode=barcode, outfile=outfile)
