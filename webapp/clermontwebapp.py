import os
import sys
from flask import Flask, render_template, request, session
from werkzeug import secure_filename
import tempfile
from argparse import Namespace
from Bio import SeqIO

from cpcr import run as clermontpcr
# check that the import works
print(clermontpcr.PcrHit)

app = Flask(__name__)
# app.secret_key = 'crytographyisnotmystrongsuit'
app.secret_key = os.urandom(24)

#  set the max file size to 20mb.  If an ecoli fasta is
#    larger than this, then its probably not an e coli
app.config['MAX_CONTENT_LENGTH'] = 20 * 1024 * 1024


def default_stream_factory(total_content_length, filename, content_type, content_length=None):
    """The stream factory that is used per default."""
    # if total_content_length > 1024 * 500:
    #    return TemporaryFile('wb+')
    # return BytesIO()
    return TemporaryFile('wb+')

def get_tmpfile_path():
    return os.path.join(tempfile.gettempdir(), "tmp.fasta")


def is_fasta(filename):
    """ verbatim from https://stackoverflow.com/questions/44293407/
    works because SeqIO returns an empty generator rather than an error
    """
    with open(filename, "r") as handle:
        fasta = SeqIO.parse(handle, "fasta")
        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file


@app.route("/", methods=['GET','POST'])
def index():
    session["FINISHED"] = False
    session["LOADED"] = False
    file_content = [""]
    # return render_template("template.html", content=content) #
    if request.method == 'POST':
        print(request.form['submit'])
        if request.form['submit'] == "Get Clermont Phylotype!" :
            # for secure filenames. Read the documentation.
            file = request.files['myfile']
            filename = secure_filename(file.filename)
            tmpdir = tempfile.gettempdir()
            print(tmpdir)
            # tmpfile = os.path.join(tmpdir, filename)
            tmpfile = get_tmpfile_path()
            file.save(tmpfile)
            with open(tmpfile) as f:
                file_content = f.read().split("\n")
            teaser = "\n".join(file_content[0:7])
            addn_lines = len(file_content) - 7
            if addn_lines < 0:
                addn_lines = 0
            header = "Here is the first bit of your file:"
            session["LOADED"] = True

            ###   Now lets run the main function
            if os.path.isfile(tmpfile) and is_fasta(tmpfile):
                print(tmpfile)
                results, profile  = runcler(
                    contigsfile=tmpfile,
                    #ignore_control=request.form.get('allowcontrolfails'),
                    # partial=request.form.get('allowpartial')
                )
            else:
                results = "Unable to read file.  Are you sure its a valid fasta?"
                profile = ""
            ###
            return render_template(
                'alt.html',
                header=header,
                content=teaser,
                results=results,
                profile=profile,
                nlines="...and %d more lines" % addn_lines
            )
        else:
            pass
    else:
        return render_template("alt.html")

def runcler(contigsfile):
    # prepare args
    args = Namespace(contigs=contigsfile,
                     no_partial=False, min_length=500)
    try:
        results, profile  = clermontpcr.main(args)
    except Exception as e:
        if e is ImportError:
            results = "Deploy error!"
        else:
            print(e)
            results = str("Control fail!  Are you sure this is E. coli? " +
                          "Please contact the support team")
        profile = ""
    return (results, profile)

if __name__ == "__main__":
    app.run(debug=True, port=5957)
