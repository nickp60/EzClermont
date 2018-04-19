import os
import sys
from flask import Flask, render_template, request, session
from werkzeug import secure_filename
import tempfile
from argparse import Namespace

# strange way to do the imports, but hey
sys.path.append(os.path.join('..', 'clermontpcr'))
import clermontpcr

app = Flask(__name__)
#  set the max file size to 20mb.  If an ecoli fasta is
#    larger than this, then its probably not an e coli
app.config['MAX_CONTENT_LENGTH'] = 20 * 1024 * 1024

###  Yucky globals
# g.TMPDIR = None



def default_stream_factory(total_content_length, filename, content_type, content_length=None):
    """The stream factory that is used per default."""
    # if total_content_length > 1024 * 500:
    #    return TemporaryFile('wb+')
    # return BytesIO()
    return TemporaryFile('wb+')

def get_tmpfile_path():
    return os.path.join(tempfile.gettempdir(), "tmp.fasta")


@app.route("/", methods=['GET','POST'])
def index():
    file_content = [""]
    # return render_template("template.html", content=content) #
    if request.method == 'POST':
        print(request.form['submit'])
        if request.form['submit'] == "Upload File" :
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
            teaser = "\n".join(file_content[0:10])
            addn_lines = len(file_content) - 10
            if addn_lines < 0:
                addn_lines = 0
            header = "Here is the first bit of your file:"

            return render_template(
                'alt.html',
                header=header,
                content=teaser,
                nlines="...and %d more lines" % addn_lines
            )
        elif request.form['submit'] == "Run" :
            tmpfile = get_tmpfile_path()
            print(tmpfile)
            if os.path.isfile(tmpfile):
                print(tmpfile)
                return runcler(contigsfile=tmpfile)
            else:
                raise ValueError("clicky clicky clicky;  select a file first")
        else:
            raise ValueError
    else:
        return render_template("alt.html")

def runcler(contigsfile, ignore_control=False, partial=False):
    # prepare args
    args = Namespace(contigs=contigsfile,
                     ignore_control=ignore_control, partial=partial)
    results = clermontpcr.main(args)
    return render_template(
        "alt.html",
        results=results)

if __name__ == "__main__":
    session['tmpdir'] = tempfile.mkdtemp("clerpcr_tmps")
    app.run(debug=True, port=5957)
