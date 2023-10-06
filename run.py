#coding: utf8

from werkzeug.utils import secure_filename
from flask import Flask,render_template,jsonify,request, send_file
import time
import os
import base64
import json

import lcmsseq
from params import CmdParams

app = Flask(__name__)
UPLOAD_FOLDER='upload'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER
basedir = os.path.abspath(os.path.dirname(__file__))
ALLOWED_EXTENSIONS = set(['csv'])

def allowed_file(filename):
    return '.' in filename and filename.rsplit('.',1)[1] in ALLOWED_EXTENSIONS

@app.route('/seq')
def upload_sample():
    return render_template('./upload.html')

@app.route('/seq/upload', methods=['POST'],strict_slashes=False)
def api_upload():
    file_dir = os.path.join(basedir,app.config['UPLOAD_FOLDER'])
    if not os.path.exists(file_dir):
        os.makedirs(file_dir)
    print(request.form)
    f = request.files['myfile']  # get file name from the POSTed form

    params = CmdParams(request.form)
    print(params)

    if f and allowed_file(f.filename):  # check the filetype 
        fname = secure_filename(f.filename)
        print(fname)
        ext = fname.rsplit('.',1)[1]  # the suffix
        unix_time = int(time.time())
        new_filename = str(unix_time)+'.'+ext  # replace the name
        f.save(os.path.join(file_dir,new_filename))  # save to 'upload' folder

        lcmsseq.read_params('default.cfg')
        final_file = os.path.join(file_dir, new_filename)
        ret, img_file = lcmsseq.process(final_file, fname, params)

        prety_list = json.dumps(ret, indent=4, sort_keys=True)
        return render_template('seq.html', seq_list=ret, img_file=img_file)
    else:
        return jsonify({"errno":1001,"errmsg":"failed"})

if __name__ == '__main__':
    app.run(host='0.0.0.0', port=8999, debug=True)
