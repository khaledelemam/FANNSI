from flask import Flask,render_template, request, redirect,url_for, session, json
from search import runSearch
import os

app = Flask(__name__)
app.config["DEBUG"] = True
app.config['JSON_SORT_KEYS'] = False


@app.route("/", methods=['POST','GET'])
def input_seq():
    if request.method == "POST":
        seq = request.form.get("seq")
        threshold = request.form.get("threshold")
        if threshold == "":
            threshold = 10
        dct = runSearch(seq,int(threshold))
        res = json.dumps(dct)
        values = json.loads(res)
        # return redirect(url_for('result', values = res))
        return render_template('result.html', result = values)
    return render_template('index.html')

# @app.route("/result", methods = ['GET'])
# def result():
#     values = json.loads(request.args.get("values"))
#     return render_template('result.html', result= values)

 


if __name__ == '__main__':
    app.secret_key = 'secres_key'
    port = int(os.environ.get('PORT', 80))
    app.run(host='0.0.0.0', port=port, debug=True)