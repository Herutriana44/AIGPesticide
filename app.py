from flask import Flask, render_template, request
from GAN import GAN
import pandas as pd
import os

app = Flask(__name__)

@app.route('/')
def dashboard():
    return render_template('dashboard.html')

@app.route('/input', methods=['GET', 'POST'])
def input_page():
    data = pd.DataFrame(columns=['Substance', 'LD50', 'Time'])
    if not os.path.exists('static/data/CompoundPersticide.csv'):
        df = pd.DataFrame(columns=['Substance', 'LD50', 'Time'])
        df.to_csv('static/data/CompoundPersticide.csv', index=False)

    dt = pd.read_csv('static/data/CompoundPersticide.csv')
    return render_template('input.html', row=dt, panjang=len(dt))

@app.route('/result', methods=['GET', 'POST'])
def result_page():
    if request.method == 'POST':
        num_molecules = int(request.form['num_molecules'])
        result = GAN.main(num_molecules)
        data_result = result
        return render_template('result.html', data=data_result, panjang=len(data_result))

if __name__ == '__main__':
    app.run(debug=True)