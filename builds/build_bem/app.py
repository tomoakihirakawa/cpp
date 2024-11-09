from flask import Flask, render_template, request, redirect, url_for, flash, jsonify
from forms import SimulationForm
import os

app = Flask(__name__)
app.secret_key = 'your_secret_key'
UPLOAD_FOLDER = 'static/uploads'
app.config['UPLOAD_FOLDER'] = UPLOAD_FOLDER

# シミュレーション用オブジェクトのリスト
simulation_objects = []

@app.route('/', methods=['GET', 'POST'])
def index():
    form = SimulationForm()
    if form.validate_on_submit():
        # シミュレーションオブジェクトの属性とパスを収集
        for obj in simulation_objects:
            obj['role'] = request.form.get(f'role_{obj["name"]}')
        
        # 入力パラメータとオブジェクトを使用してシミュレーションを開始
        params = {
            'case': form.case.data,
            'dt': form.dt.data,
            'H': form.H.data,
            'outputdir': form.outputdir.data,
            'objects': simulation_objects  # オブジェクトリストを追加
        }
        
        # 実行後の処理を追加
        flash('Simulation started with custom objects!', 'success')
        return redirect(url_for('index'))

    return render_template('index.html', form=form, objects=simulation_objects)

@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'})
    
    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No selected file'})

    # ファイル保存とオブジェクトリストの更新
    filename = file.filename
    filepath = os.path.join(app.config['UPLOAD_FOLDER'], filename)
    file.save(filepath)

    # 3Dオブジェクトリストに追加
    simulation_objects.append({
        'name': filename,
        'path': filepath,
        'role': 'undefined'  # 初期状態は未定義
    })

    return jsonify({'filename': filename, 'path': filepath})

if __name__ == '__main__':
    app.run(debug=True)