#!/bin/sh

# ParaViewの実行ファイルへのパス
PARAVIEW_EXEC="/Applications/ParaView-6.0.0-RC3.app/Contents/MacOS/paraview"

# 引数が存在するかチェック
if [ -z "$1" ]; then
    echo "使用法: $0 <シミュレーションID>"
    exit 1
fi

# ★変更点1：引数を環境変数 'SIM_ID' に設定してエクスポート
export SIM_ID="$1"

echo "環境変数 SIM_ID を '${SIM_ID}' に設定しました。"
echo "ParaViewを起動します..."

# ★変更点2：ParaView起動時の引数をシンプルにする
$PARAVIEW_EXEC --script=visualize.py