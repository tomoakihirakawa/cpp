import paraview
from paraview.simple import *

# サンプルデータセットの読み込み
cone = Cone()

# ビューの作成
renderView = GetActiveViewOrCreate('RenderView')

# コーンの表示設定
coneDisplay = Show(cone, renderView)

# 色と背景の設定
coneDisplay.Representation = 'Surface'
renderView.Background = [0.1, 0.2, 0.3]

# カメラの位置設定
renderView.ResetCamera()

# レンダリングウィンドウのサイズ設定
renderView.ViewSize = [800, 600]

# レンダリングと画像の保存
SaveScreenshot('/path/to/save/image.png', renderView)

# ParaViewのGUIを使用している場合は、以下のコマンドでGUIを更新できます。
# しかし、スクリプトから実行される場合は必要ありません。
# Interact()
