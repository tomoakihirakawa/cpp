import os
import json
import pyvista as pv
from trame.app import get_server
from trame.ui.vuetify import SinglePageLayout
from trame.widgets import vuetify, vtk, html

# --- 初期設定 ---
SCRIPT_DIRECTORY = os.path.dirname(os.path.abspath(__file__))
from input_generator_logic import get_case_data
server = get_server(client_type="vue2")
state, ctrl = server.state, server.controller

# --- PyVista プロッター設定 ---
plotter = pv.Plotter(off_screen=True)
actors = {}
last_picked_actor = None
original_colors = {}
# `remoteId` とオブジェクト名のマッピングを保存するための辞書
remote_id_map = {}
# --- アプリケーションのロジック ---

@server.state.change("selected_case")
def update_application(selected_case, **kwargs):
    if not selected_case:
        return
    global last_picked_actor
    last_picked_actor = None
    state.selected_object_info = "オブジェクトをクリックして選択"
    state.active_panel = None
    settings, objects = get_case_data(selected_case, base_dir=SCRIPT_DIRECTORY)
    state.settings, state.objects = (settings, objects) if settings else ({}, [])
    update_ui_and_plot()
    float_obj = next((obj for obj in state.objects if "float" in obj.get("name", "").lower()), None)
    if float_obj and "position" in float_obj:
        pos = float_obj["position"]
        state.float_pos_x, state.float_pos_y, state.float_pos_z = pos[0], pos[1], pos[2]
    else:
        state.float_pos_x, state.float_pos_y, state.float_pos_z = 0, 0, 0

@server.trigger("open_edit_dialog")
def open_edit_dialog(index):
    print(f"DEBUG: open_edit_dialog called with index: {index}")
    try:
        if index < len(state.objects):
            state.edited_object = dict(state.objects[index])
            state.edited_object_index = index
            state.edit_dialog = True
            print(f"DEBUG: Dialog opened for object: {state.edited_object.get('name', 'Unknown')}")
        else:
            print(f"DEBUG: Invalid index {index}, objects length: {len(state.objects)}")
    except Exception as e:
        print(f"DEBUG: Error in open_edit_dialog: {e}")

@server.trigger("save_changes")
def save_changes():
    if state.edited_object_index is not None:
        numeric_props = ["mass", "isFixed"]
        for prop in numeric_props:
            if prop in state.edited_object:
                try:
                    val = float(state.edited_object[prop])
                    state.edited_object[prop] = int(val) if val.is_integer() else val
                except (ValueError, TypeError):
                    pass
        state.objects[state.edited_object_index] = state.edited_object
        state.dirty_state("objects")
    close_edit_dialog()
    update_ui_and_plot()

@server.trigger("close_edit_dialog")
def close_edit_dialog():
    state.edit_dialog = False
    state.edited_object = {}
    state.edited_object_index = None

# https://github.com/pyvista/pyvista/discussions/5132
# https://github.com/Kitware/trame-vtk/blob/master/examples/validation/PickingLocalView.py

@server.trigger("on_click")
def on_view_click(event=None):
    print(f"DEBUG: on_view_click called with event: {event}")
    global last_picked_actor
    
    # イベントが有効かつ remoteId が存在する場合にのみ処理
    if not event or 'remoteId' not in event or not event['remoteId']:
        # オブジェクトがクリックされなかった場合の処理
        if last_picked_actor:
            last_picked_actor.prop.color = original_colors.get(last_picked_actor, "white")
            last_picked_actor = None
        state.selected_object_info = "背景がクリックされました"
        state.active_panel = None
    else:
        # 以前に選択されたアクターの色を戻す
        if last_picked_actor:
            last_picked_actor.prop.color = original_colors.get(last_picked_actor, "white")
        
        # remoteId からオブジェクト名を取得
        remote_id = event['remoteId']
        name = remote_id_map.get(remote_id)
        
        if name:
            picked_actor = actors.get(name)
            if picked_actor:
                picked_actor.prop.color = "red"
                last_picked_actor = picked_actor
                state.selected_object_info = f"選択中: {name}"
                for i, obj in enumerate(state.objects):
                    if obj.get("name") == name:
                        state.active_panel = i
                        break
        
    if hasattr(ctrl, 'view_update') and callable(ctrl.view_update):
        ctrl.view_update()


@server.state.change("float_pos_x", "float_pos_y", "float_pos_z")
def update_float_position(float_pos_x, float_pos_y, float_pos_z, **kwargs):
    float_key = next((key for key in actors if "float" in key or "box" in key), None)
    if float_key:
        actors[float_key].position = (float_pos_x, float_pos_y, float_pos_z)
        if hasattr(ctrl, 'view_update') and callable(ctrl.view_update):
            ctrl.view_update()

# 新しいアプローチ: ボタンクリックを state change で処理
@server.state.change("edit_button_clicked")
def handle_edit_button_click(edit_button_clicked, **kwargs):
    if edit_button_clicked is not None:
        print(f"DEBUG: Edit button clicked for index: {edit_button_clicked}")
        if edit_button_clicked >= 0:
            open_edit_dialog(edit_button_clicked)
        elif edit_button_clicked == -1:
            # 保存ボタンの処理
            save_changes()
        # クリック処理後にリセット
        state.edit_button_clicked = None

# def update_ui_and_plot(view=None, **kwargs):
#     print("DEBUG: update_ui_and_plot called")
#     settings_widgets = []
#     if state.settings:
#         for key, value in state.settings.items():
#             settings_widgets.append(vuetify.VCol(vuetify.VTextField(label=key, value=value, readonly=True, dense=True, outlined=True), cols=6, py_1=True))
    
#     object_widgets = []
#     if state.objects:
#         for i, obj in enumerate(state.objects):
#             with vuetify.VExpansionPanel() as panel:
#                 with vuetify.VExpansionPanelHeader():
#                     with vuetify.VRow(align="center", justify="space_between", no_gutters=True):
#                         vuetify.VCol(obj.get("name"), cols="auto")
#                         # with vuetify.VCol(cols="auto", classes="flex-grow-0"):
#                         #     # ボタンをdivで囲み、クリックイベントに .stop を追加
#                         #     with html.Div(attributes={"@click.stop": ""}):
#                         #         print(f"DEBUG: Creating panel for object {i}: {obj.get('name', 'Unknown')}")
#                         #         vuetify.VBtn(
#                         #         icon=True,
#                         #         small=True,
#                         #         # ctrl.open_edit_dialog を呼び出すように修正
#                         #         click=lambda i=i: ctrl.open_edit_dialog(i)
#                         #     ).add_child(vuetify.VIcon("mdi-pencil", small=True))

#                         with vuetify.VCol(cols="auto", classes="flex-grow-0"):
#                             print(f"DEBUG: Creating panel for object {i}: {obj.get('name', 'Unknown')}")
#                             vuetify.VBtn(
#                                 icon=True,
#                                 small=True,
#                                 click=lambda i=i: state.set("edit_button_clicked", i)
#                             ).add_child(vuetify.VIcon("mdi-pencil", small=True))

#                 with vuetify.VExpansionPanelContent():
#                     vuetify.VSheet(html.Pre(json.dumps(obj, indent=2, ensure_ascii=False)), classes="pa-2", outlined=True, rounded=True)
#             object_widgets.append(panel)

#     layout.global_settings_container.children = settings_widgets
#     layout.object_settings_container.children = object_widgets
#     plotter.clear()
#     actors.clear()
#     original_colors.clear()
#     colors = {"water": "lightblue", "tank": "grey", "float": "yellow", "wavemaker": "green", "absorber": "orange", "box": "brown"}
#     if state.objects:
#         for obj in state.objects:
#             name = obj.get("name")
#             if "objfile" in obj and os.path.exists(obj["objfile"]):
#                 mesh = pv.read(obj["objfile"])
#                 color_key = next((k for k in colors if k in name.lower()), "default")
#                 color_value = colors.get(color_key, "white")
#                 actor = plotter.add_mesh(mesh, color=color_value, opacity=0.5 if "water" in name.lower() else 1.0)
#                 actors[name] = actor
#                 original_colors[actor] = color_value
#             elif "wave gauge" in obj.get("type", "").lower() and "position" in obj:
#                 pos = obj["position"]
#                 actor = plotter.add_mesh(pv.Line(pos[:3], pos[3:]), color="red", line_width=4)
#                 actors[name] = actor
#                 original_colors[actor] = "red"
#     plotter.reset_camera()
    
#     # 条件付きで view_update を呼び出す
#     if hasattr(ctrl, 'view_update') and callable(ctrl.view_update):
#         ctrl.view_update()




# `post_ui_setup`のロジックをこの関数に統合
@server.trigger("on_click")
def update_ui_and_plot(event=None):
    print("DEBUG: event in update_ui_and_plot:", event['remoteId'] if event and 'remoteId' in event else "No event")
    settings_widgets = []
    if state.settings:
        for key, value in state.settings.items():
            settings_widgets.append(vuetify.VCol(vuetify.VTextField(label=key, value=value, readonly=True, dense=True, outlined=True), cols=6, py_1=True))
    
    object_widgets = []
    if state.objects:
        for i, obj in enumerate(state.objects):
            with vuetify.VExpansionPanel() as panel:
                with vuetify.VExpansionPanelHeader():
                    with vuetify.VRow(align="center", justify="space_between", no_gutters=True):
                        vuetify.VCol(obj.get("name"), cols="auto")
                        with vuetify.VCol(cols="auto", classes="flex-grow-0"):
                            print(f"DEBUG: Creating panel for object {i}: {obj.get('name', 'Unknown')}")
                            vuetify.VBtn(
                                icon=True,
                                small=True,
                                click=lambda i=i: state.set("edit_button_clicked", i)
                            ).add_child(vuetify.VIcon("mdi-pencil", small=True))
                with vuetify.VExpansionPanelContent():
                    vuetify.VSheet(html.Pre(json.dumps(obj, indent=2, ensure_ascii=False)), classes="pa-2", outlined=True, rounded=True)
            object_widgets.append(panel)

    layout.global_settings_container.children = settings_widgets
    layout.object_settings_container.children = object_widgets
    plotter.clear()
    # actors.clear()
    original_colors.clear()
    # remote_id_map.clear() # IDマッピングをクリア
    colors = {"water": "lightblue", "tank": "grey", "float": "yellow", "wavemaker": "green", "absorber": "orange", "box": "brown"}
    # change color of     name = event['remoteId'] if event and 'remoteId' in event else None
    if event and 'remoteId' in event:
        name = remote_id_map.get(event['remoteId'])
        print(remote_id_map)
        print(f"DEBUG: Clicked object name: {name}")
        if name and name in actors:
            actor = actors[name]
            actor.prop.color = "red"
            last_picked_actor = actor
            state.selected_object_info = f"選択中: {name}"
            for i, obj in enumerate(state.objects):
                if obj.get("name") == name:
                    state.active_panel = i
                    break
        # nameと一致するcolorsの色も変更，他は変更しない
        for key, color in colors.items():
            if key in name.lower():
                colors[key] = "red"
                    
        

        # else:   
        #     # オブジェクトが見つからない場合は背景がクリックされたとみなす
        #     if last_picked_actor:
        #         last_picked_actor.prop.color = original_colors.get(last_picked_actor, "white")
        #         last_picked_actor = None
        #     state.selected_object_info = "背景がクリックされました"
        #     state.active_panel = None               


    if state.objects:
        for obj in state.objects:
            name = obj.get("name")
            if "objfile" in obj and os.path.exists(obj["objfile"]):
                mesh = pv.read(obj["objfile"])
                color_key = next((k for k in colors if k in name.lower()), "default")
                color_value = colors.get(color_key, "white")
                actor = plotter.add_mesh(mesh, color=color_value, opacity=0.5 if "water" in name.lower() else 1.0)
                actors[name] = actor
                original_colors[actor] = color_value
            elif "wave gauge" in obj.get("type", "").lower() and "position" in obj:
                pos = obj["position"]
                actor = plotter.add_mesh(pv.Line(pos[:3], pos[3:]), color="red", line_width=4)
                actors[name] = actor
                original_colors[actor] = "red"
    plotter.reset_camera()
    
    # 描画更新後にremoteIdマッピングを生成
    if view:
        print("DEBUG: post_ui_setup logic integrated")
        for name, actor in actors.items():
            remote_id = view.get_scene_object_id(actor)
            if remote_id:
                remote_id_map[remote_id] = name

    # 最後にまとめてビューを更新
    if hasattr(ctrl, 'view_update') and callable(ctrl.view_update):
        ctrl.view_update()



# UI構築後、アクターとremoteIdのマッピングを生成する
def post_ui_setup(view):
    """UI構築後に呼び出され、remoteIdマッピングを生成する"""
    print("DEBUG: post_ui_setup called")
    for name, actor in actors.items():
        remote_id = view.get_scene_object_id(actor)
        if remote_id:
            remote_id_map[remote_id] = name
            
    if hasattr(ctrl, 'view_update') and callable(ctrl.view_update):
        ctrl.view_update()

# --- UIレイアウト定義 ---
with SinglePageLayout(server) as layout:
    layout.title.set_text("BEM-MEL Interactive Pre-processor")
    
    with vuetify.VDialog(v_model=("edit_dialog", False), max_width="600px"):
        with vuetify.VCard():
            vuetify.VCardTitle("オブジェクト編集: {{ edited_object.name }}")
            with vuetify.VCardText():
                with vuetify.VContainer():
                    # v-forの代わりに、動的にフィールドを生成する方法に変更
                    with vuetify.VRow():
                        vuetify.VCol(
                            vuetify.VTextField(
                                label="名前",
                                v_model=("edited_object.name", ""),
                            ),
                            cols=6
                        )
                        vuetify.VCol(
                            vuetify.VTextField(
                                label="タイプ",
                                v_model=("edited_object.type", ""),
                            ),
                            cols=6
                        )
            # app.py の VCardActions ブロック

            with vuetify.VCardActions():
                vuetify.VSpacer()

                # キャンセルボタン: ブラウザ側で完結させるので attributes を使用
                vuetify.VBtn(
                    "キャンセル",
                    attributes={'@click': "edit_dialog = false"}
                )

                # 保存ボタン: Python側のロジックを呼び出すので click= を使用
                vuetify.VBtn(
                    "保存",
                    color="primary",
                    click=ctrl.save_changes
                )    

    with layout.content:
        with vuetify.VContainer(fluid=True, classes="pa-0 fill-height"):
            with vuetify.VRow(classes="fill-height"):
                with vuetify.VCol(cols=4, classes="pa-2"):
                    with vuetify.VCard(classes="fill-height", outlined=True):
                        vuetify.VCardTitle("⚙️ コントロールパネル")
                        vuetify.VDivider()
                        vuetify.VSelect(label="シミュレーションケースを選択", v_model=("selected_case", None), items=("case_options", ["Ruehl2016", "Tanizawa1996"]), dense=True, outlined=True, classes="ma-2")
                        vuetify.VCardSubtitle(text="選択オブジェクト")
                        vuetify.VChip("{{ selected_object_info }}", classes="ma-2", color="primary", outlined=True)
                        vuetify.VCardSubtitle(classes="mt-2", text="グローバル設定")
                        layout.global_settings_container = vuetify.VRow(classes="px-2")
                        vuetify.VDivider(classes="my-2")
                        vuetify.VCardSubtitle(text="インタラクティブ操作 (浮体)")
                        with vuetify.VCard(outlined=True, classes="pa-2 ma-2"):
                            vuetify.VSlider(label="X", v_model=("float_pos_x", 0), min=-5, max=5, step=0.1, thumb_label=True, dense=True)
                            vuetify.VSlider(label="Y", v_model=("float_pos_y", 0), min=-5, max=5, step=0.1, thumb_label=True, dense=True)
                            vuetify.VSlider(label="Z", v_model=("float_pos_z", 0), min=-5, max=5, step=0.1, thumb_label=True, dense=True)
                        vuetify.VDivider(classes="my-2")
                        vuetify.VCardSubtitle(text="オブジェクト設定")
                        layout.object_settings_container = vuetify.VExpansionPanels(v_model=("active_panel", None))

                # with vuetify.VCol(cols=8, classes="pa-0"):
                #     view = vtk.VtkLocalView(
                #         plotter.ren_win,
                #         picking_modes=["click"],
                #         attributes={'@click': "on_view_click($event)"},
                #     )
                #     ctrl.view_update = view.update
                with vuetify.VCol(cols=8, classes="pa-0"):
                    view = vtk.VtkLocalView(
                        plotter.ren_win,
                        # picking_modes=["click"],
                        # attributes={'picking_modes': "['click']"},
                        # click=ctrl.on_click, # <<< ここが公式の書き方
                        picking_modes=("picking_modes", ["click"]),
                        click=(on_view_click, "[$event]"),   
                    )
                    ctrl.view_update = view.update
                    
    # ctrl.on_server_ready.add(lambda **kwargs: update_ui_and_plot(view=view))

# --- サーバー起動 ---
if __name__ == "__main__":
    state.settings, state.objects = {}, []
    state.float_pos_x, state.float_pos_y, state.float_pos_z = 0, 0, 0
    state.selected_case = None
    state.selected_object_info = "オブジェクトをクリックして選択"
    state.active_panel = None
    state.edit_dialog = False
    state.edit_button_clicked = None
    state.edited_object_index = None
    state.edited_object = {}
    state.case_options = ["Ruehl2016", "Tanizawa1996"]  # これを追加
    server.start()