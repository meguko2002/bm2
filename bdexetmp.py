import PySimpleGUI as sg


# sg.theme('Light Blue 2')
# # ファイル選択
filename = sg.popup_get_file('ファイル選択')

# 入力設定ウィンドウ
layout = [[sg.Text(''.join(['file :', filename]))],
          [sg.Text('データ元を選んでください')],
          [sg.Button('プリンタ'), sg.Button('ラインカメラ'), sg.Button('反射センサ')],
          ]

window1 = sg.Window('Window 1', layout)

window1.close()
