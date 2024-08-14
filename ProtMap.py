import h5py
import umap
import dearpygui.dearpygui as dpg
import numpy as np

directory = "embeddings"

embeddings_list = []
protein_ids = []
x_coords = []
y_coords = []
minindex = 0
flagPopup = False

import os
from pymol import cmd


def close_popup():
    global flagPopup
    dpg.delete_item("popup")
    flagPopup = False


def generate_png():
    dir = "pngs"
    for ProtId in protein_ids:
        cmd.reinitialize()
        pdb_file = f"pdbs/{ProtId}.pdb"
        if os.path.isfile(pdb_file):
            cmd.load(pdb_file)

            cmd.remove("resn HOH")

            cmd.show("ribbon")
            cmd.bg_color("white")
            cmd.spectrum("count", "rainbow")

            cmd.set('ray_opaque_background', 0)
            cmd.zoom(complete=1)

            # Save an image
            output_file = f"{dir}/{ProtId}.png"
            cmd.png(output_file, width=300, height=300, dpi=600, ray=1)


def show_popup():
    nearPointIndex = 0
    global flagPopup
    minDist = 9999999999999999
    global minindex
    if not flagPopup:
        with dpg.popup(plot_id, tag="popup", mousebutton=dpg.mvMouseButton_Right, modal=True):
            dpg.add_text("", tag="popup_text")
            dpg.add_button(label="Close Popup",
                           callback=lambda: [dpg.configure_item("popup", show=False), close_popup()])
        plot_pos = dpg.get_plot_mouse_pos()

        near_point = False
        minindex = 0
        for x, y in zip(x_coords, y_coords):
            valu = pow(plot_pos[0] - x, 2) + pow(plot_pos[1] - y, 2)
            if valu < minDist:
                minDist = valu
                minindex = nearPointIndex

            nearPointIndex += 1
        if minDist < 0.02:
            near_point = True

        if near_point:
            protein_id = protein_ids[minindex]
            dpg.set_value("popup_text", f"Protein name: {protein_id}")

            with dpg.drawlist(width=200, height=200, parent="popup"):
                dpg.draw_image(protein_id, (0, 0), (200, 200), tag="imageDraw")
        else:
            dpg.set_value("popup_text", "Right-click near a data point")

        dpg.configure_item("popup", show=near_point)
        flagPopup = True


for file in os.listdir(directory):
    file2 = os.path.join(directory, file)
    with h5py.File(file2, 'r') as f:
        for protein_id in f.keys():
            embedding = np.array(f[protein_id])
            embedding = np.mean(embedding, axis=0)
            embeddings_list.append(embedding)
            protein_ids.append(protein_id)

print("done sampling")
# generate_png() #use only for data gathering
reducer = umap.UMAP()
embedding_umap = reducer.fit_transform(embeddings_list)
x_coords = embedding_umap[:, 0].copy()
y_coords = embedding_umap[:, 1].copy()

print("done Umap")
print(x_coords)
alphaXCoords = []
alphaYCoords = []
betaXCoords = []
betaYCoords = []
alphaBetaXCoords = []
alphaBetaYCoords = []
membraneXCoords = []
membraneYCoords = []
naXCoords = []
naYCoords = []
classNames = []
for i in range(len(x_coords)):
    file = open("classes/" + protein_ids[i])
    m = file.read()
    classNames.append(m)
    if m == "alpha":
        alphaXCoords.append(x_coords[i])
        alphaYCoords.append(y_coords[i])
    elif m == "beta":
        betaXCoords.append(x_coords[i])
        betaYCoords.append(y_coords[i])
    elif m == "alpha-beta":
        alphaBetaXCoords.append(x_coords[i])
        alphaBetaYCoords.append(y_coords[i])
    elif m == "membrane":
        membraneXCoords.append(x_coords[i])
        membraneYCoords.append(y_coords[i])
    else:
        naXCoords.append(x_coords[i])
        naYCoords.append(y_coords[i])

for i in range(len(x_coords)):
    print(x_coords[i])
    print(y_coords[i])
    print(protein_ids[i])

dpg.create_context()
for file in os.listdir("pngs"):
    file2 = os.path.join("pngs", file)
    width, height, channels, data = dpg.load_image(file2)
    temp1 = file2.find('.')

    filename1 = file2[temp1 - 4:temp1]
    print(filename1)
    with dpg.texture_registry():
        dpg.add_static_texture(width, height, data, tag=filename1)

with dpg.handler_registry():
    dpg.add_mouse_click_handler(callback=show_popup)

with dpg.window(label="Tutorial", width=1200, height=1200, tag="mainWindow"):
    dpg.add_text("Right click near a point to see name and picture")
    with dpg.plot(label="Line Series", height=-1, width=-1) as plot_id:
        dpg.add_plot_legend()
        dpg.add_plot_axis(dpg.mvXAxis, label="x")
        dpg.add_plot_axis(dpg.mvYAxis, label="y", tag="yaxis")

        dpg.add_scatter_series(alphaXCoords, alphaYCoords, label="alpha", parent="yaxis", tag="alpha")
        dpg.add_scatter_series(betaXCoords, betaYCoords, label="beta", parent="yaxis", tag="beta")
        dpg.add_scatter_series(alphaBetaXCoords, alphaBetaYCoords, label="alpha-beta", parent="yaxis", tag="alpha-beta")
        dpg.add_scatter_series(membraneXCoords, membraneYCoords, label="membrane", parent="yaxis", tag="membrane")
        dpg.add_scatter_series(naXCoords, naYCoords, label="N/A", parent="yaxis", tag="na")

        dpg.add_button(label="Delete Series 1", parent=dpg.last_item(), callback=lambda: dpg.delete_item("series_1"))

dpg.create_viewport(title='Custom Title', width=1200, height=1200)
dpg.setup_dearpygui()
dpg.show_viewport()
dpg.start_dearpygui()
dpg.destroy_context()
