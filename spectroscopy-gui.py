""" An application for interactive treatment and analysis of spectra

This graphicaly interfaced application is conceived to use the spectroscopy
library, matplotlib and gui library to provide an easy to use application
for the treatment and analysis of spectra from different sources.

The gui opens a control window that will show the list of loaded data and can
be used to manipulate the data items via menus.

loading a spectrum adds a name to the list, and depending on the parameters
might show a plot of the spectrum in a new window or add it to the current
plot with or without a rescale of the axes.

Plots are associated with a dataframe derived from the spectra but not
necessarilly identical because of unit changes offsets etc.
"""

import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np

import spectroscopy as spec
import spectroscopy.gui as gui

def main():
# Chat GPT example to use Matplotlib

# Create a Tkinter window
    root = tk.Tk()
    root.title("Spectroscopy")
    root.geometry("400x150")

    # Chat GPT for menubars
    menu_bar = tk.Menu(root)
    file_menu = tk.Menu(menu_bar, tearoff=0)
    file_menu.add_command(label="Open", command=open_file)
    file_menu.add_command(label="Open Recent", command=save_file)
    file_menu.add_command(label="Import", command=save_file)
    file_menu.add_command(label="Save", command=save_file)
    file_menu.add_command(label="Save As", command=save_file)
    file_menu.add_command(label="Export", command=save_file)
    file_menu.add_command(label="Plot", command=save_file)
    file_menu.add_command(label="Add to Plot", command=save_file)
    file_menu.add_command(label="Close", command=save_file)
    file_menu.add_command(label="Close All", command=save_file)
    file_menu.add_command(label="Exit", command=root.destroy )
    menu_bar.add_cascade(label="File", menu=file_menu)

    edit_menu = tk.Menu(menu_bar, tearoff=0)
    edit_menu.add_command(label="Undo", command=save_file)
    edit_menu.add_command(label="Redo", command=save_file)
    edit_menu.add_command(label="Edit history", command=edit_data)
    edit_menu.add_command(label="Edit acquisition", command=save_file)
    edit_menu.add_command(label="Edit metadata", command=save_file)
    edit_menu.add_command(label="Convert x", command=save_file)
    edit_menu.add_command(label="Convert y", command=save_file)
    edit_menu.add_command(label="Shift x", command=save_file)
    edit_menu.add_command(label="Shift y", command=save_file)
    edit_menu.add_command(label="Resample", command=save_file)
    menu_bar.add_cascade(label="Edit", menu=edit_menu)

    treat_menu = tk.Menu(menu_bar, tearoff=0)
    treat_menu.add_command(label="Smooth", command=save_file)
    treat_menu.add_command(label="Differentiate", command=save_file)
    treat_menu.add_command(label="Integrate", command=save_file)
    treat_menu.add_command(label="Baseline subtract", command=save_file)
    treat_menu.add_command(label="Spectrum calculator", command=save_file)
    treat_menu.add_command(label="Sharpen peaks", command=save_file)
    treat_menu.add_command(label="Scattering correction", command=save_file)
    treat_menu.add_command(label="ATR correction", command=save_file)
    menu_bar.add_cascade(label="Treatment", menu=treat_menu)

    analysis_menu = tk.Menu(menu_bar, tearoff=0)
    analysis_menu.add_command(label="Peak fitting", command=save_file)
    analysis_menu.add_command(label="Identification", command=save_file)
    analysis_menu.add_command(label="Spectrum fitting", command=save_file)
    menu_bar.add_cascade(label="Analysis", menu=analysis_menu)

    root.config(menu=menu_bar)

# Run the GUI
    root.mainloop()


def open_file():

    file_name = tk.filedialog.askopenfilename()
    # Code to open the file and display it in the GUI

def edit_data():
    dialog = tk.Toplevel(root)
    dialog.title("Edit Data")
    dialog.geometry("300x200")
    label = tk.Label(dialog, text="Enter new value:")
    label.pack()

    entry = tk.Entry(dialog)
    entry.pack()

    button = tk.Button(dialog, text="OK", command=lambda: save_data(entry.get()))
    button.pack()

def save_file():
    # Code to save the new value to the data structure
    pass

if __name__ == '__main__':
    main()
