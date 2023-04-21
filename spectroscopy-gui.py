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

class FileItem():                       # This is the gui object
    """docstring for FileItem
    """
    def __init__(self, file_name):
        self.file_name = file_name
        self.select = False
        self.dirty = False
        self.spectrum = spec.Spectrum()
        self.plots = []

class App(tk.Tk):
    def __init__(self):
        self.file_list = []             # List of loaded files
        self.plot_list = []             # List of open plots
        self.recent_list = []           # List of recently closed files

        tk.Tk.__init__(self)            # Startup the gui

        self.title("Spectroscopy")
        self.geometry("400x150")

        # Set up the menus for the Application
        menu_bar = tk.Menu(self)
        file_menu = tk.Menu(menu_bar, tearoff=0)
        file_menu.add_command(label="Open", command=self.open_file)
        file_menu.add_command(label="Open Recent", command=save_file)
        file_menu.add_command(label="Import", command=save_file)
        file_menu.add_separator()
        file_menu.add_command(label="Save", command=save_file)
        file_menu.add_command(label="Save As", command=save_file)
        file_menu.add_command(label="Export", command=save_file)
        file_menu.add_separator()
        file_menu.add_command(label="Close", command=save_file)
        file_menu.add_command(label="Close All", command=save_file)
        file_menu.add_command(label="Exit", command=self.quit )
        menu_bar.add_cascade(label="File", menu=file_menu)

        edit_menu = tk.Menu(menu_bar, tearoff=0)
        edit_menu.add_command(label="Undo", command=save_file)
        edit_menu.add_command(label="Redo", command=save_file)
        edit_menu.add_separator()
        edit_menu.add_command(label="Edit history", command=edit_data)
        edit_menu.add_command(label="Edit acquisition", command=save_file)
        edit_menu.add_command(label="Edit metadata", command=save_file)
        edit_menu.add_separator()
        edit_menu.add_command(label="Convert x", command=save_file)
        edit_menu.add_command(label="Convert y", command=save_file)
        edit_menu.add_separator()
        edit_menu.add_command(label="Shift x", command=save_file)
        edit_menu.add_command(label="Shift y", command=save_file)
        edit_menu.add_separator()
        edit_menu.add_command(label="Resample", command=save_file)
        menu_bar.add_cascade(label="Edit", menu=edit_menu)

        plot_menu = tk.Menu(menu_bar, tearoff=0)
        plot_menu.add_command(label="New plot", command=save_file)
        plot_menu.add_command(label="Close plot", command=save_file)
        plot_menu.add_command(label="Close all plots", command=save_file )
        menu_bar.add_cascade(label="Plot", menu=plot_menu)

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

        help_menu = tk.Menu(menu_bar, tearoff=0)
        help_menu.add_command(label="About", command=save_file)
        help_menu.add_command(label="Help", command=save_file)
        menu_bar.add_cascade(label="Help", menu=help_menu)

        self.config(menu=menu_bar)

    def quit(self):
        self.close_all_plots
        self.close_all_files
        self.final_splash
        sys.exit(0)

    def open_file(self):
        file_name = tk.filedialog.askopenfilename()
        new_file = File_item(file_name)
        self.filelist.append(new_file)

    def close_all_plots(self):
        for plot_item in self.plot_list:
            pass

def edit_data():
    dialog = tk.Toplevel(self)
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
    app = App()
    app.mainloop()
