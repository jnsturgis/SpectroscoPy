# Chat GPT example to use Matplotlib

import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np

# Create a Tkinter window
root = tk.Tk()
root.title("My Plot")

# Create a frame to hold the plot
frame = tk.Frame(root)
frame.pack()

# Chat GPT for menubars

menu_bar = tk.Menu(root)
file_menu = tk.Menu(menu_bar, tearoff=0)
file_menu.add_command(label="Open", command=open_file)
file_menu.add_command(label="Save", command=save_file)
menu_bar.add_cascade(label="File", menu=file_menu)
root.config(menu=menu_bar)

# and for editing data structure dialogs

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

def save_data(new_value):
    # Code to save the new value to the data structure
    pass

# Create a Matplotlib figure and add an Axes object
fig = Figure(figsize=(5, 4), dpi=100)
ax = fig.add_subplot(111)

# Generate some data and plot it on the Axes
x = np.linspace(0, 10, 100)
y = np.sin(x)
ax.plot(x, y)

# Create a Matplotlib canvas and add it to the frame
canvas = FigureCanvasTkAgg(fig, master=frame)
canvas.draw()
canvas.get_tk_widget().pack()

# Run the GUI
root.mainloop()

# From Stackoverflow 
class App(tk.Tk):

    def __init__(self):
        tk.Tk.__init__(self)
        menubar = tk.Menu(self)
        fileMenu = tk.Menu(menubar, tearoff=False)
        menubar.add_cascade(label="File", underline=0, menu=fileMenu)
        fileMenu.add_command(label="Exit", underline=1,
                             command=quit, accelerator="Ctrl+Q")
        self.config(menu=menubar)

        self.bind_all("<Control-q>", self.quit)

    def quit(self, event):
        print("quitting...")
        sys.exit(0)

if __name__ == "__main__":
    app = App()
    app.mainloop()
