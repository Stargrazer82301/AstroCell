# Import smorgasbord
import pdb
import os
import tkinter
import tkinter.filedialog


class OptionsGUI:
    """ Class that operates the tkinter GUI for AstroCell """

    def __init__(self, master):
        """ Initialise the GUI """

        self.master = master
        master.title("AstroCell")

        self.label = tkinter.Label(master, text="\n Good lord, a GUI!")
        self.label.grid(columnspan=1, rowspan=10, sticky=tkinter.W)

        self.greet_button = tkinter.Button(master, text="Hello", font=('helvetica',15,"bold"), command=self.greet)
        self.greet_button.grid(row=2, column=1)

        self.close_button = tkinter.Button(master, text="Close", font=('helvetica',15,"bold"), command=master.destroy)
        self.close_button.grid(row=4, column=2)

        self.CenterWindow()

    def CenterWindow(self):
        """ Method that places the window in the centre of the screen """

        # State the desired width and height of the window
        window_w = 500
        window_h = 250

        # Find the dimensions of the screen
        screen_w = self.master.winfo_screenwidth()
        screen_h = self.master.winfo_screenheight()

        # Calculate correct positioning for window to be centred
        window_x = (screen_w - window_w)/2
        window_y = (screen_h - window_h)/2

        # Define window positioning geometry using determined values
        pdb.set_trace()
        self.master.geometry('%dx%d+%d+%d' % (window_w, window_h, window_x, window_y))


    def greet(self):
        print("Greetings!")



# Run an initial tk instance, to request path to directory the contents of which are to be processed by AstroCell
root = tkinter.Tk()
root.withdraw()
root.dir_name = tkinter.filedialog.askdirectory(initialdir=os.getcwd(), title="Select folder to be processed")
print(root.dir_name)

# Run options GUI, to provide setup for AstroCell
root.deiconify()
gui = OptionsGUI(root)
root.mainloop()
