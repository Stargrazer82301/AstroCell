# Import smorgasbord
import pdb
import os
import tkinter
import tkinter.filedialog





class OptionsGUI:
    """ Class that operates the tkinter GUI for AstroCell """

    def __init__(self, master):
        """ Initialise the GUI """

        # Incorporate master tk instance into object, and title window 'AstroCell'
        self.master = master
        self.master.title("AstroCell")

        # Add widgets for user to select target directory
        self.SelectDir()

        # Add widgets for user to select number of cell colour types to be included
        self.SelectColours()

        # Add widgets for user to select MC factor
        self.SelectMC()

        # Button to cancel AstroCell
        self.cancel_button = tkinter.Button(self.master, text='Cancel', font=('helvetica',12), command=self.master.destroy)
        self.cancel_button.grid(row=4, column=1, pady=5)

        # Button to run AstroCell with the chosen settings
        self.run_button = tkinter.Button(self.master, text='Run', font=('helvetica',12,"bold"), command=self.PrintDir)
        self.run_button.grid(row=4, column=0, pady=5)

        # Centre window, and have widgets react dynamically to window resizing
        self.CenterWindow()



    def SelectDir(self):
        """ Method that adds widgets for user to select target directory """

        # Entry box for choosing target directory
        self.dir_path = os.getcwd()
        self.dir_entry = tkinter.Entry(self.master, font=('helvetica',12,"bold"), width=50)
        self.dir_entry.insert(0, os.getcwd())
        self.dir_entry.grid(row=0, column=0, columnspan=2, sticky='W', padx=10, pady=(10,2))

        # Button for choosing target directory
        self.dir_button = tkinter.Button(self.master, text='Browse', font=('helvetica',12,"bold"), command=self.GetDir)
        self.dir_button.grid(row=1, column=1, pady=(2,10))



    def SelectColours(self):
        """ Method that adds widgets for user to state how many cell types are to be classified """

        # Create main label for entry
        self.cell_colours_label = tkinter.Label(self.master, text='Number of cell types to classify:', font=('helvetica',12,"bold"))
        self.cell_colours_label.grid(row=2, column=0, sticky='W')

        # Create entry box
        self.cell_colours_entry = tkinter.Entry(self.master, font=('helvetica',12,"bold"), width=10)
        self.cell_colours_entry.insert(0, '2')
        self.cell_colours_entry.grid(row=2, column=1, sticky='W')



    def SelectMC(self):
        """ Method that adds widgets for user to control Monte-Carlo factor governing cmputational intensity of deblending """

        # Create main label for factor entry
        self.mc_factor_label = tkinter.Label(self.master, text='Computational intensity:', font=('helvetica',12,"bold"))
        self.mc_factor_label.grid(row=3, column=0, sticky='W')

        # Create entry box
        self.mc_factor_entry = tkinter.Entry(self.master, font=('helvetica',12,"bold"), width=10)
        self.mc_factor_entry.insert(0, '1')
        self.mc_factor_entry.grid(row=3, column=1, sticky='W')



    def GetDir(self):
        """ Method to bring up file browser dialogue box to allow user to select target directory """

        # Bring up dialogue box
        self.dir_path = tkinter.filedialog.askdirectory(initialdir=os.getcwd(), title='Select folder to be processed')

        # If file browser dialogue was cancelled, reset target directory to current directory
        if self.dir_path == ():
            self.dir_path = os.getcwd()

        # Update text in directory entry box
        self.dir_entry.delete(0, 'end')
        self.dir_entry.insert(0, self.dir_path)



    def CenterWindow(self):
        """ Method that places the window in the centre of the screen """

        # Determine current window dimensions
        self.master.update()
        window_dims = tuple(int(_) for _ in self.master.geometry().split('+')[0].split('x'))
        window_w = window_dims[0]#650
        window_h = window_dims[1]#350

        # Find the dimensions of the screen
        screen_w = self.master.winfo_screenwidth()
        screen_h = self.master.winfo_screenheight()

        # Calculate correct positioning for window to be centred
        window_x = (screen_w - window_w)/2
        window_y = (screen_h - window_h)/2

        # Define window positioning geometry using determined values
        self.master.geometry('%dx%d+%d+%d' % (window_dims+(window_x,window_y)))

        # Weight grid so that widgets relocate to fill window as it resizes
        for r in range(0, self.master.grid_size()[0]):
            self.master.rowconfigure(r, weight=1)
        for c in range(0, self.master.grid_size()[1]):
            self.master.columnconfigure(c, weight=1)



    def PrintDir(self):
        """ Dummy method to print current target directory """

        # Print target directory
        print(self.dir_path)





# Commence main task of program
if __name__ == "__main__":

    # Initialise tkinter instance, and run GUI
    root = tkinter.Tk()
    root.deiconify()
    gui = OptionsGUI(root)
    root.mainloop()
