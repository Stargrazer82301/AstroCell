# Import smorgasbord
import pdb
import os
import sys
import tkinter
import tkinter.filedialog
import AstroCell.Main

# Various imports to force PyInstaller to work
import six
import packaging
import packaging.version
import packaging.specifiers
makework_list = [six, packaging, packaging.version, packaging.specifiers]
makework_list = [ makework.__class__ for makework in makework_list ]





class OptionsGUI:
    """ Class that operates the tkinter GUI for AstroCell """



    def __init__(self, master):
        """ Initialise the GUI """

        # Incorporate master tk instance into object, and title window 'AstroCell'
        self.master = master
        self.master.title("AstroCell")

        # By default, AstroCell will not run when window closed (ie, uncless specifically instructed to run)
        self.run_bool = False

        # Add widgets for user to select target directory
        self.SelectDir()

        # Add widgets for user to select number of cell colour types to be included
        self.SelectColours()

        # Add widgets for selecting paralle operation
        self.SelectParallel()

        # Add widgets for user to select MC factor
        self.SelectMC()

        # Add widgets to run AstroCell
        self.SelectRun()

        # Add widgets to cancel AstroCell
        self.SelectCancel()

        # Centre window, and have widgets react dynamically to window resizing
        self.CenterWindow()



    def SelectDir(self):
        """ Method that adds widgets for user to select target directory """

        # Create main label for entry
        self.dir_label = tkinter.Label(self.master, text='Select target directory:', font=('helvetica',12,"bold"))
        self.dir_label.grid(row=0, column=0, sticky='w', padx=10, pady=(2,0))

        # Button for choosing target directory
        self.dir_button = tkinter.Button(self.master, text='Browse', font=('helvetica',12,"bold"), command=self.GetDir)
        self.dir_button.grid(row=0, column=1, sticky='e', padx=10, pady=(20,0))

        # Entry box for choosing target directory
        self.dir_path = os.getcwd()
        self.dir_entry = tkinter.Entry(self.master, font=('helvetica',12,"bold"), width=75)
        self.dir_entry.insert(0, os.getcwd())
        self.dir_entry.grid(row=1, column=0, columnspan=2, sticky='e', padx=10, pady=(2,2))

        # Add explanatory text for entry
        self.dir_text = tkinter.Label(self.master, font=('helvetica',11),
                                      text='All image files found inside the target diretory will be processed by AstroCell.')
        self.dir_text.grid(row=2, column=0, columnspan=2, sticky='w', padx=10, pady=(2,20))



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



    def SelectColours(self):
        """ Method that adds widgets for user to state how many cell types are to be classified """

        # Create main label for entry
        self.cell_colours_label = tkinter.Label(self.master, text='Number of cell types to classify:', font=('helvetica',12,"bold"))
        self.cell_colours_label.grid(row=3, column=0, sticky='w', padx=10, pady=(20,2))

        # Create entry box
        self.cell_colours_entry = tkinter.Entry(self.master, font=('helvetica',12,"bold"), width=10)
        self.cell_colours_entry.insert(0, '2')
        self.cell_colours_entry.grid(row=3, column=1, sticky='e', padx=10, pady=(20,0))

        # Add explanatory text for entry
        self.dir_text = tkinter.Label(self.master, font=('helvetica',11), anchor='w', justify='left', wraplength=650,
                                      text='Number of cell colour types AstroCell should try to classify. Default is 2 (ie, blue-green, blue-brown).')
        self.dir_text.grid(row=4, column=0, columnspan=2, sticky='w', padx=10, pady=(2,20))



    def SelectParallel(self):
        """ Method that adds widgets for user to choose whether AstroCell should operate in parallel """

        # Create main label for entry
        self.parallel_label = tkinter.Label(self.master, text='Operate in parallel:', font=('helvetica',12,"bold"))
        self.parallel_label.grid(row=5, column=0, sticky='w', padx=10, pady=(20,2))

        parallel_menu_list = ['Yes','No']
        self.parallel_menu_var = tkinter.StringVar()
        self.parallel_menu_var.set('Yes')
        self.parallel_dropdown = tkinter.OptionMenu(self.master, self.parallel_menu_var, *parallel_menu_list).grid(row=5, column=1, sticky='e', padx=10, pady=(20,2))

        # Add explanatory text for menu
        self.dir_text = tkinter.Label(self.master, font=('helvetica',11), anchor='w', justify='left', wraplength=650,
                                      text='Select whether AstroCell should operate in parallel. If use, AstroCell will use all but one of this computer\'s CPU cores. Default is yes.')
        self.dir_text.grid(row=6, column=0, columnspan=2, sticky='w', padx=10, pady=(2,20))



    def SelectMC(self):
        """ Method that adds widgets for user to control Monte-Carlo factor governing cmputational intensity of deblending """

        # Create main label for factor entry
        self.mc_factor_label = tkinter.Label(self.master, text='Computational intensity:', font=('helvetica',12,"bold"))
        self.mc_factor_label.grid(row=7, column=0, sticky='w', padx=10, pady=(20,2))

        # Create entry box
        self.mc_factor_entry = tkinter.Entry(self.master, font=('helvetica',12,"bold"), width=10)
        self.mc_factor_entry.insert(0, '1')
        self.mc_factor_entry.grid(row=7, column=1, sticky='e', padx=10, pady=(20,2))

        # Add explanatory text for entry
        self.dir_text = tkinter.Label(self.master, font=('helvetica',11), anchor='w', justify='left', wraplength=650,
                                      text='How much processing to perform when de-blending cells (the most computationally-intensive stage of AstroCell). Defult of 1 typically provides good balance of effectiveness and speed; setting to 2 will run de-blending for twice as long, 0.5 will run for half as long, etc.')
        self.dir_text.grid(row=8, column=0, columnspan=2, sticky='w', padx=10, pady=(2,20))



    def SelectRun(self):
        """ Method that adds widgets to start AstroCell """

        # Button to run AstroCell with the chosen settings
        self.run_button = tkinter.Button(self.master, text='Run', font=('helvetica',12,"bold"), command=self.AstroCellRun)
        self.run_button.grid(row=11, column=0, padx=10, pady=(20,20))



    def SelectCancel(self):
        """ Method that adds widgets to cancel AstroCell """

        # Button to cancel AstroCell
        self.cancel_button = tkinter.Button(self.master, text='Cancel', font=('helvetica',12), command=self.master.destroy)
        self.cancel_button.grid(row=11, column=1, padx=10, pady=(20,20))



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



    def AstroCellRun(self):
        """ Method to record that AstroCell is to run, then exit GUI """

        # Extract cell colour entry
        self.cell_colours = int(self.cell_colours_entry.get())
        self.parallel = self.parallel_menu_var.get()

        # Extract paralell option
        if self.parallel == 'Yes':
            self.parallel = True
        elif self.parallel == 'No':
            self.parallel = False

        # Extract Monte-Carlo factor
        self.mc_factor = float(self.mc_factor_entry.get())

        # Change run boolean to be True (from defaut of False), then destroy the whole window
        self.run_bool = True
        self.master.destroy()





# Commence main task of program
if __name__ == "__main__":

    # Initialise tkinter instance, and run GUI
    root = tkinter.Tk()
    root.deiconify()
    gui = OptionsGUI(root)
    root.mainloop()

    # Evalulate if AstroCell was instructed to run; if so, run using the parameters provided by user
    if gui.run_bool:
        AstroCell.Main.Run(in_dir=gui.dir_path,
                           cell_colours=gui.cell_colours,
                           substructure_flag=False,
                           parallel=gui.parallel,
                           mc_factor=gui.mc_factor,
                           dill_dir=False)
