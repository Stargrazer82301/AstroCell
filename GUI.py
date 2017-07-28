# Import smorgasbord
import tkinter

class MyFirstGUI:

    def __init__(self, master):
        self.master = master
        master.title("AstroCell")

        self.label = tkinter.Label(master, text="\n Good lord, a GUI!")
        self.label.pack()

        self.greet_button = tkinter.Button(master, text="Hello", font=('helvetica',15,"bold"), command=self.greet)
        self.greet_button.pack()

        self.close_button = tkinter.Button(master, text="Close", command=master.destroy)
        self.close_button.pack()

        self.centerWindow()

    def centerWindow(self):

        w = 500
        h = 250

        sw = self.master.winfo_screenwidth()
        sh = self.master.winfo_screenheight()

        x = (sw - w)/2
        y = (sh - h)/2
        self.master.geometry('%dx%d+%d+%d' % (w, h, x, y))


    def greet(self):
        print("Greetings!")

# Run GUI
root = tkinter.Tk()
my_gui = MyFirstGUI(root)
root.mainloop()
