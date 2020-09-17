from tkinter import *
from tkinter import ttk
from PIL import ImageTk, Image
from tkinter import filedialog
import NGS_app
from tkinter import messagebox
#import pandas._libs.tslibs.base
import sklearn.utils._cython_blas
try:
    from ctypes import windll
    windll.shcore.SetProcessDpiAwareness(1)
except:
    pass


root = Tk()
root.title("Amplicon sequences data visualization")
root.geometry("1550x700")
root.iconbitmap("dna.ico")

style = ttk.Style(root)
style.theme_use("clam")

style.map("CustomButton.TButton", foreground=[("pressed", "red")])


def button_hover(e):
    upload_tax["bg"] = "white"
    status_label.config(text="This table must contain the OTUs taxonomy")


def button_hover_1(e):
    upload_otu["bg"] = "white"
    status_label.config(text="This table must contain the subsumpled OTUs file")


def button_hover_leave(e):
    upload_tax["bg"] = "SystemButtonFace"
    status_label.config(text="")


def button_hover_leave_1(e):
    upload_otu["bg"] = "SystemButtonFace"
    status_label.config(text="")


def button_hover_concat(e):
    pars_button["bg"] = "white"
    status_label.config(text="Both taxonomy and OTUs tables will be concatenated and organised")


def button_hover_concat_leave(e):
    pars_button["bg"] = "SystemButtonFace"
    status_label.config(text="")


def button_hover_calc(e):
    calc_button["bg"] = "white"
    status_label.config(text="Calculations based on the chosen taxonomic level and percentage")


def button_hover_calc_leave(e):
    calc_button["bg"] = "SystemButtonFace"
    status_label.config(text="")


def button_heat_map_enter(e):
    heatmap_button["bg"] = "white"
    status_label.config(text="Heatmap is expected")


def button_heat_map_leave(e):
    heatmap_button["bg"] = "SystemButtonFace"
    status_label.config(text="")


def button_barplot_enter(e):
    bar_button["bg"] = "white"
    status_label.config(text="Barplot is expected")


def button_barplot_leave(e):
    bar_button["bg"] = "SystemButtonFace"
    status_label.config(text="")


def pca_enter(e):
    pca_button["bg"] = "white"
    status_label.config(text="Principal component analyses graph")


def pca_leave(e):
    pca_button["bg"] = "SystemButtonFace"
    status_label.config(text="")


def create_bar():
    frame2.destroy()
    frame3 = Frame(root, bd=1, background="white", relief="sunken", height=490, width=980)
    frame3.grid(row=0, column=2, pady=10, padx=10, rowspan=3)
    image_bar = NGS_app.barplot(dnew)
    image_bar_var = ImageTk.PhotoImage(Image.open('Output_1.jpg'))
    im_label = Label(frame3, image=image_bar_var)
    im_label.image = image_bar_var
    im_label.pack()


def upload_tax_f():
    root.filename_tax = filedialog.askopenfilename(initialdir="/",
                                               title="Open a text file",
                                               filetypes=(("TXT File", "*.txt"),))
    global O
    O = NGS_app.pars_tax(root.filename_tax)
    # print(O)


def upload_otu_f():
    root.filename_otu = filedialog.askopenfilename(initialdir="/",
                                               title="Open a text file",
                                               filetypes=(("TXT File", "*.txt"),))
    global T
    T = NGS_app.pars_otu(root.filename_otu)
    # print(T)


def concat_tables():
    global data
    data = NGS_app.concat_tables(T, O)
    data.to_csv("amplicon_sequencing_file.csv")
    popup()


# create pop up function
def popup():
    messagebox.showinfo("message box", "The file has been downloaded")


def create_heat():
    image_heatmap = NGS_app.heatmap_plot(dnew)
    image_heatmap_var = ImageTk.PhotoImage(Image.open('Output.jpg'))
    #frame2.destroy()
    im_label = Label(frame2, image=image_heatmap_var)
    im_label.image = image_heatmap_var
    im_label.pack()




def open_window():
    window = Toplevel()
    window.title("Principal components analyses")
    window.geometry("1130x600")
    window.iconbitmap("dna.ico")
    options = ["kingdom", "phylum", "class", "order", "family", "genus"]
    combo_taxa = ttk.Combobox(window, value=options, state="readonly")
    combo_taxa.current(5)
    combo_taxa.grid(row=0, column=0)
    options_perc = [0, 5, 10, 15]
    combo_perc = ttk.Combobox(window, value=options_perc, state="readonly")
    combo_perc.current(0)
    combo_perc.grid(row=0, column=1)
    frame_b_pca = Frame(window, bd=1, bg="white", relief="sunken", height=505, width=1005)
    frame_b_pca.grid(row=1, column=0, columnspan=2)


    def show_pca():
        global data
        global dnew
        PCA_pic = NGS_app.PCA_f(combo_taxa.get(), int(combo_perc.get()), data, dnew)
        image_PCA_var = ImageTk.PhotoImage(Image.open('Output_2.jpg'))
        #frame_b_pca.destroy()
        im_label = Label(frame_b_pca, image=image_PCA_var)
        im_label.image = image_PCA_var
        im_label.pack()
    show_pca_button = Button(window, text="Show", height=2, width=10,
                             activebackground="gray", command=show_pca)
    show_pca_button.grid(row=0, column=2)


expl_label = Label(root, text="This application is built on the tables generated by the mothur software,"
                              " in order to proceed to the visualization the two files .subsample.shared"
                              " and .cons.taxonomy are needed",
                   bg="white", font=("Helvtica", 10),
                   relief="raised", state="disabled")
expl_label.grid(row=4, column=0, columnspan=3, pady=5, padx=5)

frame = Frame(root, bd=1, bg="white", relief="sunken", height=385, width=150)
frame.grid(row=0, column=1)

frame_Label = Label(frame, text="Choose the desired taxonomic level", bg="white")
frame_Label.pack(padx=20, pady=20)

frame1 = Frame(root, bd=1, bg="white", relief="sunken", height=170, width=240)
frame1.grid(row=1, column=1)

frame1_label = Label(frame1, text="This field contains explanations about the results:  \nThe samples names\nand the umber of distinct taxa\nabove the treshhold",
                     bg="white", fg="gray")
frame1_label.pack(padx=20, pady=20)

frame2 = Frame(root, bd=1, bg="white", relief="sunken", height=480, width=645)
frame2.grid(row=0, column=2, rowspan=3, padx=20, pady=20)

# list of taxa combo box
options = ["kingdom", "phylum", "class", "order", "family", "genus"]
combo_taxa = ttk.Combobox(frame, value=options, state="readonly")
combo_taxa.current(5)
combo_taxa.pack()

frame_Label = Label(frame, text="Choose the desired percentage", bg="white")
frame_Label.pack(padx=20, pady=20)

# list of percentage
options_perc = [0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
combo_perc = ttk.Combobox(frame, value=options_perc, state="readonly")
combo_perc.current(0)
combo_perc.pack()
# print(type(combo_perc.get()))
# print(combo_taxa.get())


# upload taxonomy file button
upload_tax = Button(root, text="Upload tax file", height=10, width=20,
                    activebackground="gray", command=upload_tax_f)
upload_tax.grid(row=0, column=0, sticky="W")
upload_tax.bind("<Enter>", button_hover)
upload_tax.bind("<Leave>", button_hover_leave)

# upload otu file button
upload_otu = Button(root, text="Upload otu file", height=10, width=20,
                    activebackground="gray", command=upload_otu_f)
upload_otu.grid(row=1, column=0, sticky="W")
upload_otu.bind("<Enter>", button_hover_1)
upload_otu.bind("<Leave>", button_hover_leave_1)

pars_button = Button(root, text="concatenate both files", height=3, width=20,
                     activebackground="gray", command=concat_tables)
pars_button.grid(row=3, column=0, sticky="W")
pars_button.bind("<Enter>", button_hover_concat)
pars_button.bind("<Leave>", button_hover_concat_leave)

bar_button = Button(root, text="Bar plot", height=3, width=9,
                    activebackground="gray", command=create_bar)
bar_button.grid(row=3, column=1)
bar_button.bind("<Enter>", button_barplot_enter)
bar_button.bind("<Leave>", button_barplot_leave)

heatmap_button = Button(root, text="Heat map", height=3, width=9,
                        activebackground="gray", command=create_heat)

heatmap_button.grid(row=3, column=1, sticky="W")
heatmap_button.bind("<Enter>", button_heat_map_enter)
heatmap_button.bind("<Leave>", button_heat_map_leave)


pca_button = Button(root, text="PCA", height=3, width=9,
                    activebackground="gray", command=open_window)

pca_button.grid(row=3, column=1, sticky="E")
pca_button.bind("<Enter>", pca_enter)
pca_button.bind("<Leave>", pca_leave)


def calculate():
    global  dnew
    dnew = NGS_app.keep_more_than(combo_taxa.get(), float(combo_perc.get()), data)
    inf = NGS_app.info(dnew, combo_taxa.get())
    label_info = Label(frame1, text=inf, bg="white")
    frame1_label.destroy()
    label_info.pack(padx=20, pady=20)


calc_button = Button(root, text="Calculate", height=2, width=11,
                     activebackground="gray", command=calculate)
calc_button.grid(row=1, column=1, sticky="s")

calc_button.bind("<Enter>", button_hover_calc)
calc_button.bind("<Leave>", button_hover_calc_leave)

#status Labels
status_label = Label(root, text="", bd=1,
                     relief="sunken")
status_label.grid(row=5, column=2)

root.mainloop()

