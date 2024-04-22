Ti_label = tk.Text(inner_frame, height=1, width=20, bg="white", bd=0, font='Helvetica 11 bold')
Ti_label.grid(row=4, column=0)

# Inserting text
Ti_label.insert(tk.END, "Starting Temperature, ")
Ti_label.insert(tk.END, "T", "bold_italic")
Ti_label.insert(tk.END, "0", "subscript_bold_italic")  # Marking "i" as subscript, bold and italic

# Tag configuration for bold and italic text
Ti_label.tag_configure("bold_italic", font=("Helvetica", 11, "bold italic"))

# Tag configuration for subscript, bold and italic text
Ti_label.tag_configure("subscript_bold_italic", font=("Helvetica", 11, "bold italic"), offset=-4)

# Disable editing
Ti_label.configure(state="disabled")