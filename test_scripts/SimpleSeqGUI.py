#!/usr/bin/env python3

import tkinter as tk
from tkinter import filedialog

class DirectorySelector:
    def __init__(self, root):
        self.root = root
        self.root.title("Directory Selection App")
        self.root.geometry("640x480")
        self.root.resizable(False, False)
        self.selected_directory = tk.StringVar()
        self.label = tk.Label(root, text="Select input directory:")
        self.label.pack(ipadx=10, ipady=15, side=tk.LEFT, fill=tk.BOTH, expand=False)
        self.select_button = tk.Button(self.root, text="Browse", command=self.select_directory)
        self.select_button.pack(ipadx=30, ipady=10, side=tk.LEFT, fill=tk.BOTH, expand=False)
        self.directory_label = tk.Label(self.root, textvariable=self.selected_directory, background="white", relief="sunken", width=40)
        self.directory_label.pack(ipadx=60, ipady=15, side=tk.LEFT, fill=tk.BOTH, expand=False)

    def select_directory(self):
        directory_path = filedialog.askdirectory()
        if directory_path:
            self.selected_directory.set(directory_path)

class DirectorySelector:
    def __init__(self, root):
        self.root = root
        self.root.title("Directory Selection App")
        self.root.geometry("640x480")
        self.root.resizable(False, False)

        self.selected_directory = tk.StringVar()

        self.label = tk.Label(root, text="Select input directory:")
        self.label.pack(ipadx=10, ipady=15, side=tk.LEFT, fill=tk.BOTH, expand=False)

        self.select_button = tk.Button(self.root, text="Browse", command=self.select_directory)
        self.select_button.pack(ipadx=30, ipady=10, side=tk.LEFT, fill=tk.BOTH, expand=False)

        self.directory_label = tk.Label(self.root, textvariable=self.selected_directory, background="white", relief="sunken", width=40)
        self.directory_label.pack(ipadx=60, ipady=15, side=tk.LEFT, fill=tk.BOTH, expand=False)

    def select_directory(self):
        directory_path = filedialog.askdirectory()
        if directory_path:
            self.selected_directory.set(directory_path)

if __name__ == "__main__":
    root = tk.Tk()
    app = DirectorySelector(root)
    root.mainloop()
