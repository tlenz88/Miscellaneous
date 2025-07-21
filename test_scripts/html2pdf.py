import pdfkit
import sys
import os

pdfkit.from_file(sys.argv[1], f"{os.path.basename(sys.argv[1])[:-5]}.pdf")
