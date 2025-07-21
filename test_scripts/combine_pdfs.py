
from pypdf import PdfMerger
import sys

pdfs = []
for i in sys.argv[2:]:
	pdfs.append(i)
merger = PdfMerger()

for pdf in pdfs:
	merger.append(pdf)

merger.write(sys.argv[1])
merger.close()
"""
from pypdf import PdfReader, PdfWriter
import sys

# Function to extract specific pages from a PDF
def extract_pages(pdf_path, pages):
    reader = PdfReader(pdf_path)
    writer = PdfWriter()
    
    for page_num in pages:
        writer.add_page(reader.pages[page_num])
    
    return writer

def combine_pdfs(pdfs_and_pages, output_path):
    combined_writer = PdfWriter()
    
    for pdf, pages in pdfs_and_pages:
        writer = extract_pages(pdf, pages)
        
        for page in writer.pages:
            combined_writer.add_page(page)
    
    with open(output_path, "wb") as f_out:
        combined_writer.write(f_out)

# List of PDFs and the pages to extract (0-indexed)
pdfs_and_pages = [
    (sys.argv[1], [0, 1, 2]),
    (sys.argv[2], [4]),
    (sys.argv[1], [4, 3]),
    (sys.argv[2], [5]),
    (sys.argv[1], [5, 6]),
    (sys.argv[3], [0, 1, 2]),
    (sys.argv[4], [4]),
    (sys.argv[3], [4, 3]),
    (sys.argv[4], [5]),
    (sys.argv[3], [5, 6])
]

combine_pdfs(pdfs_and_pages, "combined_output.pdf")
"""
