#!/usr/bin/env python3 

import yaml

# 定义一个字典
data = {
    'name': 'John Doe',
    'age': 30,
    'children': ['Jane Doe', 'Doe Jr.'],
    'married': True,
    'salary': None
}

# 将字典写入 YAML 文件
with open('output.yaml', 'w') as file:
    yaml.dump(data, file, default_flow_style=False, allow_unicode=True)

from reportlab.lib.pagesizes import letter
from reportlab.pdfgen import canvas
from reportlab.lib.utils import ImageReader

def create_pdf(image_paths, descriptions, output_path):
    c = canvas.Canvas(output_path, pagesize=letter)
    width, height = letter  # Size of the page
    
    text_height = 750
    for img_path, desc in zip(image_paths, descriptions):
        # Add image
        img = ImageReader(img_path)
        img_width, img_height = img.getSize()
        aspect = img_width / img_height
        c.drawImage(img, 100, 500, width=400, height=400 / aspect)
        
        # Add description
        c.drawString(100, text_height, desc)
        text_height -= 50 + 400 / aspect  # Adjust for next image and text
        
        # Move to next page if necessary
        if text_height < 100:
            c.showPage()
            text_height = 750

    c.save()

# Image paths and their descriptions
images = [
    'path_to_image1.png',
    'path_to_image2.png',
    'path_to_image3.png'
]
descriptions = [
    'Description of the first image.',
    'Description of the second image.',
    'Description of the third image.'
]

# Output path for the PDF
output_pdf_path = 'output.pdf'

# Create the PDF
create_pdf(images, descriptions, output_pdf_path)
