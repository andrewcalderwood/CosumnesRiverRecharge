rem cd "C:\Users\ajcalder\Box\Classwork\ECI_240\Eutrophication HW\"
call activate geosp
jupyter nbconvert --to=html Key_output_review_html.ipynb --TemplateExporter.exclude_input=True --TemplateExporter.exclude_output_prompt=True

PAUSE