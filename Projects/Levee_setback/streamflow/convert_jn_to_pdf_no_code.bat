rem cd "C:\Users\ajcalder\Box\Classwork\ECI_240\Eutrophication HW\"
call activate geosp
jupyter nbconvert --to=pdf Key_output_review_html.ipynb --TemplateExporter.exclude_input=True

PAUSE