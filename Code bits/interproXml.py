import xml.etree.ElementTree as ET

file = "output0.xml"
url = "{http://www.ebi.ac.uk/interpro/resources/schemas/interproscan5}"

tree = ET.parse(file)
root = tree.getroot()

#print(root.attrib)

for element in root:
    print(element.tag, '\n',  element.attrib)

for thing in tree.iter(url + 'sequence'):
    print(thing)
