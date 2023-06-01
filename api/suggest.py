import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup


def process_element(element):
    soup = BeautifulSoup(ET.tostring(element), 'xml')
    # Perform your search and data extraction here using Beautiful Soup
    # For example, to find a specific tag:
    tag = soup.find('name')
    if tag:
        print(tag.text)


def parse_large_xml(file_path):
    context = ET.iterparse(file_path, events=('end',))
    for event, element in context:
        process_element(element)
        element.clear()


file_path = 'data/drugs.xml'
parse_large_xml(file_path)
