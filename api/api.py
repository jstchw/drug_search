from flask import Flask, request, jsonify
from flask_cors import CORS
from lxml import etree as ET
import sys
import csv
import re

app = Flask(__name__)
CORS(app)

def process_element(element, drug_id):
    ns = {'db': 'http://www.drugbank.ca'}
    drug_tag = element.find('db:drugbank-id', ns)
    if drug_tag is not None and drug_tag.text == drug_id:
        drug_name = element.find('db:name', ns).text
        half_life = element.find('db:half-life', ns)
        classification = element.find('db:classification/db:class', ns)
        products = [prod.find('db:name', ns).text for prod in element.findall('db:products/db:product', ns)]

        return {
            'name': drug_name,
            'half_life': half_life.text if half_life is not None else None,
            'classification': classification.text if classification is not None else None,
            'products': list(set(products))
        }
    return None



def search_csv(column_name, search_term):
    with open('data/db_vocab.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row[column_name].lower() == search_term.lower():
                return row['DrugBank ID']
        return None





@app.route('/api/get_info', methods=['GET'])
def get_info():
    drug_name = request.args.get('drug_name')
    drug_id = search_csv('Common name', drug_name)
    context = ET.iterparse('data/drugs.xml', events=('end', ))
    results = []

    for event, element in context:
        if element.tag == '{http://www.drugbank.ca}drug' and element.getparent().tag == '{http://www.drugbank.ca}drugbank':
            drug = process_element(element, drug_id)
            if drug:
                results.append(drug)
                break
            element.clear()
    return jsonify(results)


if __name__ == '__main__':
    app.run(debug=True)
