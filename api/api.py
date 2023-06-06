from flask import Flask, request, jsonify, send_file
from flask_cors import CORS
from lxml import etree as ET
import sys
import csv
from rdkit import Chem
from rdkit.Chem import Draw
import io

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
            'groups': [group.text for group in element.findall('db:groups/db:group', ns)],
            # Find brand XPath: /drugbank/drug/international-brands/international-brand/name
            'brands': [brand.find('db:name', ns).text for brand in element.findall('db:international-brands/db:international-brand', ns)],
            'iupac': element.find('db:calculated-properties/db:property[db:kind="IUPAC Name"]/db:value', ns).text,
            'formula': element.find('db:calculated-properties/db:property[db:kind="Molecular Formula"]/db:value', ns).text,
            'indication': element.find('db:indication', ns).text,
            'description': element.find('db:state', ns).text,

        }
    return None



def search_csv(column_name, search_term):
    with open('data/suggestion_db.csv', 'r') as file:
        reader = csv.DictReader(file)
        for row in reader:
            if row[column_name].lower() == search_term.lower():
                return row['DrugBank ID']
        return None, 404


@app.route('/api/get_molecule', methods=['GET'])
def get_molecule():
    drug_name = request.args.get('drug_name')
    drug_id = search_csv('Common name', drug_name)
    suppl = Chem.SDMolSupplier('data/molecules.sdf')

    for mol in suppl:
        if mol is not None and mol.GetProp('DRUGBANK_ID') == drug_id:
            img = Draw.MolToImage(mol)
            img_bytes = io.BytesIO()
            img.save(img_bytes, format='PNG')
            img_bytes = img_bytes.getvalue()
            return send_file(io.BytesIO(img_bytes), mimetype='image/png')
    return "No molecule found", 404




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

@app.route('/api/get_suggestions', methods=['GET'])
def get_suggestions():
    with open('data/suggestion_db.csv', 'r') as file:
        reader = csv.DictReader(file)
        return jsonify([row['Common name'] for row in reader])


if __name__ == '__main__':
    app.run(debug=True)
