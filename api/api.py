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


def process_element(element, drug_id=None, drug_name=None, drug_brand=None):
    ns = {'db': 'http://www.drugbank.ca'}

    drug_tag = None

    if drug_id is not None:
        drug_tag = element.find('db:drugbank-id[@primary="true"]', ns)
        if drug_tag is not None and drug_tag.text != drug_id:
            return None

    product_name = None
    if drug_brand is not None:
        products = element.findall('db:products/db:product', ns)
        for product in products:
            if product.find('db:name', ns).text.lower() == drug_brand.lower():
                product_name = product.find('db:name', ns).text
                break
    if product_name is None:
            product_name = drug_tag.text

    return {
        'name': element.find('db:name', ns).text,
        'half_life': element.find('db:half-life', ns).text,
        'classification': element.find('db:classification/db:class', ns).text,
        'groups': [group.text for group in element.findall('db:groups/db:group', ns)],
        # Find brand XPath: /drugbank/drug/international-brands/international-brand/name
        'brands': [brand.find('db:name', ns).text for brand in element.findall('db:international-brands/db:international-brand', ns)],
        'iupac': element.find('db:calculated-properties/db:property[db:kind="IUPAC Name"]/db:value', ns).text,
        'formula': element.find('db:calculated-properties/db:property[db:kind="Molecular Formula"]/db:value', ns).text,
        'indication': element.find('db:indication', ns).text,
        'description': element.find('db:state', ns).text,
        'product': product_name,
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
    search_type = request.args.get('search_type')

    results = []

    if search_type == 'patient.drug.activesubstance.activesubstancename':
        drug_id = search_csv('Common name', drug_name)
        for element in root.findall('{http://www.drugbank.ca}drug'):
                if element.tag == '{http://www.drugbank.ca}drug' and element.getparent().tag == '{http://www.drugbank.ca}drugbank':
                     drug = process_element(element=element, drug_id=drug_id)
                     if drug:
                         results.append(drug)
                         break
    elif search_type == 'patient.drug.openfda.brand_name':
        for element in root.findall('{http://www.drugbank.ca}drug'):
                    if element.tag == '{http://www.drugbank.ca}drug' and element.getparent().tag == '{http://www.drugbank.ca}drugbank':
                        products_element = element.find('{http://www.drugbank.ca}products')
                        if products_element is not None:
                            for product_element in products_element.findall('{http://www.drugbank.ca}product'):
                                name_element = product_element.find('{http://www.drugbank.ca}name')
                                if name_element is not None and name_element.text.lower() == drug_name.lower():
                                    drug = process_element(element=element, drug_name=name_element.text,
                                                           drug_brand=drug_name.lower())
                                    if drug:
                                        results.append(drug)
                                        break
    return jsonify(results)

@app.route('/api/get_suggestions', methods=['GET'])
def get_suggestions():
    with open('data/suggestion_db.csv', 'r') as file:
        reader = csv.DictReader(file)
        return jsonify([row['Common name'] for row in reader])


def on_startup():
    global root
    tree = ET.parse('data/drugs.xml')
    root = tree.getroot()


on_startup()

if __name__ == '__main__':
    app.run(debug=True)
