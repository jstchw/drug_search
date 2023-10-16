import csv
import io
import openai
import sys
import json

from flask import Flask, request, jsonify, send_file, abort
from flask_cors import CORS
from lxml import etree as ET
from rdkit import Chem
from rdkit.Chem import Draw

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
        'name': ((temp := element.find('db:name', ns)) is not None and temp.text) or None,
        'half_life': ((temp := element.find('db:half-life', ns)) is not None and temp.text) or None,
        'classification': ((temp := element.find('db:classification/db:class', ns)) is not None and temp.text) or None,
        'groups': ((temp := element.findall('db:groups/db:group', ns)) is not None and [group.text for group in
                                                                                        temp]) or None,
        'brands': ((temp := element.findall('db:international-brands/db:international-brand', ns)) is not None and [
            brand.find('db:name', ns).text for brand in temp]) or None,
        'iupac': ((temp := element.find('db:calculated-properties/db:property[db:kind="IUPAC Name"]/db:value',
                                        ns)) is not None and temp.text) or None,
        'formula': ((temp := element.find('db:calculated-properties/db:property[db:kind="Molecular Formula"]/db:value',
                                          ns)) is not None and temp.text) or None,
        'indication': ((temp := element.find('db:indication', ns)) is not None and temp.text) or None,
        'description': ((temp := element.find('db:state', ns)) is not None and temp.text) or None,
        'product': product_name,
    }


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

    if search_type == 'patient.drug.openfda.generic_name':
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
    if not results:
        abort(404, description="No results found")
    return jsonify(results)


@app.route('/api/get_suggestions', methods=['GET'])
def get_suggestions():
    with open('data/suggestion_db.csv', 'r') as file:
        reader = csv.DictReader(file)
        return jsonify([row['Common name'] for row in reader])


@app.route('/api/get_summary', methods=['POST'])
def get_summary():
    openai.api_key_path = 'keys/openai.txt'
    data = request.get_json()
    instruction = 'You are an AI trained to analyze and summarize data from charts. ' \
                  'Your task is to derive new conclusions and facts from the given data, ' \
                  'not just rephrase the information. You should provide a concise summary of the key insights ' \
                  'focusing on the most valuable information and comparing the data. Your response should be limited ' \
                  'to 300-400 characters and should not include any information other than the summary. ' \
                  'Here is the data: '

    chat_completion = openai.ChatCompletion.create(model="gpt-3.5-turbo",
                                                   temperature=1,
                                                   messages=[{"role": "system", "content": instruction},
                                                             {"role": "user", "content": json.dumps(data)}])
    return jsonify(chat_completion)


@app.route('/api/get_causes', methods=['GET'])
def get_causes():
    with open('data/causal/root.csv', 'r') as file:
        reader = csv.DictReader(file)
        filtered_rows = [
            {"feature": row["Feature"], "value": row["value"], "z_score": float(row["z score"])}
            for row in reader
        ]
        return jsonify(filtered_rows)



def on_startup():
    global root
    tree = ET.parse('data/drugs_cleaned.xml')
    root = tree.getroot()


on_startup()

if __name__ == '__main__':
    app.run(debug=True)
