from flask import Flask, request, jsonify
from flask_cors import CORS
import xml.etree.ElementTree as ET
from bs4 import BeautifulSoup

app = Flask(__name__)
CORS(app)

def process_element(element, drug_name):
    soup = BeautifulSoup(ET.tostring(element), 'xml')
    drug_tag = soup.find('drug/name')
    if drug_tag and drug_name.lower() in drug_tag.text.lower():
        return drug_tag.text

@app.route('/api/get_info', methods=['GET'])
def get_info():
    drug_name = request.args.get('drug_name')
    context = ET.iterparse('api/data/drugs.xml', events=('end', ))
    results = []
    for event, element in context:
        drug = process_element(element, drug_name)
        if drug:
            drugs.append(drug)
        element.clear()
    return jsonify(drugs)

if __name__ == '__main__':
    app.run(debug=True)